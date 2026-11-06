function S = read_moco_inverse_solution(filePath, muscleNames)
%READ_MOCO_INVERSE_SOLUTION 从 MocoInverse 的 .sto 中抽取肌肉/坐标列。
%
% 这个函数解决什么问题
% --------------------
% 本项目的 MocoInverse 输出有几种“形态”，但下游接口必须用统一的 S 结构：
%
%   (a) 完整解（未加 --sagittal-cost，或 Moco 自己写的完整 trajectory）：
%         /forceset/<muscle>/activation   肌肉激活 state
%         /forceset/<muscle>              肌肉兴奋 control（裸列名）
%         /jointset/<joint>/<coord>/value 坐标
%         /jointset/<joint>/<coord>/speed 坐标速度
%
%   (b) 只有 activation 的解（本仓库初步 --sagittal-cost 的常见形态）：
%         只有 18 个 /forceset/<muscle>/activation，没有坐标、没有 excitation；
%         坐标需要另读 ik_radians.mot。
%         注意：本 parser 仍能读出这种文件（缺失列保持 NaN 作为内部
%         sentinel），但 moco_inverse_to_lasso_dataset.m 会因为没有
%         excitation 而报 MissingExcitation；正式流程同时要求完整的
%         excitation 和 activation。
%
%   (c) 只有 excitation 的辅助文件（moco_inverse_excitations_*.sto）：
%         只有 /forceset/<muscle>。
%
% 本函数只做“列名解析 + 矩阵组装”，不做单位转换、滤波、相位判断；
% 这些属于 moco_inverse_to_lasso_dataset.m 的职责。
%
% 列名匹配规则（按顺序，先匹配 activation，再匹配 excitation，最后坐标）
% --------------------------------------------------------------------
%   activation : ^/?forceset/([^/]+)/activation$
%   excitation : ^/?forceset/([^/]+)$
%                但必须排除 residual_* / reserve_*，否则会把骨盆残差/储备
%                作动器误当成肌肉兴奋。
%   坐标       : ^/?jointset/[^/]+/([^/]+)/(value|speed)$
%
% 输入
% ----
%   filePath    : .sto/.mot 路径；
%   muscleNames : 1×nMus cellstr，按 OpenSim 全名给出，例如
%                 {'hamstrings_r', ..., 'tib_ant_l'}。
%                 返回矩阵的列严格按这个顺序；文件里缺失的肌肉保持 NaN
%                 作为内部 sentinel。接口会立即检查并报错，不会把 NaN
%                 带入 LASSO 数据集。
%
% 输出结构 S
% ----------
%   S.time          nT×1    时间 [s]
%   S.activation    nT×nMus 肌肉激活（缺失为 NaN）
%   S.excitation    nT×nMus 肌肉兴奋（缺失为 NaN）
%   S.hasActivation 逻辑，是否存在任何有限 activation
%   S.hasExcitation 逻辑，是否存在任何有限 excitation
%   S.coordNames    1×nQ   坐标名，顺序与 S.q 的列一致
%   S.q             nT×nQ  坐标值（文件原单位；Moco 内部是 rad / m）
%   S.qdot          nT×nQ  坐标速度（缺失列保持 NaN）
%   S.hasCoordinates 逻辑，是否存在 /jointset/.../value
%   S.hasCoordinateVel 逻辑，文件里是否出现过 /speed 列
%   S.labels        原始列名（1×nC cellstr）
%   S.file          文件路径
%   S.muscleNames   调用者传入的肌肉顺序（便于下游核对）
%
% 重要约定
% --------
% * 本函数不根据 inDegrees 做 deg2rad：Moco 解里的 /jointset/.../value
%   是 SimTK 内部单位（rad / m），不是度；IK 的 deg->rad 在接口里处理。
% * 如果文件里只有 activation 或只有 excitation，另一个矩阵保持 NaN
%   作为 sentinel；S.hasActivation / S.hasExcitation 让调用者判断。
%   注意：moco_inverse_to_lasso_dataset.m 要求两者完整，缺失会直接报错。
% * 坐标 /value 和 /speed 可能在不同位置出现，甚至只出现一部分；
%   最后统一按 coordNames 重排 qdot，保证与 q 的列一一对应。

% -------------------------------------------------------------------------
% 0. 参数校验与规范化
% -------------------------------------------------------------------------
if nargin < 1 || isempty(filePath)
    error('read_moco_inverse_solution:badInput', 'filePath is required.');
end
if nargin < 2 || isempty(muscleNames)
    error('read_moco_inverse_solution:badInput', 'muscleNames is required.');
end
if ischar(muscleNames) || isstring(muscleNames)
    muscleNames = cellstr(muscleNames);
end
% 统一成 1×nMus 行向量，后面 find(strcmp(...)) 的返回值是标量列号。
muscleNames = muscleNames(:).';
nMus = numel(muscleNames);

% 纯 MATLAB 读取，不依赖 OpenSim；返回 labels/data/time。
T = read_sto_mot_table(filePath);
nT = numel(T.time);

% -------------------------------------------------------------------------
% 1. 预分配输出
% -------------------------------------------------------------------------
% 肌肉矩阵先填 NaN：文件里没有的肌肉保持 NaN，下游可据此判断“缺失”而不是
% 误以为肌肉兴奋为 0。坐标矩阵先留空，遇到 /jointset 列时再动态追加。
S = struct();
S.file          = T.file;
S.time          = T.time(:);
S.labels        = T.labels;
S.muscleNames   = muscleNames;
S.activation    = nan(nT, nMus);
S.excitation    = nan(nT, nMus);
S.coordNames    = {};
S.q             = [];
S.qdot          = [];

% qMap/qdMap: 坐标名 -> 在 S.q / 原始 S.qdot 中的列号。
% 用 containers.Map 而不是 find(strcmp) 是为了 O(1) 查找；坐标数量只有 19，
% 但写成通用形式便于以后换 full-body 模型。
qMap   = containers.Map('KeyType', 'char', 'ValueType', 'int32');
qdMap  = containers.Map('KeyType', 'char', 'ValueType', 'int32');

% -------------------------------------------------------------------------
% 2. 逐列解析
% -------------------------------------------------------------------------
% 从第 2 列开始（第 1 列固定是 time）。每个分支用 continue 跳出，避免一个
% 列名同时被多个正则匹配（例如 /forceset/x/activation 也会匹配 excitation
% 的宽松模式，所以 activation 必须先判断）。
for k = 2:numel(T.labels)
    lab = T.labels{k};
    vals = T.data(:, k);

    % ---- 2.1 /forceset/<muscle>/activation ----------------------------
    % 这是肌肉激活 state。Moco 的完整解和 activation-only 解都有这种列。
    tok = regexp(lab, '^/?forceset/([^/]+)/activation$', 'tokens', 'once');
    if ~isempty(tok)
        j = find(strcmp(muscleNames, tok{1}), 1);
        if ~isempty(j)
            S.activation(:, j) = vals;
        end
        continue;
    end

    % ---- 2.2 /forceset/<muscle>（肌肉兴奋 control） -------------------
    % 裸列名就是 excitation control。必须排除 residual_/reserve_：
    % 它们也是 /forceset/<name> 形式，但不是肌肉，不能进入 LASSO 目标。
    tok = regexp(lab, '^/?forceset/([^/]+)$', 'tokens', 'once');
    if ~isempty(tok)
        name = tok{1};
        if startsWith(name, 'residual_') || startsWith(name, 'reserve_')
            continue;
        end
        j = find(strcmp(muscleNames, name), 1);
        if ~isempty(j)
            S.excitation(:, j) = vals;
        end
        continue;
    end

    % ---- 2.3 /jointset/<joint>/<coordinate>/<value|speed> ------------
    % Moco 的坐标列路径固定为 5 段；第 4 段是坐标名，第 5 段是 value/speed。
    % value 追加到 S.q，speed 追加到 S.qdot；两者最后按名字对齐。
    tok = regexp(lab, '^/?jointset/[^/]+/([^/]+)/(value|speed)$', 'tokens', 'once');
    if ~isempty(tok)
        cname = tok{1};
        kind  = tok{2};
        if strcmpi(kind, 'value')
            if ~isKey(qMap, cname)
                qMap(cname) = size(S.q, 2) + 1;
                % 动态追加一列；nT 可能为 0 的极端情况不会出现在 Moco 输出。
                S.q(:, end+1) = nan(nT, 1); %#ok<AGROW>
                S.coordNames{end+1} = cname; %#ok<AGROW>
            end
            S.q(:, qMap(cname)) = vals;
        else
            if ~isKey(qdMap, cname)
                qdMap(cname) = size(S.qdot, 2) + 1;
                S.qdot(:, end+1) = nan(nT, 1); %#ok<AGROW>
            end
            S.qdot(:, qdMap(cname)) = vals;
        end
    end
end

% -------------------------------------------------------------------------
% 3. 状态标志
% -------------------------------------------------------------------------
% 用 any(isfinite(...)) 而不是 isempty：矩阵可能被预分配成 NaN，但文件里
% 一列都没有；只有出现有限值才算“存在”。
S.hasActivation    = any(isfinite(S.activation(:)));
S.hasExcitation    = any(isfinite(S.excitation(:)));
S.hasCoordinates   = ~isempty(S.q);
S.hasCoordinateVel = ~isempty(S.qdot);

% -------------------------------------------------------------------------
% 4. 把 qdot 对齐到 q / coordNames 的顺序
% -------------------------------------------------------------------------
% Moco 通常把所有 /value 写在前面、/speed 写在后面，但列顺序不是契约。
% 这里按名字重建 Qd：缺失的 speed 填 NaN，接口随后会用 q 的数值微分兜底。
if S.hasCoordinates
    Qd = nan(nT, numel(S.coordNames));
    for c = 1:numel(S.coordNames)
        if isKey(qdMap, S.coordNames{c})
            Qd(:, c) = S.qdot(:, qdMap(S.coordNames{c}));
        end
    end
    S.qdot = Qd;
end

% -------------------------------------------------------------------------
% 5. 至少要有一种肌肉列
% -------------------------------------------------------------------------
% 如果既没有 activation 也没有 excitation，这个文件不可能是 MocoInverse
% 肌肉结果（可能只是 IK/GRF），直接报错，避免下游拿到全 NaN 后才困惑。
if ~S.hasActivation && ~S.hasExcitation
    error('read_moco_inverse_solution:noMuscleColumns', ...
        ['No muscle activation or excitation columns found in %s. ' ...
         'Expected "/forceset/<muscle>" or "/forceset/<muscle>/activation".'], ...
        filePath);
end
end
