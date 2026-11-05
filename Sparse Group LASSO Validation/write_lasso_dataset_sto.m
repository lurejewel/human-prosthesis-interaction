function write_lasso_dataset_sto(outFile, ds)
%WRITE_LASSO_DATASET_STO 写出与 SCONE UN.sto 同布局的 LASSO-ready .sto。
%
% 为什么要“模仿 UN.sto”
% ---------------------
% 旧 LASSO 流程（Adaptive_LASSO_StanceSwing.m + extractStoMuscleFeatures.m）
% 已经围绕 SCONE 的列名写好：
%
%   <muscle>_r.activation
%   <muscle>_r.excitation
%   <muscle>_r.fiber_length_norm
%   <muscle>_r.mtu_force_norm
%   pelvis_tilt / pelvis_tilt_u
%   hip_flexion_r / _u
%   knee_extension_r / _u
%   ankle_dorsiflexion_r / _u
%   leg1_r.state        (0/1 stance, 2/3/4 swing)
%
% 与其修改 reader，不如让 Moco 接口输出同样的列名。这样：
%   * 旧 UN.sto 路径和新 Moco 路径共用同一套特征提取/拟合代码；
%   * 可以直接用旧脚本做对照实验；
%   * 以后换模型/换数据时，LASSO 侧不需要改列名。
%
% 本函数还会写左腿列和 leg0_l.state，方便以后做双侧分析；当前 LASSO 只用
% 右腿（与旧流程的双侧对称假设一致）。
%
% 输入 ds 的字段
% --------------
%   ds.time        nT×1
%   ds.muscleNames 1×nMus，顺序为 right 9 + left 9（或至少与矩阵列一致）
%   ds.activation  nT×nMus  清洗后的 activation
%   ds.excitation  nT×nMus  清洗后的 excitation target
%   ds.lCEN        nT×nMus  归一化纤维长度
%   ds.fATN        nT×nMus  归一化 MTU 力
%   ds.q, ds.qdot, ds.qNames  坐标/速度（用于写关节列）
%   ds.phase.stateR/stateL    nT×1，0..4
%
% 输出
% ----
%   标准 OpenSim .sto：
%     LASSO dataset generated from MocoInverse
%     version=1
%     nRows=...
%     nColumns=...
%     inDegrees=no
%     endheader
%     time<TAB>...
%
% 注意：
% * nColumns 包含 time 列；OpenSim TimeSeriesTable 接受这种写法；
% * 所有列必须是有限数值；NaN 会让 OpenSim reader 报错，接口在调用前
%   已经用 validMask 做过检查，但这里不做静默替换；
% * 列顺序中肌肉块是 interleaved（activation/excitation/fiber/force），
%   extractStoMuscleFeatures 会用 endsWith 过滤，不依赖连续块。

% -------------------------------------------------------------------------
% 0. 最小一致性检查
% -------------------------------------------------------------------------
% 只检查 activation 的尺寸，因为它是 LASSO 的 target 来源；其他矩阵如果
% 尺寸不对，后面的列拼接/写入会自然报错（或产生明显错误的列数）。
if nargin < 2
    error('write_lasso_dataset_sto:badInput', 'outFile and ds are required.');
end
nT = numel(ds.time);
nMus = numel(ds.muscleNames);
if size(ds.activation, 1) ~= nT || size(ds.activation, 2) ~= nMus
    error('write_lasso_dataset_sto:badActivation', ...
        'ds.activation must be nT x nMus.');
end

% labels/data 用“边加边扩”的方式构建：总列数约 88，规模很小，可读性
% 比预分配更重要。
labels = {'time'};
data = ds.time(:);

% -------------------------------------------------------------------------
% 1. 肌肉列
% -------------------------------------------------------------------------
% 每块肌肉连续写 4 列，列名后缀严格使用 OpenSim/SCONE 的命名。
% 注意 fiber_length_norm / mtu_force_norm 是 extractStoMuscleFeatures
% 用 endsWith 识别的关键列，不能改名。
for i = 1:nMus
    nm = ds.muscleNames{i};
    labels{end+1} = [nm '.activation'];          %#ok<AGROW>
    data(:, end+1) = ds.activation(:, i);        %#ok<AGROW>
    labels{end+1} = [nm '.excitation'];          %#ok<AGROW>
    data(:, end+1) = ds.excitation(:, i);        %#ok<AGROW>
    labels{end+1} = [nm '.fiber_length_norm'];   %#ok<AGROW>
    data(:, end+1) = ds.lCEN(:, i);              %#ok<AGROW>
    labels{end+1} = [nm '.mtu_force_norm'];      %#ok<AGROW>
    data(:, end+1) = ds.fATN(:, i);              %#ok<AGROW>
end

% -------------------------------------------------------------------------
% 2. 矢状面关节列（extractStoMuscleFeatures 真正使用的特征）
% -------------------------------------------------------------------------
% 旧 reader 会按 side 查找：
%   pelvis_tilt / pelvis_tilt_u
%   hip_flexion_r / _u
%   knee_extension_r / _u
%   ankle_dorsiflexion_r / _u
% 这里把左右都写；当前 LASSO 只用右腿。get_coord 在坐标缺失时返回 0，
% 不会因为某个模型没有 pelvis_tx 而中断整张表。
jointBases = {'pelvis_tilt', ...
              'hip_flexion_r', 'knee_extension_r', 'ankle_dorsiflexion_r', ...
              'hip_flexion_l', 'knee_extension_l', 'ankle_dorsiflexion_l'};
for b = 1:numel(jointBases)
    nm = jointBases{b};
    [v, vd] = get_coord(ds, nm);
    labels{end+1} = nm;    data(:, end+1) = v;   %#ok<AGROW>
    labels{end+1} = [nm '_u']; data(:, end+1) = vd; %#ok<AGROW>
end

% -------------------------------------------------------------------------
% 3. 骨盆平移（LASSO 不用，但 QC/画图/排查时有用）
% -------------------------------------------------------------------------
% 只有在 ds.qNames 里确实存在时才写；否则 get_coord 会返回全 0，写一列
% 假 0 会误导 QC。当前 human0918/SQR 都有 pelvis_tx/ty，所以不影响主流程。
for nm = {'pelvis_tx', 'pelvis_ty'}
    if isfield(ds, 'qNames') && any(strcmp(ds.qNames, nm{1}))
        [v, vd] = get_coord(ds, nm{1});
        labels{end+1} = nm{1};       data(:, end+1) = v;  %#ok<AGROW>
        labels{end+1} = [nm{1} '_u']; data(:, end+1) = vd; %#ok<AGROW>
    end
end

% -------------------------------------------------------------------------
% 4. 步态相位标签
% -------------------------------------------------------------------------
% 旧 reader 的约定：
%   right side -> 'leg1_r.state'
%   left  side -> 'leg0_l.state'
% 值必须为 0..4；extractStoMuscleFeatures 会检查这个范围。
% 我们的 0/1 与 2/3/4 子划分是启发式，但 stance/swing 的二分来自 GRF。
labels{end+1} = 'leg1_r.state'; data(:, end+1) = ds.phase.stateR(:); %#ok<AGROW>
labels{end+1} = 'leg0_l.state'; data(:, end+1) = ds.phase.stateL(:); %#ok<AGROW>

% 内部一致性检查：labels 和 data 必须同步增长。如果这里报错，说明上面某个
% 分支只加了 label 没加 data（或反过来），属于代码 bug，不是数据问题。
nCol = numel(labels);
if size(data, 2) ~= nCol
    error('write_lasso_dataset_sto:internal', ...
        'Column bookkeeping error: %d labels vs %d data columns.', ...
        nCol, size(data, 2));
end

% -------------------------------------------------------------------------
% 5. 写文件
% -------------------------------------------------------------------------
% 头部只使用标准 OpenSim 字段（version/nRows/nColumns/inDegrees/endheader），
% 第一行是自由文本名字。OpenSim TimeSeriesTable 接受这种布局（已在
% OpenSim 4.6 Python 绑定上验证过）。
fid = fopen(outFile, 'w');
if fid < 0
    error('write_lasso_dataset_sto:cannotOpen', 'Cannot write %s.', outFile);
end
cleaner = onCleanup(@() fclose(fid));

fprintf(fid, 'LASSO dataset generated from MocoInverse\n');
fprintf(fid, 'version=1\n');
fprintf(fid, 'nRows=%d\n', nT);
fprintf(fid, 'nColumns=%d\n', nCol);
fprintf(fid, 'inDegrees=no\n');
fprintf(fid, 'endheader\n');
fprintf(fid, '%s\n', strjoin(labels, '\t'));

fmt = ['%.8f', repmat('\t%.8f', 1, nCol - 1), '\n'];
for k = 1:nT
    fprintf(fid, fmt, data(k, :));
end

fprintf('[write] %s  (%d rows x %d columns)\n', outFile, nT, nCol);
end

% ========================================================================
% 局部函数：按坐标名取 value / speed
% ========================================================================
% 返回值约定：
%   * 坐标存在 -> 返回 ds.q / ds.qdot 的对应列；
%   * 坐标不存在 -> 返回全 0（不是 NaN），这样写出的 .sto 仍是数值表；
%   * ds.qdot 缺失或列数不足 -> speed 返回全 0。
% 对 LASSO 必需的 7 个矢状坐标，接口保证它们存在；pelvis_tx/ty 是可选列，
% 调用方会先检查 qNames 再决定是否写。
function [v, vd] = get_coord(ds, name)
%GET_COORD  Return a coordinate and its speed, or zeros if absent.
v  = nan(size(ds.time));
vd = nan(size(ds.time));
if ~isfield(ds, 'qNames') || isempty(ds.qNames)
    v(:) = 0;
    vd(:) = 0;
    return;
end
j = find(strcmp(ds.qNames, name), 1);
if isempty(j)
    v(:) = 0;
    vd(:) = 0;
    return;
end
v = ds.q(:, j);
if isfield(ds, 'qdot') && ~isempty(ds.qdot) && size(ds.qdot, 2) >= j
    vd = ds.qdot(:, j);
else
    vd(:) = 0;
end
end
