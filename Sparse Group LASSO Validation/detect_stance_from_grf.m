function P = detect_stance_from_grf(grfFile, time, comX, opts)
%DETECT_STANCE_FROM_GRF 按 cal_gait_phase.m 的规则重建 0..4 步态相位。
%
% 设计原则
% --------
% 这里不再使用“滞回 + 短 run 清理 + 60/40、20/60/20 启发式”那一套。
% 相位标签必须严格遵循预测仿真里的 `functions/core/cal_gait_phase.m`：
%
%   0 Early Stance : GRF >= th 且 calcnCOM >= pelvisCOM
%   1 Late Stance  : GRF >= th 且 calcnCOM <  pelvisCOM 且 opposite GRF <  th
%   2 Liftoff      : GRF >= th 且 calcnCOM <  pelvisCOM 且
%                    (opposite GRF > th 或 calcnCOM + 1 < pelvisCOM)
%   3 Swing        : GRF <  th 且 calcnCOM <  pelvisCOM
%   4 Landing      : GRF <  th 且 calcnCOM >= pelvisCOM
%
% 后续帧的转移规则（与 cal_gait_phase.m 完全一致）：
%
%   0 -> 1 : calcnCOM < pelvisCOM
%   1 -> 2 : opposite GRF >= th 或 calcnCOM + 1 < pelvisCOM
%   2 -> 3 : own GRF < th
%   3 -> 4 : calcnCOM > pelvisCOM
%   4 -> 0 : own GRF >= th
%
% 其中：
%   * `th = stanceThresholdBW * bodyWeight`，默认 `0.02 * BW`，对应
%     forward_simulation.m 中的 `simCache.stanceTh = 0.02`；
%   * `calcnCOM + 1` 里的 1 是**米**，与 cal_gait_phase.m 源码一致；
%   * GRF 是竖直分量（N），不是归一化值；
%   * COM 由 compute_com_x_from_kinematics.m 用模型计算，不是从 GRF 猜的。
%
% 重要：`P.stanceR/P.stanceL` 是“LASSO 的 stance 组”，即 state 0/1；
% `state 2 (Liftoff)` 在 extractStoMuscleFeatures.m 里被归入 swing 组
% （它用 state>=2 作为 swing）。如果想知道“脚是否还接触地面”，看
% `P.contactR/P.contactL`（GRF >= th）。
%
% 输入
% ----
%   grfFile : .mot/.sto GRF 文件
%   time    : 目标时间网格 nT×1
%   comX    : compute_com_x_from_kinematics 的输出：
%             .pelvis, .calcnR, .calcnL（nT×1 全局 X [m]）
%   opts    : 可选
%       .grfXmlFile        匹配的 ExternalLoads XML
%       .bodyWeight        [N]，默认 78*9.81
%       .stanceThresholdBW 默认 0.02（必须与 forward_simulation.m 一致）
%       .rightColumn/.leftColumn 显式列名覆盖
%
% 输出 P
% ------
%   .stateR/.stateL     nT×1，严格 0..4
%   .stanceR/.stanceL   logical，state 0/1（LASSO stance 组）
%   .contactR/.contactL logical，GRF >= th（实际接触，包含 phase 2）
%   .eventsR/.eventsL   phase 0 onset 的样本索引
%   .forceR/.forceL     插值后的竖直 GRF [N]
%   .thresholds         实际使用的 bodyWeight/threshold
%   .comX               传入的 COM 结构
%   .file/.xmlFile/.rightColumn/.leftColumn

% -------------------------------------------------------------------------
% 0. 参数与输入检查
% -------------------------------------------------------------------------
if nargin < 4 || isempty(opts)
    opts = struct();
end
opts = apply_defaults(opts);

if nargin < 1 || isempty(grfFile) || ~isfile(char(grfFile))
    error('detect_stance_from_grf:noFile', ...
        'A valid GRF .mot/.sto file is required to reconstruct gait phases.');
end
if nargin < 3 || isempty(comX)
    error('detect_stance_from_grf:noComX', ...
        ['comX is required: the 0..4 phase rules need pelvis/calcn COM X. ' ...
         'Call compute_com_x_from_kinematics first.']);
end
grfFile = char(grfFile);
time = time(:);
nT = numel(time);
required = {'pelvis','calcnR','calcnL'};
for i = 1:numel(required)
    if ~isfield(comX, required{i}) || numel(comX.(required{i})) ~= nT
        error('detect_stance_from_grf:badComX', ...
            'comX.%s must exist and have %d entries.', required{i}, nT);
    end
end

% -------------------------------------------------------------------------
% 1. 读 GRF 并找到右/左竖直力列
% -------------------------------------------------------------------------
T = read_sto_mot_table(grfFile);

% 如果用户没显式给 XML，就在同目录找同名 .xml（Moco grf_window_*.xml）。
xmlFile = char(opts.grfXmlFile);
if isempty(xmlFile)
    [d, base] = fileparts(grfFile);
    cand = fullfile(d, [base, '.xml']);
    if isfile(cand)
        xmlFile = cand;
    end
end
[rightCol, leftCol] = resolve_force_columns(T.labels, xmlFile, opts);

jR = find(strcmp(T.labels, rightCol), 1);
jL = find(strcmp(T.labels, leftCol), 1);
if isempty(jR) || isempty(jL)
    error('detect_stance_from_grf:columnNotFound', ...
        'Could not find GRF columns "%s" / "%s" in %s.', ...
        rightCol, leftCol, grfFile);
end

% 线性插值到目标网格；接口已经用 GRF 时间范围裁剪过 tOut，外推只在边界。
fR = interp1(T.time, T.data(:, jR), time, 'linear', 'extrap');
fL = interp1(T.time, T.data(:, jL), time, 'linear', 'extrap');
fR(~isfinite(fR)) = 0;
fL(~isfinite(fL)) = 0;

% -------------------------------------------------------------------------
% 2. 严格按 cal_gait_phase.m 计算 0..4
% -------------------------------------------------------------------------
stanceThN = opts.stanceThresholdBW * opts.bodyWeight;
[stateR, stateL] = phase_from_grf_com( ...
    fR, fL, comX.pelvis, comX.calcnR, comX.calcnL, ...
    stanceThN, time);

% LASSO 的 stance 组：state 0/1；phase 2 Liftoff 在旧 reader 里归 swing。
P.stanceR = stateR <= 1;
P.stanceL = stateL <= 1;

% 真实接触（GRF >= th），保留 phase 2 也算 contact。
P.contactR = fR >= stanceThN;
P.contactL = fL >= stanceThN;

% phase 0 onset：当前是 0、前一个是非 0（或第一帧就是 0）。
P.eventsR = find(stateR == 0 & [true; stateR(1:end-1) ~= 0]);
P.eventsL = find(stateL == 0 & [true; stateL(1:end-1) ~= 0]);

P.time     = time;
P.forceR   = fR;
P.forceL   = fL;
P.stateR   = stateR;
P.stateL   = stateL;
P.thresholds = struct( ...
    'bodyWeight', opts.bodyWeight, ...
    'stanceThresholdBW', opts.stanceThresholdBW, ...
    'stanceThresholdN', stanceThN);
P.comX     = comX;
P.file     = grfFile;
P.xmlFile  = xmlFile;
P.rightColumn = rightCol;
P.leftColumn  = leftCol;
end

% ========================================================================
% 局部函数 1：默认参数
% ========================================================================
% stanceThresholdBW 必须与 forward_simulation.m 的 simCache.stanceTh
% 保持一致（当前 0.02）。不要再引入 enter/exit 两个阈值。
function opts = apply_defaults(opts)
defaults = struct( ...
    'grfXmlFile', '', ...
    'bodyWeight', 78 * 9.81, ...
    'stanceThresholdBW', 0.02, ...
    'rightColumn', '', ...
    'leftColumn', '');
fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}) || isempty(opts.(fn{i}))
        opts.(fn{i}) = defaults.(fn{i});
    end
end
end

% ========================================================================
% 局部函数 2：识别右/左竖直 GRF 列名
% ========================================================================
% 优先级：
%   1) 用户显式 rightColumn/leftColumn；
%   2) ExternalLoads XML 的 <ExternalForce name="right/left"> +
%      force_identifier + 'y'；
%   3) SQR 约定 force2=right, force1=left；
%   4) 列名含 right/left 或 _r/_l；
%   5) 第一个/第二个 _vy，并 warning。
function [rightCol, leftCol] = resolve_force_columns(labels, xmlFile, opts)
if ~isempty(opts.rightColumn) && ~isempty(opts.leftColumn)
    rightCol = opts.rightColumn;
    leftCol  = opts.leftColumn;
    return;
end

rightCol = '';
leftCol  = '';

if ~isempty(xmlFile) && isfile(xmlFile)
    txt = fileread(xmlFile);
    blocks = regexp(txt, '(?s)<ExternalForce\s+name="([^"]+)">(.*?)</ExternalForce>', ...
        'tokens');
    for b = 1:numel(blocks)
        nm = blocks{b}{1};
        body = blocks{b}{2};
        fi = regexp(body, '<force_identifier>\s*([^<]+?)\s*</force_identifier>', ...
            'tokens', 'once');
        if isempty(fi)
            continue;
        end
        pref = strtrim(fi{1});
        cand = [pref, 'y'];
        if contains(lower(nm), 'right')
            rightCol = cand;
        elseif contains(lower(nm), 'left')
            leftCol = cand;
        end
    end
end
if ~isempty(rightCol) && ~isempty(leftCol)
    return;
end

vy = labels(endsWith(labels, '_vy'));
if isempty(vy)
    error('detect_stance_from_grf:noVerticalForce', ...
        'No column ending in "_vy" was found in the GRF table.');
end
r2 = vy(contains(vy, 'force2'));
l1 = vy(contains(vy, 'force1'));
if ~isempty(r2) && ~isempty(l1)
    rightCol = r2{1};
    leftCol  = l1{1};
    return;
end

rr = vy(contains(lower(vy), 'right') | contains(lower(vy), '_r'));
ll = vy(contains(lower(vy), 'left')  | contains(lower(vy), '_l'));
if ~isempty(rr), rightCol = rr{1}; end
if ~isempty(ll), leftCol  = ll{1}; end
if ~isempty(rightCol) && ~isempty(leftCol)
    return;
end

rightCol = vy{1};
if numel(vy) >= 2
    leftCol = vy{2};
else
    leftCol = vy{1};
end
warning('detect_stance_from_grf:columnGuess', ...
    'Could not identify right/left GRF reliably; using "%s" as right and "%s" as left.', ...
    rightCol, leftCol);
end

% ========================================================================
% 局部函数 3：严格按 cal_gait_phase.m 的规则计算 0..4
% ========================================================================
% 这个函数逐帧复现 cal_gait_phase.m 的 if/else 和 switch，不做任何
% 平滑/滞回/短 run 清理。如果某一帧无法落入 0..4，会报错并给出时间，
% 而不是静默猜测。
function [phaseR, phaseL] = phase_from_grf_com( ...
    grfR, grfL, pelvisCOM, calcnRCOM, calcnLCOM, stanceTh, time)

nT = numel(grfR);
phaseR = nan(nT, 1);
phaseL = nan(nT, 1);

% ---- 第一帧：cal_gait_phase.m 的初始化规则 ----------------------------
% 注意比较符号严格照抄：opposite GRF 用 > / <，不是 >= / <=。
if grfR(1) >= stanceTh
    if calcnRCOM(1) >= pelvisCOM(1)
        phaseR(1) = 0;
    elseif grfL(1) > stanceTh || calcnRCOM(1) + 1 < pelvisCOM(1)
        phaseR(1) = 2;
    elseif grfL(1) < stanceTh
        phaseR(1) = 1;
    end
elseif calcnRCOM(1) < pelvisCOM(1)
    phaseR(1) = 3;
elseif calcnRCOM(1) >= pelvisCOM(1)
    phaseR(1) = 4;
end

if grfL(1) >= stanceTh
    if calcnLCOM(1) >= pelvisCOM(1)
        phaseL(1) = 0;
    elseif grfR(1) > stanceTh || calcnLCOM(1) + 1 < pelvisCOM(1)
        phaseL(1) = 2;
    elseif grfR(1) < stanceTh
        phaseL(1) = 1;
    end
elseif calcnLCOM(1) < pelvisCOM(1)
    phaseL(1) = 3;
elseif calcnLCOM(1) >= pelvisCOM(1)
    phaseL(1) = 4;
end

assert_phase_valid(phaseR(1), phaseL(1), time(1));

% ---- 后续帧：状态转移规则 ---------------------------------------------
for k = 2:nT
    phaseR(k) = phaseR(k - 1);
    phaseL(k) = phaseL(k - 1);

    % right leg
    switch phaseR(k - 1)
        case 0
            if calcnRCOM(k) < pelvisCOM(k)
                phaseR(k) = 1;
            end
        case 1
            if grfL(k) >= stanceTh || calcnRCOM(k) + 1 < pelvisCOM(k)
                phaseR(k) = 2;
            end
        case 2
            if grfR(k) < stanceTh
                phaseR(k) = 3;
            end
        case 3
            if calcnRCOM(k) > pelvisCOM(k)
                phaseR(k) = 4;
            end
        case 4
            if grfR(k) >= stanceTh
                phaseR(k) = 0;
            end
    end

    % left leg
    switch phaseL(k - 1)
        case 0
            if calcnLCOM(k) < pelvisCOM(k)
                phaseL(k) = 1;
            end
        case 1
            if grfR(k) >= stanceTh || calcnLCOM(k) + 1 < pelvisCOM(k)
                phaseL(k) = 2;
            end
        case 2
            if grfL(k) < stanceTh
                phaseL(k) = 3;
            end
        case 3
            if calcnLCOM(k) > pelvisCOM(k)
                phaseL(k) = 4;
            end
        case 4
            if grfL(k) >= stanceTh
                phaseL(k) = 0;
            end
    end

    assert_phase_valid(phaseR(k), phaseL(k), time(k));
end
end

% ========================================================================
% 局部函数 4：相位合法性检查
% ========================================================================
function assert_phase_valid(phR, phL, t)
if ~(isscalar(phR) && isfinite(phR) && phR >= 0 && phR <= 4 && ...
     isscalar(phL) && isfinite(phL) && phL >= 0 && phL <= 4)
    error('detect_stance_from_grf:invalidPhase', ...
        ['cal_gait_phase rules produced an invalid phase at t = %.6f s ' ...
         '(right=%g, left=%g). Check GRF/COM data.'], t, phR, phL);
end
end
