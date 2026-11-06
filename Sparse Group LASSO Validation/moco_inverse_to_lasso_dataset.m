% =========================================================================
% moco_inverse_to_lasso_dataset.m
%
% 作用
% ----
% 把“实测 Vicon/GRF/sEMG -> MocoInverse”得到的肌肉结果，转换成旧 LASSO
% 流程可以直接读取的数据集。输出包括：
%   * moco_lasso_dataset.sto : 与 SCONE UN.sto 同列名，供
%     extractStoMuscleFeatures.m / Adaptive_LASSO_MocoInverse.m 使用；
%   * moco_lasso_dataset.mat : 原始/清洗后信号、特征、相位、权重、来源、
%     QC 信息，便于复现和画图；
%   * moco_lasso_dataset_qc_*.png / .csv / _report.txt : 质量检查。
%
% 为什么需要这个接口
% ------------------
% 旧 LASSO 输入 UN.sto 是 SCONE 仿真输出，已经包含：
%     <muscle>.excitation
%     <muscle>.fiber_length_norm
%     <muscle>.mtu_force_norm
%     leg1_r.state
% 但 MocoInverse（尤其是初步 --sagittal-cost 解）只有：
%     /forceset/<muscle>/activation
% 没有坐标（在 ik_radians.mot）、没有 fiber_length_norm/mtu_force_norm、
% 没有 SCONE 的 leg1_r.state。因此本脚本补三类东西：
%   1. 运动学：优先用 Moco 解里的 /jointset/.../value|speed，否则读 IK；
%   2. 肌肉特征：用预测仿真模型 human0918，沿给定 q/qdot/activation
%      积分 fiber_length，得到 lCEN/fATN；
%   3. 相位：从实测竖直 GRF 重建 stance/swing（并写近似 0..4 标签）。
%
% 处理链（对应下面 [1]–[5]）
% --------------------------
%   [1] 读 Moco 解，取 activation/excitation；可选额外 excitation/activation
%       文件；确定所有输入文件共同覆盖的时间范围，重采样到 targetDt；
%   [2] 清洗 muscle drives（despike + 零相位低通；target 始终是 raw excitation），
%       详见 clean_excitation_signals；
%   [3] 用 human0918 计算 18 块肌肉的 lCEN/fATN，并对特征再做低通；
%   [4] 从 GRF 重建 stance/swing 与 0..4 标签；
%   [5] 组装 ds 结构，写 .sto/.mat/QC CSV/图，并给 LASSO main 留 sampleWeight。
%
% 重要设计选择
% ------------
% * 特征必须用“预测仿真实际使用的模型”计算（默认 human0918），而不是
%   Moco 的 SQR_simbody。两者的 Fmax 现已对齐（human0918 已改为原来的
%   2 倍），但 lopt/tendon slack/几何仍有差异；fATN 是模型相关量，不能混用。
% * 默认 targetDt=0.01 s：100 Hz 对线性反射模型足够，也让 30–60 s 数据的
%   纤维长度积分和 LASSO 拟合可控。
% * validMask/sampleWeight：首尾 edgeMargin、非有限样本、用户指定
%   excludeTimeWindows 会被标 0；LASSO 用权重排除，不删行，保证 lambda
%   path 的时间连续性。
% * 旧 UN.sto -> Adaptive_LASSO_StanceSwing.m 路径完全不动；本脚本只新增。
%
% 运行方法
% --------
%   1) 修改下面 CONFIG 里的 solutionFile/ikFile/grfFile/modelFile 等；
%   2) 在 MATLAB 里运行：
%          moco_inverse_to_lasso_dataset
%   3) 再运行：
%          Adaptive_LASSO_MocoInverse
%
% 输出目录
% --------
%   Sparse Group LASSO Validation/moco_dataset/
% =========================================================================

% -------------------------------------------------------------------------
% 0. 路径与日志
% -------------------------------------------------------------------------
% 脚本会从自身位置推导 repoRoot，因此无论从仓库根目录还是从
% Sparse Group LASSO Validation 目录运行，默认路径都正确。
% diary 同时把控制台输出写到 report 文件，便于把一次接口运行归档。
close all; clear; clc;

scriptFullPath = mfilename('fullpath'); % Sparse Group LASSO Validation\
if isempty(scriptFullPath)
    scriptDir = pwd;
else
    scriptDir = fileparts(scriptFullPath);
end
addpath(scriptDir); 
repoRoot = fileparts(scriptDir); % Muscle-Reflex-Based Predictive Simulation\
outDir = fullfile(scriptDir, 'moco_dataset');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
outStem = 'moco_lasso_dataset';
reportFile = fullfile(outDir, [outStem '_report.txt']);
if isfile(reportFile)
    delete(reportFile);
end
diary(reportFile);
diary on;
logCleanup = onCleanup(@() diary('off'));

fprintf('=====================================================================\n');
fprintf(' MocoInverse -> LASSO dataset builder\n');
fprintf(' repository : %s\n', repoRoot);
fprintf(' output dir : %s\n', outDir);
fprintf(' started    : %s\n', char(datetime));
fprintf('=====================================================================\n');

%% ======================= CONFIG (edit here) ============================
% 所有与“这次数据”相关的路径/参数都集中在这里。换受试者/换 Moco 结果时
% 通常只需要改这一块，不需要动下面的处理逻辑。
%
% 重要：solutionFile 必须同时包含 excitation（/forceset/<muscle> 控制列）
% 和 activation（/forceset/<muscle>/activation state）。推荐直接指向正式
% MocoInverse 的完整解：
%   * 完整解（activation + excitation + 可选 /jointset 坐标）；
%   * 或“完整解 + 单独 excitationFile/activationFile”的组合。
%
% 注意：下面默认的 B30 preliminary 解只有 activation、没有 excitation，
% 在新逻辑下会直接报 MissingExcitation。运行正式数据前请把这里改成包含
% 完整 excitation + activation 的最终解。
solutionFile = fullfile(repoRoot, 'experiments', 'moco', 'outputs', ...
    'moco_inverse_solution_B30_right_final.sto');

% Optional separate files.  Leave empty when the solution already contains
% the corresponding columns.
excitationFile = '';
activationFile = '';

% Kinematics.  Required when the solution has no /jointset/.../value columns.
% ik_radians.mot is filtered at 15 Hz and already in radians.
ikFile = fullfile(repoRoot, 'experiments', 'moco', 'outputs', 'ik_radians.mot');

% Measured GRF, used only for stance/swing detection.  The matching
% ExternalLoads XML is auto-detected next to the .mot file when present.
grfFile = fullfile(repoRoot, 'experiments', 'moco', 'outputs', ...
    'grf_window_209p5_240.mot');
grfXmlFile = '';

% Model used by the predictive simulation (feature definition must match it).
modelFile = fullfile(repoRoot, 'model', 'human0918.osim');

% Output time step for the LASSO dataset [s].  Moco runs at 1-10 ms; 10 ms is
% enough for a linear reflex model and keeps the fiber-length integration and
% the subsequent LASSO fit fast.
targetDt = 0.01;

% Muscle order: canonical right leg first, then left leg.  This is also the
% column order used by human0918.osim / lasso beta columns.
muscleNames = { ...
    'hamstrings_r','bifemsh_r','glut_max_r','iliopsoas_r','rect_fem_r', ...
    'vasti_r','gastroc_r','soleus_r','tib_ant_r', ...
    'hamstrings_l','bifemsh_l','glut_max_l','iliopsoas_l','rect_fem_l', ...
    'vasti_l','gastroc_l','soleus_l','tib_ant_l'};

% ---- excitation cleaning ----------------------------------------------
% 正式 MocoInverse 结果一定包含 /forceset/<muscle> 的 excitation control；
% LASSO 的 target 始终是“清洗后的 raw excitation”。
% 不再提供 activation_consistent / activation 等回退方案：
%   * 如果解里没有 excitation，接口直接报错；
%   * aRaw 只用于计算 fATN/lCEN，不参与 target 构造。
cleanOpts = struct();
cleanOpts.lowpassHz      = 6;       % 0 disables
cleanOpts.filterOrder    = 4;
cleanOpts.doDespike      = true;
cleanOpts.medianWindow   = 5;
cleanOpts.hampelSigma    = 4;
cleanOpts.maxBurstSamples = 5;
cleanOpts.floor          = 0.01;    % must match LASSO censoring floor
cleanOpts.ceiling        = 1.0;
cleanOpts.tauAct         = 0.01;    % DGF / Millard activation time constant
cleanOpts.tauDeact       = 0.04;

% Low-pass applied to the computed muscle features before LASSO.  Features
% are smoother than excitations, but fATN inherits activation noise.
featureFilterHz = 6;

% ---- GRF / phase settings ---------------------------------------------
phaseOpts = struct();
phaseOpts.bodyWeight        = 78 * 9.81;
% 必须与 forward_simulation.m 的 simCache.stanceTh 一致：0.02 * BW。
% 不再使用 enter/exit 双阈值或 minStance/minSwingGap 清理。
phaseOpts.stanceThresholdBW = 0.02;

% Samples to exclude from the LASSO fit (seconds).  Useful for Moco chunk
% seams (e.g. B30 chunks at 217.5 / 225 / 232.5 s) or known bad windows.
excludeTimeWindows = [];   % e.g. [217.45 217.55; 224.95 225.05]
edgeMargin = 0.05;         % [s] trimmed at both ends

% Feature integration options.
featureOpts = struct();
featureOpts.substepMax = 0.002;
featureOpts.equilibrateInitial = true;
featureOpts.verbose = true;

makeQcPlots = true;
%% =======================================================================

%% ---- 1. read MocoInverse solution ---------------------------------------
if ~isfile(solutionFile)
    error('MOSolutionNotFound', 'MocoInverse solution not found: %s', solutionFile);
end
S = read_moco_inverse_solution(solutionFile, muscleNames);
fprintf('\n[1] solution : %s\n', solutionFile);
fprintf('    time range in solution : %.4f .. %.4f s (%d samples)\n', ...
    S.time(1), S.time(end), numel(S.time));
fprintf('    activation present     : %d\n', S.hasActivation);
fprintf('    excitation present     : %d\n', S.hasExcitation);
fprintf('    coordinates in solution: %d\n', S.hasCoordinates);

% Optional extra files (only used when the solution lacks a signal).
extraExc = [];
extraAct = [];
if ~isempty(excitationFile)
    E = read_moco_inverse_solution(excitationFile, muscleNames);
    if E.hasExcitation
        extraExc = E;
        fprintf('    extra excitation file  : %s\n', excitationFile);
    end
end
if ~isempty(activationFile)
    A = read_moco_inverse_solution(activationFile, muscleNames);
    if A.hasActivation
        extraAct = A;
        fprintf('    extra activation file  : %s\n', activationFile);
    end
end

% -------------------------------------------------------------------------
% 1b. 确定公共时间网格
% -------------------------------------------------------------------------
% 所有输入（Moco 解、可选 excitation/activation 文件、IK、GRF）必须取交集：
% 否则 interp1 会在某段外推，而外推出来的 GRF/运动学没有物理意义。
% 先求最大起点、最小终点，再按 targetDt 生成均匀网格。
%
% 注意：如果 Moco 解自己带 /jointset 坐标，就不需要把 IK 纳入交集；
% 这样即使 IK 只覆盖部分时间，也不会无故裁掉 Moco 解。
starts = S.time(1);
ends   = S.time(end);
if ~isempty(extraExc), starts = max(starts, extraExc.time(1)); ends = min(ends, extraExc.time(end)); end
if ~isempty(extraAct), starts = max(starts, extraAct.time(1)); ends = min(ends, extraAct.time(end)); end
if ~S.hasCoordinates && isfile(ikFile)
    IK0 = read_sto_mot_table(ikFile);
    starts = max(starts, IK0.time(1));
    ends   = min(ends, IK0.time(end));
end
if isfile(grfFile)
    G0 = read_sto_mot_table(grfFile);
    starts = max(starts, G0.time(1));
    ends   = min(ends, G0.time(end));
end
if ~isfinite(starts) || ~isfinite(ends) || ends <= starts
    error('BadTimeRange', 'No overlapping time range for solution/IK/GRF.');
end
tOut = (starts:targetDt:ends).';
nT = numel(tOut);
fprintf('    LASSO time grid        : %.4f .. %.4f s, dt = %.4f s, n = %d\n', ...
    tOut(1), tOut(end), targetDt, nT);
if nT < 50
    warning('ShortDataset', 'Only %d samples; LASSO fit may be unstable.', nT);
end

%% ---- 2. activation / excitation onto the common grid -------------------
% 正式 MocoInverse 解必须同时提供 excitation 和 activation：
%   * excitation 是 LASSO target；
%   * activation 用于计算 fATN/lCEN。
% 两者都不允许缺失；缺失/NaN 会在下面 assert_complete_muscle_signal
% 直接报错，不用 NaN 填充，也不从另一个信号反演/积分补全。
% 这里用 interp_matrix 而不是直接赋值，是因为 Moco 输出网格（1–10 ms）
% 可能和 targetDt 不完全一致；INTERP 逐列做。
aRaw = [];
xRaw = [];
% 在插值之前先检查源矩阵：如果主解/替代文件本身缺列或含 NaN/Inf，
% 直接报错。这样 interp_matrix 不会把“部分 NaN”悄悄插值补掉。
if S.hasActivation
    assert_complete_muscle_signal(S.activation, muscleNames, 'Activation');
    aRaw = interp_matrix(S.time, S.activation, tOut);
elseif ~isempty(extraAct)
    assert_complete_muscle_signal(extraAct.activation, muscleNames, 'Activation');
    aRaw = interp_matrix(extraAct.time, extraAct.activation, tOut);
end
if S.hasExcitation
    assert_complete_muscle_signal(S.excitation, muscleNames, 'Excitation');
    xRaw = interp_matrix(S.time, S.excitation, tOut);
elseif ~isempty(extraExc)
    assert_complete_muscle_signal(extraExc.excitation, muscleNames, 'Excitation');
    xRaw = interp_matrix(extraExc.time, extraExc.excitation, tOut);
end
% 再检查插值后的目标网格；这里是最终防线，确保没有 NaN 进入 LASSO。
assert_complete_muscle_signal(xRaw, muscleNames, 'Excitation');
assert_complete_muscle_signal(aRaw, muscleNames, 'Activation');

%% ---- 3. kinematics ------------------------------------------------------
% 两条分支：
%   (a) Moco 解里有 /jointset/.../value|speed：直接用，它已经是 Moco
%       实际约束/优化的坐标轨迹，单位是 rad/m；
%   (b) 没有坐标（初步 --sagittal-cost）：读 IK。
%       IK 可能是 inDegrees=yes，必须只对转动坐标做 deg2rad；
%       qdot 用中心差分（gradient），因为 IK 通常只有坐标没有速度。
%
% 无论哪条分支，最后都统一成 qOut/qdOut，且列名 qNamesOut 与列一一对应；
% build_lasso_features 会按模型坐标名去匹配，不要求列顺序完全一样。
if S.hasCoordinates
    qOut  = interp_matrix(S.time, S.q, tOut);
    qNamesOut = S.coordNames;
    if S.hasCoordinateVel && any(isfinite(S.qdot(:)))
        qdOut = interp_matrix(S.time, S.qdot, tOut);
    else
        qdOut = derivative_matrix(tOut, qOut);
    end
    % If a speed channel was missing, fall back to the numerical derivative
    % of the corresponding value channel.
    badQd = all(~isfinite(qdOut), 1);
    if any(badQd)
        qdFallback = derivative_matrix(tOut, qOut);
        qdOut(:, badQd) = qdFallback(:, badQd);
    end
    fprintf('    kinematics source      : Moco solution /jointset values\n');
else
    if ~isfile(ikFile)
        error('NoKinematics', ...
            'Solution has no coordinates and ikFile is empty/not found: %s', ikFile);
    end
    IK = read_sto_mot_table(ikFile);
    qNamesOut = IK.labels(2:end);
    qIK = IK.data(:, 2:end);
    if IK.inDegrees
        isTrans = endsWith(qNamesOut, {'_tx','_ty','_tz'});
        for c = 1:numel(qNamesOut)
            if ~isTrans(c)
                qIK(:, c) = deg2rad(qIK(:, c));
            end
        end
    end
    qdIK = derivative_matrix(IK.time, qIK);
    qOut  = interp_matrix(IK.time, qIK, tOut);
    qdOut = interp_matrix(IK.time, qdIK, tOut);
    fprintf('    kinematics source      : %s (filtered IK)\n', ikFile);
end

if any(~isfinite(qOut(:)))
    error('BadKinematics', 'Interpolated coordinate values contain NaN/Inf.');
end
if any(~isfinite(qdOut(:)))
    error('BadKinematics', 'Interpolated coordinate speeds contain NaN/Inf.');
end

%% ---- 4. clean muscle drives --------------------------------------------
% clean_excitation_signals 负责：
%   * 非有限值插值；
%   * Hampel 短 spike 替换；
%   * 零相位低通；
%   * 输出 xClean（target，来自 raw excitation）和 aClean（只用于特征）。
% 这里不再从 activation 反演 excitation；如果 xRaw 缺失，clean 函数会报错。
fprintf('\n[2] cleaning muscle drives (target = raw excitation, low-pass = %.2f Hz)\n', ...
    cleanOpts.lowpassHz);
[xClean, aClean, qcClean] = clean_excitation_signals(tOut, xRaw, aRaw, cleanOpts);
fprintf('    raw excitation present : %d\n', qcClean.hasRawExcitation);
fprintf('    raw activation present : %d\n', qcClean.hasRawActivation);
fprintf('    activation recon RMSE  : mean %.4f, max %.4f\n', ...
    mean(qcClean.activationReconRMSE), max(qcClean.activationReconRMSE));
fprintf('    flagged excitation %%   : mean %.2f, max %.2f\n', ...
    mean(qcClean.pctFlaggedExcitation), max(qcClean.pctFlaggedExcitation));


%% ---- 5. muscle features (lCEN, fATN) -----------------------------------
% 用 aClean（而不是 raw activation）计算 fATN，保持“清洗后的 activation ->
% 特征 -> target”三者一致。模型默认 human0918，原因见文档：
% 虽然 human0918 与 SQR 的 Fmax 已对齐，且 fATN 对均匀 Fmax 缩放不敏感，
% 但两者的 lopt/tendon slack/几何仍不同，fATN 仍可能有差异；控制器最终在
% human0918 里运行，所以特征必须由 human0918 计算。
%
% 特征出来后再次低通：fATN 对 activation 的高频很敏感，6 Hz 低通能显著
% 降低 LASSO 设计矩阵的噪声，同时不改变步态主频（~1 Hz）和反射延迟结构。
fprintf('\n[3] computing muscle features with %s\n', modelFile);
F = build_lasso_features_from_kinematics(modelFile, tOut, qOut, qdOut, ...
    aClean, muscleNames, qNamesOut, featureOpts);

if featureFilterHz > 0
    fs = 1 / targetDt;
    F.lCEN = zero_phase_lowpass(F.lCEN, fs, featureFilterHz, cleanOpts.filterOrder);
    F.fATN = zero_phase_lowpass(F.fATN, fs, featureFilterHz, cleanOpts.filterOrder);
    fprintf('    feature low-pass       : %.2f Hz\n', featureFilterHz);
end

%% ---- 6. phase labels from GRF ------------------------------------------
% 严格按 functions/core/cal_gait_phase.m 的规则重建 0..4：
%   * 使用实测竖直 GRF（阈值 0.02*BW，与 forward_simulation.m 一致）；
%   * 使用 pelvis / calcn_r / calcn_l 质心的全局 X 坐标；
%   * 不再使用滞回、短 run 清理或 60/40、20/60/20 启发式。
% COM 由 compute_com_x_from_kinematics 用同一个模型和 q/qdot 计算。
fprintf('\n[4] reconstructing gait phases from GRF + pelvis/calcn COM\n');
phaseOpts.grfXmlFile = grfXmlFile;
comX = compute_com_x_from_kinematics(modelFile, tOut, qOut, qdOut, qNamesOut);
P = detect_stance_from_grf(grfFile, tOut, comX, phaseOpts);
fprintf('    GRF file               : %s\n', P.file);
fprintf('    right vertical column  : %s\n', P.rightColumn);
fprintf('    left  vertical column  : %s\n', P.leftColumn);
fprintf('    stance threshold       : %.4f BW = %.2f N\n', ...
    P.thresholds.stanceThresholdBW, P.thresholds.stanceThresholdN);
fprintf('    right stance fraction  : %.3f\n', mean(P.stanceR));
fprintf('    left  stance fraction  : %.3f\n', mean(P.stanceL));
fprintf('    right phase-0 onsets   : %d\n', numel(P.eventsR));

%% ---- 7. valid mask / sample weights ------------------------------------
% 为什么要 sampleWeight 而不是删行：
%   * LASSO 的 lambda path / censoring 假设样本在时间上连续；
%   * 删除 Moco chunk seam 会让前后样本突然相邻，特征差分/相位索引会错位；
%   * 用 0 权重可以“逻辑上排除”这些样本，同时保留数组长度和索引。
%
% 这里排除三类：
%   1) 首尾 edgeMargin：fiber_length 初始 transient / 末端边界效应；
%   2) 任何 xClean/aClean/lCEN/fATN 非有限的样本；
%   3) 用户显式给的 excludeTimeWindows（例如 B30 的 217.5/225/232.5 s）。
validMask = true(nT, 1);
validMask(1:min(nT, max(1, round(edgeMargin / targetDt)))) = false;
validMask(max(1, nT - round(edgeMargin / targetDt) + 1):end) = false;

allFinite = all(isfinite([xClean, aClean, F.lCEN, F.fATN]), 2);
validMask = validMask & allFinite;

if ~isempty(excludeTimeWindows)
    for r = 1:size(excludeTimeWindows, 1)
        validMask(tOut >= excludeTimeWindows(r, 1) & ...
                  tOut <= excludeTimeWindows(r, 2)) = false;
    end
end
if nnz(validMask) < 100
    warning('FewValidSamples', ...
        'Only %d valid samples; check excludeTimeWindows / edgeMargin.', nnz(validMask));
end
sampleWeight = double(validMask);

%% ---- 8. assemble and save ----------------------------------------------
% ds 是接口的“完整档案”：既包含写 .sto 需要的列，也包含 raw/clean 对照、
% 相位、权重、来源、QC。LASSO main 会从 .mat 里取 sampleWeight 和 meta，
% 从 .sto 里取特征/target（走旧 reader），两者互补。
ds = struct();
ds.meta = struct( ...
    'created', char(datetime), ...
    'solutionFile', solutionFile, ...
    'excitationFile', excitationFile, ...
    'activationFile', activationFile, ...
    'ikFile', ikFile, ...
    'grfFile', grfFile, ...
    'grfXmlFile', P.xmlFile, ...
    'modelFile', modelFile, ...
    'targetDt', targetDt, ...
    'targetSource', qcClean.targetSource, ...
    'cleanOpts', cleanOpts, ...
    'featureFilterHz', featureFilterHz, ...
    'phaseOpts', phaseOpts, ...
    'excludeTimeWindows', excludeTimeWindows, ...
    'edgeMargin', edgeMargin, ...
    'unitConvention', 'radians, metres, seconds');
ds.time        = tOut;
ds.muscleNames = muscleNames;
ds.rawActivation = aRaw;   % pre-cleaning, complete 18-muscle matrix
ds.rawExcitation = xRaw;   % pre-cleaning, complete 18-muscle matrix
ds.activation  = aClean;
ds.excitation  = xClean;
ds.lCEN        = F.lCEN;
ds.fATN        = F.fATN;
ds.q           = qOut;
ds.qdot        = qdOut;
ds.qNames      = qNamesOut;
ds.phase = struct( ...
    'stanceR', P.stanceR, 'stanceL', P.stanceL, ...
    'contactR', P.contactR, 'contactL', P.contactL, ...
    'stateR', P.stateR, 'stateL', P.stateL, ...
    'forceR', P.forceR, 'forceL', P.forceL, ...
    'thresholds', P.thresholds, ...
    'eventsR', P.eventsR, 'eventsL', P.eventsL, ...
    'comX', P.comX);
ds.validMask    = validMask;
ds.sampleWeight = sampleWeight;
ds.qc = struct('cleaning', qcClean, 'features', F.diagnostics, 'phase', P.thresholds);

stoFile = fullfile(outDir, [outStem '.sto']);
matFile = fullfile(outDir, [outStem '.mat']);
write_lasso_dataset_sto(stoFile, ds);
save(matFile, 'ds', '-v7.3');

% ---- QC CSV + figures ---------------------------------------------------
qcTable = array2table([qcClean.stats, (1:numel(muscleNames)).'], ...
    'VariableNames', [qcClean.statNames, {'muscleIndex'}]);
qcTable = addvars(qcTable, muscleNames(:), 'Before', 1, 'NewVariableNames', 'muscle');
writetable(qcTable, fullfile(outDir, [outStem '_qc.csv']));

if makeQcPlots
    try
        plot_moco_lasso_dataset_qc(ds, outDir, outStem);
    catch ME
        warning('QcPlotFailed', 'QC plotting failed: %s', ME.message);
    end
end

fprintf('\n[5] saved\n');
fprintf('    dataset .sto : %s\n', stoFile);
fprintf('    dataset .mat : %s\n', matFile);
fprintf('    QC csv       : %s\n', fullfile(outDir, [outStem '_qc.csv']));
fprintf('    report       : %s\n', reportFile);
fprintf('\nNext step: run Adaptive_LASSO_MocoInverse.m\n');
fprintf('=====================================================================\n');

diary off;

%% ======================= local helpers ==================================
% 局部函数 1：检查 excitation/activation 是否覆盖全部肌肉
% -------------------------------------------------------------------------
% 用户要求：如果 Moco 结果缺少部分肌肉，不要用 NaN 填充，直接报错。
% 这里同时检查：
%   * M 是否为空；
%   * 列数是否等于 muscleNames 的数量；
%   * 哪些肌肉整列为 NaN（完全缺失）；
%   * 哪些肌肉含 NaN/Inf（部分缺失/求解失败）。
% 报错信息会列出具体肌肉名，方便回去重导 Moco 解。
function assert_complete_muscle_signal(M, names, label)
if isempty(M)
    error(['Missing' label], ...
        '%s is empty. Expected columns for: %s.', ...
        label, strjoin(names, ', '));
end
if size(M, 2) ~= numel(names)
    error(['Missing' label], ...
        '%s has %d columns, expected %d (%s).', ...
        label, size(M, 2), numel(names), strjoin(names, ', '));
end

missing = names(~any(isfinite(M), 1));
if ~isempty(missing)
    error(['Missing' label], ...
        '%s is completely missing for: %s.', ...
        label, strjoin(missing, ', '));
end

bad = names(any(~isfinite(M), 1));
if ~isempty(bad)
    error(['NonFinite' label], ...
        '%s contains NaN/Inf for: %s.', ...
        label, strjoin(bad, ', '));
end
end

% 局部函数 2：多列线性插值
% -------------------------------------------------------------------------
% 为什么不用 interp1(tSrc, M, tOut)：
%   * interp1 对矩阵默认按列插值，但要求每列有相同的时间点；
%   * 我们的 M 可能包含整列 NaN（缺失通道），interp1 会报错；
%   * 这里逐列插值并保留全 NaN 列，随后由
%     assert_complete_muscle_signal 直接报错，不做 NaN 填充。
%
% 输入：
%   tSrc  nTsrc×1
%   M     nTsrc×nCol
%   tOut  nTout×1
% 输出：
%   Y     nTout×nCol
function Y = interp_matrix(tSrc, M, tOut)
%INTERP_MATRIX  Column-wise linear interpolation, NaN-safe at edges.
if isempty(M)
    Y = [];
    return;
end
tSrc = tSrc(:);
nCol = size(M, 2);
Y = nan(numel(tOut), nCol);
[tSrc, ord] = sort(tSrc);
M = M(ord, :);
for c = 1:nCol
    v = M(:, c);
    good = isfinite(v);
    if ~any(good)
        continue;
    end
    if nnz(good) == 1
        Y(:, c) = v(find(good, 1));
    else
        Y(:, c) = interp1(tSrc(good), v(good), tOut, 'linear', 'extrap');
    end
end
end

% 局部函数 3：逐列中心差分
% -------------------------------------------------------------------------
% gradient 对非均匀 t 也按中心差分处理（端点单边）。逐列调用是为了避免
% gradient(Q,t) 在矩阵/向量 t 组合下的维度歧义。
% 输出单位：Q 是 rad/m，t 是 s -> D 是 rad/s 或 m/s。
function D = derivative_matrix(t, Q)
%DERIVATIVE_MATRIX  Central-difference speed, one column at a time.
nT = size(Q, 1);
D = nan(size(Q));
for c = 1:size(Q, 2)
    D(:, c) = gradient(Q(:, c), t);
end
end
