function [fit, lasso] = fit_adaptive_lasso_phase_model(stoFilePath, options)
%FIT_ADAPTIVE_LASSO_PHASE_MODEL  按步态相位拟合 adaptive LASSO 肌肉反射模型。
%
% 这个函数做什么
% --------------
% 输入一个 LASSO-ready 的 .sto（右腿 9 块肌肉的 activation/excitation +
% fiber_length_norm/mtu_force_norm + 关节运动学 + leg1_r.state），对每个
% 相位（默认 stance/swing）分别拟合：
%
%     y_m(t) = A(t) * beta_{p,m} + b_{p,m} + noise
%
% 其中 A(t) 是 26 维特征（9 lCEN + 9 fATN + 8 关节角度/角速度），
% y 是目标 excitation，m 是肌肉。最终输出可直接给
% load_lasso_reflex_controller.m 的 grouped lasso 结构。
%
% 算法与旧 main 的关系
% --------------------
% 基本流程完全对齐 Adaptive_LASSO_StanceSwing.m：
%   1. 按相位取样本；
%   2. 在 uncensored 样本上做（默认 Huber-robust）加权 OLS 初估；
%   3. 自适应权重 w_j = 1 / max(|beta_OLS_j|, eps)^gamma；
%   4. 列缩放 A* = A ./ w，把 adaptive LASSO 变成标准 LASSO；
%   5. 沿 lambda path 用 sgl_fit 求解，并在每个 lambda 上做 support refit；
%   6. 按 refit RMSE 选 lambda（最小 RMSE 附近取最稀疏）；
%   7. 在选定 support 上做 post-selection OLS（去偏）。
%
% 针对实测 Moco 数据新增两点：
%   * sampleWeight：把 Moco 分段接缝、首尾边缘、无效样本标 0，但不删行；
%     sgl_fit 的 uncensored 平方项和 censored hinge 项都会乘权重。
%   * Huber IRLS：初始 OLS（为 adaptive 权重服务）和最终 refit 都用
%     Huber 权重迭代，避免少量 excitation burst 的大残差把系数/权重拉偏。
%
% 输入
% ----
%   stoFilePath : LASSO-ready .sto；列名要求见 write_lasso_dataset_sto.m，
%                 或任意 SCONE 风格、含 *_r.activation /
%                 *.fiber_length_norm / *.mtu_force_norm / leg1_r.state 的文件。
%   options     : struct，字段及默认值见 apply_defaults 和下面的说明。
%
% 输出
% ----
%   fit   : 完整拟合结果：
%           .Beta_refit / .Intercept_refit      最终（robust）系数
%           .YPred_refit / .RMSE_refit           预测与误差
%           .LambdaPath                         每个肌肉/相位的 lambda path
%           .S / .sampleWeight / .options        结构信息
%   lasso : grouped 控制器：
%           .groups(g).beta/bias/mask/phases/label
%           以及 featureNames/muscleNames/source/preprocessing 元数据。
%
% options 主要字段（默认值）
% -------------------------
%   gamma 1, epsilonOLS 1e-8, nzTol 1e-5, lambda1 0
%   lambdaPathConfig.nLambda 60
%   lambdaPathConfig.lambdaMinRatio 1e-4
%   lambdaPathConfig.tolSelect 0.05
%   c 0.01, rho 50, delaySteps 1, excitationThreshold 0.01
%   sampleWeight []（全 1）
%   robustInitIter 2（初始 OLS 的 Huber IRLS；0 关闭）
%   robustIter 2（post-selection refit 的 Huber IRLS）
%   huberDeltaFactor 1.345
%   verbose true
%   groupLabels {'stance','swing'}
%   groupPhases {[0 1], [2 3 4]}
%   source 'MocoInverse'
%
% 数值/实现注意
% ------------
% * c 是 censoring floor，必须与接口 cleanOpts.floor 一致；
% * delaySteps 在 preprocessExcitationAndFeatures 中把 target 前移、
%   features 后移；sampleWeight 和 phaseIdx 必须做同样的延迟裁剪；
% * lambda path 的每个点都做 refit，计算量较大，但这是旧流程的设定；
% * 保存的 fit 里 lambda path 的 stats 只保留 compact 版本，避免几十 MB
%   的 FISTA objHist 把 .mat 撑爆。

if nargin < 1 || isempty(stoFilePath)
    error('fit_adaptive_lasso_phase_model:badInput', 'stoFilePath is required.');
end
if nargin < 2 || isempty(options)
    options = struct();
end
opts = apply_defaults(options);
stoFilePath = char(stoFilePath);
if ~isfile(stoFilePath)
    error('fit_adaptive_lasso_phase_model:notFound', ...
        'STO file not found: %s', stoFilePath);
end

%% ===== STEP 1: DATA LOADING & PREPROCESSING =====
% 1a. 从 .sto header 自动识别右腿肌肉
% -------------------------------------------------------------------------
% 旧流程假定“右腿 9 块肌肉”按 canonical 顺序出现；这里仍然只拟合右腿，
% 并假定双侧对称（与旧 main / demo 的假设一致）。
% 用 endsWith('_r.activation') 而不是写死肌肉名，是为了兼容 SCONE UN.sto
% 和新接口生成的 moco_lasso_dataset.sto；两者的列名相同。
fid = fopen(stoFilePath, 'r');
if fid == -1
    error('fit_adaptive_lasso_phase_model:cannotOpen', ...
        'Cannot open STO file: %s', stoFilePath);
end
ln = '';
while ischar(ln) && ~contains(ln, 'endheader')
    ln = fgetl(fid);
end
if ~ischar(ln)
    fclose(fid);
    error('fit_adaptive_lasso_phase_model:noEndHeader', ...
        'endheader not found in STO file: %s', stoFilePath);
end
headerLine = fgetl(fid);
while ischar(headerLine) && isempty(strtrim(headerLine))
    headerLine = fgetl(fid);
end
fclose(fid);
if ~ischar(headerLine)
    error('fit_adaptive_lasso_phase_model:noHeader', ...
        'No column-label line found in %s.', stoFilePath);
end
colNames = regexp(strtrim(headerLine), '\s+', 'split');
colNames = colNames(~cellfun(@isempty, colNames));

isRightActivation = endsWith(colNames, '_r.activation');
muscleIdx = find(isRightActivation);
if isempty(muscleIdx)
    error('fit_adaptive_lasso_phase_model:noMuscles', ...
        'No "<muscle>_r.activation" columns found in %s.', stoFilePath);
end
muscleNames = cell(1, numel(muscleIdx));
for i = 1:numel(muscleIdx)
    muscleNames{i} = extractBefore(colNames{muscleIdx(i)}, '.activation');
end
numTargets = numel(muscleNames);
fprintf('[LASSO] %d right-leg muscles detected in %s\n', numTargets, stoFilePath);
fprintf('        %s\n', strjoin(muscleNames, ', '));

gamma = opts.gamma;
epsilon_ols = opts.epsilonOLS;
nzTol = opts.nzTol;
lambda1 = opts.lambda1;
lambdaPathConfig = opts.lambdaPathConfig;
c = opts.c;
rho = opts.rho;
delaySteps = opts.delaySteps;
excitationThreshold = opts.excitationThreshold;

assert(lambda1 == 0, ...
    'Adaptive LASSO column-scaling only adapts L1 penalty. lambda1 must be 0.');

% 1b. 用第一块肌肉建立共享特征矩阵 X 和 mask M
% -------------------------------------------------------------------------
% extractStoMuscleFeatures 会读：
%   <muscle>_r.excitation（target）
%   <muscle>_r.fiber_length_norm / mtu_force_norm（特征）
%   9 个右腿肌肉的 *_r 列 + 8 个关节/骨盆列
% 后面的肌肉会复用同一个 X/M；如果某个肌肉的 X/M 与 base 不完全一致，
% 说明 .sto 数据有错位，直接报错比“悄悄用错特征”安全。
[yAllBase, XAll, featureToIndexMap, M, phase] = ...
    extractStoMuscleFeatures(stoFilePath, muscleNames{1});
[yBase, X] = preprocessExcitationAndFeatures(yAllBase, XAll, delaySteps, excitationThreshold);

Y = zeros(numel(yBase), numTargets);
Y(:, 1) = yBase;
for m = 2:numTargets
    [yAllM, XAllM, featureToIndexMapM, MM, ~] = ...
        extractStoMuscleFeatures(stoFilePath, muscleNames{m});
    if ~isequal(size(XAllM), size(XAll)) || max(abs(XAllM(:) - XAll(:))) > 1e-12
        error('fit_adaptive_lasso_phase_model:FeatureMismatch', ...
            'Feature matrix for %s is not identical to base muscle %s.', ...
            muscleNames{m}, muscleNames{1});
    end
    if ~isequal(MM, M)
        error('fit_adaptive_lasso_phase_model:MaskMismatch', ...
            'Mask matrix for %s is not identical to base muscle %s.', ...
            muscleNames{m}, muscleNames{1});
    end
    if ~isequal(sort(keys(featureToIndexMapM)), sort(keys(featureToIndexMap)))
        error('fit_adaptive_lasso_phase_model:FeatureMapMismatch', ...
            'Feature map for %s differs from base muscle %s.', ...
            muscleNames{m}, muscleNames{1});
    end
    [yM, XM] = preprocessExcitationAndFeatures(yAllM, XAllM, ...
        delaySteps, excitationThreshold);
    if ~isequal(size(XM), size(X)) || max(abs(XM(:) - X(:))) > 1e-12
        error('fit_adaptive_lasso_phase_model:PreprocessMismatch', ...
            'Preprocessed feature matrix for %s differs from base muscle %s.', ...
            muscleNames{m}, muscleNames{1});
    end
    Y(:, m) = yM;
end

% 1c. mask 后的设计矩阵 A = X*M
% -------------------------------------------------------------------------
% M 目前是 identity（buildMask），但保留矩阵乘法的形式，方便以后加入
% 结构化 mask（例如只允许某些 feature 进入某些肌肉）。
% 标准化用 population std（std(A,1,1)），与旧 main 一致；系数最后会
% 反变换回原始特征尺度，因此标准化只影响优化条件数，不改变预测。
A = X * M;
A_orig = A;
muA = mean(A_orig, 1);
sA = std(A_orig, 1, 1);
sA(sA == 0) = 1;
A = (A_orig - muA) ./ sA;
p = size(A, 2);

nMus = (p - 8) / 2;
if nMus < 1 || abs(nMus - round(nMus)) > 1e-12
    error('fit_adaptive_lasso_phase_model:InvalidPredictorCount', ...
        'size(A,2) must satisfy p = 2*#mus + 8. Got p = %d.', p);
end

%% ===== STEP 2: OPTIMIZATION SETTINGS =====
optsCommon = struct();
optsCommon.maxIter = 20000;
optsCommon.tol = 1e-10;
optsCommon.verbose = false;
optsCommon.doWarmStart = false;
optsCommon.backtrackBeta = 0.5;
optsCommon.L0 = 1000;
optsCommon.tolCensor = 1e-12;

%% ===== STEP 2B: MAP PHASE INDICES TO DELAYED DATA =====
% preprocessExcitationAndFeatures 把 target 前移 delaySteps、features 后移，
% 所以原来的 phase.stanceIdx/phase.swingIdx（基于未延迟的 .sto）必须减去
% delaySteps 才能和 Y/A 的行号对齐。差的样本直接丢弃。
% 这一步如果和 sampleWeight 的裁剪不一致，会出现“权重和样本错位”的
% 隐蔽 bug，因此两者都使用同一个 delaySteps。
nSamples = size(A, 1);
rawStanceIdx = phase.stanceIdx(:);
rawSwingIdx = phase.swingIdx(:);

stanceIdx = rawStanceIdx(rawStanceIdx > delaySteps) - delaySteps;
swingIdx = rawSwingIdx(rawSwingIdx > delaySteps) - delaySteps;
stanceIdx = stanceIdx(stanceIdx >= 1 & stanceIdx <= nSamples);
swingIdx = swingIdx(swingIdx >= 1 & swingIdx <= nSamples);

stanceMask = false(nSamples, 1);
stanceMask(stanceIdx) = true;
swingMask = false(nSamples, 1);
swingMask(swingIdx) = true;

if any(stanceMask & swingMask)
    error('fit_adaptive_lasso_phase_model:PhaseOverlap', ...
        'Stance and swing indices overlap after delay mapping.');
end
if any(~(stanceMask | swingMask))
    warning('fit_adaptive_lasso_phase_model:PhaseGap', ...
        'Phase mapping leaves %d samples unlabeled after delay.', ...
        nnz(~(stanceMask | swingMask)));
end

phaseNames = {'stance', 'swing'};
phaseIdxMap = struct('stance', stanceIdx, 'swing', swingIdx);

% ---- sample weights ------------------------------------------------------
% 接口给的是“延迟前”的 nT 长度；LASSO 的 Y/A 是 nT-delay 长度。
% 如果长度正好是 nSamples+delaySteps，就裁掉前 delaySteps 个，保证
% 权重和 target 逐行对应。也允许调用者直接传入已裁剪的 nSamples 长度。
wAll = opts.sampleWeight;
if isempty(wAll)
    wAll = ones(nSamples, 1);
else
    wAll = double(wAll(:));
    % The interface stores weights on the full (pre-delay) time grid.  The
    % LASSO response/features are trimmed by delaySteps inside
    % preprocessExcitationAndFeatures, so trim here identically.
    if numel(wAll) == nSamples + delaySteps
        wAll = wAll((delaySteps + 1):end);
    end
    if numel(wAll) ~= nSamples
        error('fit_adaptive_lasso_phase_model:BadWeights', ...
            'options.sampleWeight has %d entries, expected %d.', ...
            numel(wAll), nSamples);
    end
    if any(~isfinite(wAll)) || any(wAll < 0)
        error('fit_adaptive_lasso_phase_model:BadWeights', ...
            'options.sampleWeight must be finite and nonnegative.');
    end
end

%% ===== PRE-ALLOCATIONS =====
beta_OLS = struct();
b_OLS = struct();
W = struct();
Beta_adaptive = struct();
Intercept_adaptive = struct();
YPred_adaptive = struct();
RMSE_adaptive = struct();
stats_ada = struct();
S = struct();
LambdaPath = struct();
Beta_refit = struct();
Intercept_refit = struct();
YPred_refit = struct();
RMSE_refit = struct();
stats_refit = struct();

for phaseIter = 1:numel(phaseNames)
    phaseName = phaseNames{phaseIter};
    beta_OLS.(phaseName) = zeros(p, numTargets);
    b_OLS.(phaseName) = zeros(numTargets, 1);
    W.(phaseName) = zeros(p, numTargets);

    Beta_adaptive.(phaseName) = zeros(p, numTargets);
    Intercept_adaptive.(phaseName) = zeros(numTargets, 1);
    YPred_adaptive.(phaseName) = cell(numTargets, 1);
    RMSE_adaptive.(phaseName) = zeros(numTargets, 1);
    stats_ada.(phaseName) = cell(numTargets, 1);
    S.(phaseName) = cell(numTargets, 1);
    LambdaPath.(phaseName) = cell(numTargets, 1);

    Beta_refit.(phaseName) = zeros(p, numTargets);
    Intercept_refit.(phaseName) = zeros(numTargets, 1);
    YPred_refit.(phaseName) = cell(numTargets, 1);
    RMSE_refit.(phaseName) = zeros(numTargets, 1);
    stats_refit.(phaseName) = cell(numTargets, 1);
end

YPred_adaptive_full = nan(nSamples, numTargets);
YPred_refit_full = nan(nSamples, numTargets);

%% ===== STEP 3-6: OLS -> ADAPTIVE LASSO -> ROBUST REFIT =====
% 外层：每块肌肉独立拟合（共享 X/A，但 target/系数不同）。
% 内层：每个相位独立（stance/swing 的反射律不同）。
% 注意：不同肌肉的特征矩阵相同，但 censoring 的 U/C 划分不同，因为
% 每块肌肉的 y 不同；因此不能把 9 块肌肉合并成一个大回归。
for m = 1:numTargets
    fprintf('\n===== Fitting target muscle: %s (%d/%d) =====\n', ...
        muscleNames{m}, m, numTargets);

    for phaseIter = 1:numel(phaseNames)
        phaseName = phaseNames{phaseIter};
        idx = phaseIdxMap.(phaseName);
        A_phase = A(idx, :);
        A_phase_orig = A_orig(idx, :);
        Y_phase = Y(idx, m);
        wPhase = wAll(idx);
        if sum(wPhase) <= 0
            error('fit_adaptive_lasso_phase_model:NoValidSamples', ...
                'All samples are excluded for %s (%s).', ...
                muscleNames{m}, phaseName);
        end

        % --- STEP 3: INITIAL WEIGHTED OLS ESTIMATE ---
        % 只用“uncensored”样本（y > c+tol）做初始回归：
        %   * 这些样本的 y 是真实观测，不是被 floor 截断的；
        %   * 用它们估计 beta_OLS，再用 |beta_OLS| 构造 adaptive 权重。
        %
        % 初始 OLS 默认也做 Huber IRLS（robustInitIter=2）：
        %   * 这里只跑一次/每肌肉相位，计算量远小于 lambda path；
        %   * 如果初始 OLS 被 burst 拉偏，adaptive 权重会跟着错，
        %     所以先 robust 一下更稳；
        %   * robustInitIter=0 时退回普通加权 OLS。
        c_m = c;
        tolC = optsCommon.tolCensor;
        Uidx = find(Y_phase > c_m + tolC);
        nU = numel(Uidx);

        if nU == 0
            warning('fit_adaptive_lasso_phase_model:NoUncensored', ...
                'No uncensored samples for %s (%s). Using zero OLS.', ...
                muscleNames{m}, phaseName);
            beta_OLS_std = zeros(p, 1);
            b_OLS_std = c_m;
        else
            A_U = A_phase(Uidx, :);
            y_U = Y_phase(Uidx);
            w_U = max(wPhase(Uidx), 0);
            augA = [A_U, ones(nU, 1)];
            rankDef = rank(augA) < min(size(augA));
            if rankDef
                warning('fit_adaptive_lasso_phase_model:RankDeficient', ...
                    'OLS for %s (%s) is rank-deficient (nU=%d, p=%d). Using pinv.', ...
                    muscleNames{m}, phaseName, nU, p);
            end

            % 第一步：普通加权最小二乘
            coeffs = solve_weighted_ls(augA, y_U, w_U, rankDef);

            % 第二步：Huber IRLS（只对 uncensored 样本；这些样本没有 censoring）
            if opts.robustInitIter > 0 && nU >= 5 && sum(w_U) > 0
                for it = 1:opts.robustInitIter
                    res = y_U - augA * coeffs;
                    fitH = w_U > 0;
                    if ~any(fitH)
                        break;
                    end
                    rH = res(fitH);
                    s = 1.4826 * median(abs(rH - median(rH)));
                    if ~isfinite(s) || s < 1e-8
                        s = 1e-8;
                    end
                    delta = opts.huberDeltaFactor * s;
                    wH = ones(nU, 1);
                    wH(fitH) = min(1, delta ./ max(abs(rH), eps));
                    coeffs = solve_weighted_ls(augA, y_U, w_U .* wH, rankDef);
                end
            end

            beta_OLS_std = coeffs(1:p);
            b_OLS_std = coeffs(end);
        end

        beta_OLS_orig = beta_OLS_std ./ sA.';
        b_OLS_orig = b_OLS_std - muA * beta_OLS_orig;
        beta_OLS.(phaseName)(:, m) = beta_OLS_orig;
        b_OLS.(phaseName)(m) = b_OLS_orig;

        % --- STEP 4: ADAPTIVE WEIGHTS ---
        % Zou (2006) adaptive LASSO：对 OLS 系数大的特征给较小惩罚，
        % 对 OLS 系数接近 0 的特征给很大惩罚。epsilon_ols 防止除零。
        % 这里用标准化空间的 beta_OLS_std；后续 A* = A ./ w 也在这个空间，
        % 因此 w 的量级和 A 的列尺度一致。
        w = 1 ./ max(abs(beta_OLS_std), epsilon_ols).^gamma;
        if any(~isfinite(w))
            error('fit_adaptive_lasso_phase_model:BadWeights', ...
                'Non-finite adaptive weights for %s (%s).', ...
                muscleNames{m}, phaseName);
        end
        W.(phaseName)(:, m) = w;

        % --- STEP 5: ADAPTIVE LASSO VIA sgl_fit (LAMBDA PATH) ---
        % 列缩放：A_star = A_phase .* (1./w)'。
        % 这样对 beta_star 的普通 L1 惩罚等价于对 beta_orig 的加权 L1
        %   lambda * sum_j w_j |beta_j|。
        % lambda2Max 是“beta=0 仍最优”的最小 lambda（censored 版本），
        % 用它作为 path 起点，再按 alphaGrid 从 1 降到 1e-4。
        A_star = A_phase .* (1 ./ w)';
        [lambda2Max, b0AtZero, lambdaMaxInfo] = ...
            compute_lambda2_max_censored(A_star, Y_phase, c_m, rho, tolC);

        if lambda2Max > 0
            alphaGridThis = lambdaPathConfig.alphaGrid;
            lambdaGridThis = lambda2Max * alphaGridThis;
        else
            alphaGridThis = 1;
            lambdaGridThis = 0;
        end

        nPath = numel(lambdaGridThis);
        betaStarStdPath = zeros(p, nPath);
        betaAdaStdPath = zeros(p, nPath);
        interceptAdaStdPath = zeros(1, nPath);
        betaRefitStdPath = zeros(p, nPath);
        interceptRefitStdPath = zeros(1, nPath);
        trainErrorAdaptive = nan(nPath, 1);
        trainErrorRefit = nan(nPath, 1);
        supportSize = zeros(nPath, 1);
        supportPath = cell(nPath, 1);
        statsAdaPath = cell(nPath, 1);
        statsRefitPath = cell(nPath, 1);

        prevBetaStarStd = zeros(p, 1);
        prevInterceptStd = b0AtZero;

        % 沿 lambda path 从强惩罚到弱惩罚：
        %   * 每个点都调用 sgl_fit 解 adaptive LASSO；
        %   * 用上一点的解 warm start，显著减少 FISTA 迭代；
        %   * 每个点还在当前 support 上做一次 refit（去偏），用于选 lambda。
        for kk = 1:nPath
            lambda2_k = lambdaGridThis(kk);
            if lambda2Max <= 0
                betaStarStdK = zeros(p, 1);
                interceptAdaStdK = b0AtZero;
                statsAdaK = struct('exitType', 'degenerate_lambda2max_zero', ...
                    'lambda1', lambda1, 'lambda2', lambda2_k, ...
                    'converged', true, 'iters', 0);
            else
                optsM = optsCommon;
                optsM.muscleName = [muscleNames{m}, '_', phaseName, ...
                    '_path_', num2str(kk)];
                optsM.c = c_m;
                optsM.rho = rho;
                optsM.beta0 = prevBetaStarStd;
                optsM.b0 = prevInterceptStd;
                optsM.doWarmStart = true;
                optsM.sampleWeight = wPhase;
                [betaStarStdK, interceptAdaStdK, ~, statsAdaK] = ...
                    sgl_fit(A_star, Y_phase, lambda1, lambda2_k, optsM);
            end

            betaAdaStdK = betaStarStdK ./ w;
            S_k = find(abs(betaAdaStdK) > nzTol);

            if isempty(S_k)
                betaRefitStdK = zeros(p, 1);
                interceptRefitStdK = b0AtZero;
                statsRefitK = struct('exitType', 'intercept_only_empty_support', ...
                    'lambda1', 0, 'lambda2', 0, 'converged', true, 'iters', 0);
            else
                A_refit = A_phase;
                nonS = setdiff(1:p, S_k);
                A_refit(:, nonS) = 0;
                beta0Refit = betaAdaStdK;
                beta0Refit(nonS) = 0;

                optsRefit = optsCommon;
                optsRefit.muscleName = [muscleNames{m}, '_', phaseName, ...
                    '_refit_path_', num2str(kk)];
                optsRefit.c = c_m;
                optsRefit.rho = rho;
                optsRefit.doWarmStart = true;
                optsRefit.beta0 = beta0Refit;
                optsRefit.b0 = interceptAdaStdK;
                optsRefit.sampleWeight = wPhase;
                [betaRefitStdK, interceptRefitStdK, ~, statsRefitK] = ...
                    sgl_fit(A_refit, Y_phase, 0, 0, optsRefit);
                betaRefitStdK(nonS) = 0;
            end

            yPredAdaK = max(A_phase * betaAdaStdK + interceptAdaStdK, c_m);
            trainErrorAdaptive(kk) = weighted_rmse(Y_phase, yPredAdaK, wPhase);

            yPredRefitK = max(A_phase * betaRefitStdK + interceptRefitStdK, c_m);
            trainErrorRefit(kk) = weighted_rmse(Y_phase, yPredRefitK, wPhase);

            betaStarStdPath(:, kk) = betaStarStdK;
            betaAdaStdPath(:, kk) = betaAdaStdK;
            interceptAdaStdPath(kk) = interceptAdaStdK;
            betaRefitStdPath(:, kk) = betaRefitStdK;
            interceptRefitStdPath(kk) = interceptRefitStdK;
            supportPath{kk} = S_k;
            supportSize(kk) = numel(S_k);
            % Store compact diagnostics only: the full FISTA objHist of every
            % lambda path would make the saved fit diagnostics unnecessarily
            % large (hundreds of MB for 18 muscles x 2 phases x 60 lambdas).
            statsAdaPath{kk} = compact_stats(statsAdaK);
            statsRefitPath{kk} = compact_stats(statsRefitK);

            prevBetaStarStd = betaStarStdK;
            prevInterceptStd = interceptAdaStdK;
        end

        % 选择规则（与旧 main 一致）：
        %   1) 找 refit RMSE 最小值 minErr；
        %   2) 允许 (1+tolSelect)*minErr 的容差；
        %   3) 在容差内取“最稀疏”的 lambda（path 从大到小，first feasible）。
        % 这样不会为了极小 RMSE 牺牲稀疏性；tolSelect=0.05 表示允许 5% 误差。
        validErr = isfinite(trainErrorRefit);
        if ~any(validErr)
            error('fit_adaptive_lasso_phase_model:NoValidPathRMSE', ...
                'No finite refit errors for %s (%s).', muscleNames{m}, phaseName);
        end
        [minErr, localMinPos] = min(trainErrorRefit(validErr));
        validIdx = find(validErr);
        minIdx = validIdx(localMinPos);
        thresholdErr = (1 + lambdaPathConfig.tolSelect) * minErr;
        feasibleIdx = find(validErr & ...
            trainErrorRefit <= thresholdErr + 10 * eps(max(1, thresholdErr)));
        if isempty(feasibleIdx)
            selectedIdx = minIdx;
        else
            selectedIdx = feasibleIdx(1);
        end

        beta_star_std = betaStarStdPath(:, selectedIdx);
        interceptAdaStd = interceptAdaStdPath(selectedIdx);
        betaAdaStd = betaAdaStdPath(:, selectedIdx);
        betaRefitStd = betaRefitStdPath(:, selectedIdx);
        interceptRefitStd = interceptRefitStdPath(selectedIdx);
        S.(phaseName){m} = supportPath{selectedIdx};
        stats_ada.(phaseName){m} = statsAdaPath{selectedIdx};
        stats_refit.(phaseName){m} = statsRefitPath{selectedIdx};

        % --- STEP 6: ROBUST HUBER IRLS ON THE SELECTED SUPPORT ---
        % 先用 path 里选定的 support 和系数作为初值；然后迭代：
        %   1) 用当前权重做一次 support 上的 WLS（lambda=0）；
        %   2) 只在 uncensored 样本上计算残差，用 1.4826*MAD 估计尺度；
        %   3) Huber 权重 wH = min(1, delta/|r|)，delta=1.345*scale；
        %   4) 组合 base weight（接缝/无效样本为 0）得到新的 wCur。
        %
        % 为什么只对 uncensored 样本降权：
        %   censored 样本的 y_i=c，若 latent prediction A_i beta+b > c，
        %   这是要惩罚的 hinge 违反，不是“y 的离群点”。如果也按
        %   |y-max(pred,c)| 给它降权，会削弱 rho 对 censoring 的约束。
        %   因此 censored 样本保持 base weight（通常为 1），只对
        %   uncensored 样本做 Huber 降权。
        S_k = supportPath{selectedIdx};
        nonS = setdiff(1:p, S_k);
        A_refit = A_phase;
        A_refit(:, nonS) = 0;
        wCur = wPhase;
        betaRefitStdK = betaRefitStd(:, selectedIdx);
        interceptRefitStdK = interceptRefitStd(selectedIdx);
        statsRefitFinal = statsRefitPath{selectedIdx};

        for it = 1:opts.robustIter
            optsRefit = optsCommon;
            optsRefit.muscleName = [muscleNames{m}, '_', phaseName, ...
                '_robust_', num2str(it)];
            optsRefit.c = c_m;
            optsRefit.rho = rho;
            optsRefit.doWarmStart = true;
            optsRefit.beta0 = betaRefitStdK;
            optsRefit.b0 = interceptRefitStdK;
            optsRefit.sampleWeight = max(wCur, 0);

            [betaRefitStdK, interceptRefitStdK, ~, statsRefitFinal] = ...
                sgl_fit(A_refit, Y_phase, 0, 0, optsRefit);
            betaRefitStdK(nonS) = 0;

            latent = A_phase * betaRefitStdK + interceptRefitStdK;
            isCensored = Y_phase <= c_m + optsCommon.tolCensor;
            fitU = ~isCensored & wPhase > 0;

            if any(fitU)
                % 只用“uncensored 且 base weight>0”的样本估计残差尺度；
                % 被 sampleWeight=0 排除的接缝样本不应影响 Huber 尺度。
                resU = Y_phase(fitU) - latent(fitU);
                s = 1.4826 * median(abs(resU - median(resU)));
                if ~isfinite(s) || s < 1e-8
                    s = 1e-8;
                end
                delta = opts.huberDeltaFactor * s;

                % 默认所有样本 wH=1；只对可用的 uncensored 样本降权。
                wH = ones(size(Y_phase));
                wH(fitU) = min(1, delta ./ max(abs(resU), eps));
                wCur = wPhase .* wH;
            else
                % 没有可用的 uncensored 样本时没有可靠残差尺度，停止 robust。
                wCur = wPhase;
            end
        end
        stats_refit.(phaseName){m} = compact_stats(statsRefitFinal);

        % pathInfo 保留完整路径（用于诊断/画 lambda path），但 stats 已
        % compact，避免把每次 FISTA 的 20k 步 objHist 都写进 .mat。
        pathInfo = struct();
        pathInfo.alphaGrid = alphaGridThis;
        pathInfo.lambdaGrid = lambdaGridThis;
        pathInfo.lambda2Max = lambda2Max;
        pathInfo.b0AtZero = b0AtZero;
        pathInfo.lambdaMaxInfo = lambdaMaxInfo;
        pathInfo.trainErrorAdaptive = trainErrorAdaptive;
        pathInfo.trainErrorRefit = trainErrorRefit;
        pathInfo.supportSize = supportSize;
        pathInfo.supportPath = supportPath;
        pathInfo.selectedIdx = selectedIdx;
        pathInfo.selectedAlpha = alphaGridThis(selectedIdx);
        pathInfo.selectedLambda2 = lambdaGridThis(selectedIdx);
        pathInfo.selectedSupport = supportPath{selectedIdx};
        pathInfo.selectionMetric = lambdaPathConfig.selectionMetric;
        pathInfo.tolSelect = lambdaPathConfig.tolSelect;
        pathInfo.minTrainErrorRefit = minErr;
        pathInfo.thresholdTrainErrorRefit = thresholdErr;
        pathInfo.selectedTrainErrorRefit = trainErrorRefit(selectedIdx);
        pathInfo.selectedTrainErrorAdaptive = trainErrorAdaptive(selectedIdx);
        pathInfo.betaStarStdPath = betaStarStdPath;
        pathInfo.betaAdaStdPath = betaAdaStdPath;
        pathInfo.interceptAdaStdPath = interceptAdaStdPath;
        pathInfo.betaRefitStdPath = betaRefitStdPath;
        pathInfo.interceptRefitStdPath = interceptRefitStdPath;
        pathInfo.statsAdaPath = statsAdaPath;
        pathInfo.statsRefitPath = statsRefitPath;
        pathInfo.sampleWeight = wPhase;
        pathInfo.robustIter = opts.robustIter;
        pathInfo.robustFinalWeights = wCur;
        LambdaPath.(phaseName){m} = pathInfo;

        fprintf(['[%s|%s] lambda2Max=%.6g | selectedIdx=%d/%d | ' ...
            'selectedLambda2=%.6g | minRefitRMSE=%.6g | ' ...
            'selectedRefitRMSE=%.6g | nnz=%d\n'], ...
            muscleNames{m}, phaseName, lambda2Max, selectedIdx, nPath, ...
            lambdaGridThis(selectedIdx), minErr, ...
            weighted_rmse(Y_phase, max(A_phase * betaRefitStdK + ...
                interceptRefitStdK, c_m), wPhase), numel(S_k));

        betaAdaOrig = betaAdaStd ./ sA.';
        interceptAdaOrig = interceptAdaStd - muA * betaAdaOrig;
        betaRefitOrig = betaRefitStdK ./ sA.';
        interceptRefitOrig = interceptRefitStdK - muA * betaRefitOrig;

        Beta_adaptive.(phaseName)(:, m) = betaAdaOrig;
        Intercept_adaptive.(phaseName)(m) = interceptAdaOrig;
        Beta_refit.(phaseName)(:, m) = betaRefitOrig;
        Intercept_refit.(phaseName)(m) = interceptRefitOrig;

        yPredAda = max(A_phase_orig * betaAdaOrig + interceptAdaOrig, c_m);
        YPred_adaptive.(phaseName){m} = yPredAda;
        RMSE_adaptive.(phaseName)(m) = weighted_rmse(Y_phase, yPredAda, wPhase);

        yPredRefit = max(A_phase_orig * betaRefitOrig + interceptRefitOrig, c_m);
        YPred_refit.(phaseName){m} = yPredRefit;
        RMSE_refit.(phaseName)(m) = weighted_rmse(Y_phase, yPredRefit, wPhase);

        YPred_adaptive_full(idx, m) = yPredAda;
        YPred_refit_full(idx, m) = yPredRefit;
    end
end

if any(isnan(YPred_adaptive_full(:))) || any(isnan(YPred_refit_full(:)))
    warning('fit_adaptive_lasso_phase_model:PredictionGap', ...
        'Stitched predictions contain NaN values. Check phase mapping.');
end

RMSE_adaptive_full = sqrt(sum(wAll .* (Y - YPred_adaptive_full).^2, 1) ./ ...
    max(sum(wAll), eps)).';
RMSE_refit_full = sqrt(sum(wAll .* (Y - YPred_refit_full).^2, 1) ./ ...
    max(sum(wAll), eps)).';

%% ===== SUMMARY =====
fprintf('\n===== MULTI-MUSCLE ADAPTIVE LASSO SUMMARY (PHASE-SPECIFIC) =====\n');
for m = 1:numTargets
    for phaseIter = 1:numel(phaseNames)
        phaseName = phaseNames{phaseIter};
        pathInfo = LambdaPath.(phaseName){m};
        fprintf(['[%s|%s] OLS_nnz=%d | Ada_nnz=%d | Refit_nnz=%d | ' ...
            'RMSE_ada=%.4g | RMSE_refit=%.4g | lambda2Sel=%.4g\n'], ...
            muscleNames{m}, phaseName, ...
            nnz(abs(beta_OLS.(phaseName)(:, m)) > nzTol), ...
            nnz(abs(Beta_adaptive.(phaseName)(:, m)) > nzTol), ...
            nnz(abs(Beta_refit.(phaseName)(:, m)) > nzTol), ...
            RMSE_adaptive.(phaseName)(m), RMSE_refit.(phaseName)(m), ...
            pathInfo.selectedLambda2);
    end
end

%% ===== ASSEMBLE OUTPUTS =====
fit = struct();
fit.source       = opts.source;
fit.stoFile      = stoFilePath;
fit.muscleNames  = muscleNames;
fit.A            = A_orig;
fit.A_std        = A;
fit.featureScaling = struct('muA', muA, 'sA', sA);
fit.Y            = Y;
fit.phase        = struct('stanceIdx', stanceIdx, 'swingIdx', swingIdx, ...
    'stanceMask', stanceMask, 'swingMask', swingMask);
fit.sampleWeight = wAll;
fit.robustInitIter = opts.robustInitIter;
fit.robustIter   = opts.robustIter;
fit.beta_OLS     = beta_OLS;
fit.b_OLS        = b_OLS;
fit.W            = W;
fit.Beta_adaptive = Beta_adaptive;
fit.Intercept_adaptive = Intercept_adaptive;
fit.YPred_adaptive = YPred_adaptive_full;
fit.YPred_adaptive_phase = YPred_adaptive;
fit.RMSE_adaptive = RMSE_adaptive;
fit.RMSE_adaptive_full = RMSE_adaptive_full;
fit.stats_ada    = stats_ada;
fit.S            = S;
fit.Beta_refit   = Beta_refit;
fit.Intercept_refit = Intercept_refit;
fit.YPred_refit  = YPred_refit_full;
fit.YPred_refit_phase = YPred_refit;
fit.RMSE_refit   = RMSE_refit;
fit.RMSE_refit_full = RMSE_refit_full;
fit.stats_refit  = stats_refit;
fit.lambdaPathConfig = lambdaPathConfig;
fit.options      = opts;
fit.LambdaPath   = LambdaPath;

indexToFeatureName = buildInverseFeatureMap(featureToIndexMap, size(X, 2));
fit.featureNames = indexToFeatureName;

% ---- grouped controller struct ----
% 为什么要把 <= nzTol 的 beta 强制清零：
%   load_lasso_reflex_controller.m 会检查
%       abs(beta(~mask)) > 0
%   并报错。mask 由 abs(beta)>nzTol 得到；如果不先清零，像 1e-12 这样的
%   数值残渣会落在 mask=false 区域，loader 会拒绝加载。
% 清零点不影响预测（<1e-5），但保证结构零和 mask 严格一致。
% Force an exact structural-zero pattern: load_lasso_reflex_controller.m
% rejects any beta entry that is nonzero where mask is false.
betaStance = fit.Beta_refit.stance;
betaStance(abs(betaStance) <= nzTol) = 0;
biasStance = fit.Intercept_refit.stance(:)';
maskStance = abs(betaStance) > nzTol;

betaSwing  = fit.Beta_refit.swing;
betaSwing(abs(betaSwing) <= nzTol) = 0;
biasSwing  = fit.Intercept_refit.swing(:)';
maskSwing  = abs(betaSwing) > nzTol;

lasso = struct();
lasso.format    = 'grouped';
lasso.nPhases   = 5;
lasso.phaseIds  = 0:4;
lasso.nGroups   = 2;
lasso.groups    = struct();
lasso.groups(1).label  = opts.groupLabels{1};
lasso.groups(1).phases = opts.groupPhases{1};
lasso.groups(1).beta   = betaStance;
lasso.groups(1).bias   = biasStance;
lasso.groups(1).mask   = maskStance;
lasso.groups(2).label  = opts.groupLabels{2};
lasso.groups(2).phases = opts.groupPhases{2};
lasso.groups(2).beta   = betaSwing;
lasso.groups(2).bias   = biasSwing;
lasso.groups(2).mask   = maskSwing;

lasso.featureNames = indexToFeatureName;
lasso.muscleNames  = muscleNames;
lasso.source       = opts.source;
lasso.datasetFile  = stoFilePath;
lasso.created      = char(datetime);
lasso.preprocessing = opts.preprocessing;

fprintf('\n===== LASSO controller assembled (grouped format) =====\n');
fprintf('  stance nnz(beta) = %d, swing nnz(beta) = %d\n', ...
    nnz(maskStance), nnz(maskSwing));
fprintf('  total optimizable parameters = %d\n', ...
    nnz(maskStance) + 9 + nnz(maskSwing) + 9);
fprintf('=======================================================\n');
end

%% ========================================================================
% 局部函数：加权最小二乘（供初始 OLS / Huber IRLS 共用）
% ========================================================================
% 用 row-scaling 实现加权：
%   minimize || sqrt(w) .* (y - augA*coeffs) ||^2
% 等价于 (augA .* sw) \ (y .* sw)，其中 sw=sqrt(w)。
% 若所有权重为 0（极端排除情况），退化为普通 OLS。
function coeffs = solve_weighted_ls(augA, y, w, usePinv)
w = max(w(:), 0);
sw = sqrt(w);
if sum(w) <= 0
    sw = ones(size(w));
end
if usePinv
    coeffs = pinv(augA .* sw) * (y .* sw);
else
    coeffs = (augA .* sw) \ (y .* sw);
end
end

%% ========================================================================
% 局部函数：默认参数与嵌套 lambdaPathConfig 合并
% ========================================================================
% 采用“用户字段优先、缺省/空值才填默认”的策略。
% lambdaPathConfig 是嵌套 struct，单独处理：
%   * 用户给了 nLambda 但没给 alphaGrid 时，按新的 nLambda 重新生成网格；
%   * 用户显式给 alphaGrid 时，完全尊重用户网格（即使 nLambda 不一致）。
function opts = apply_defaults(options)
defaults = struct();
defaults.gamma = 1;
defaults.epsilonOLS = 1e-8;
defaults.nzTol = 1e-5;
defaults.lambda1 = 0;
defaults.lambdaPathConfig = struct( ...
    'nLambda', 60, ...
    'lambdaMinRatio', 1e-4, ...
    'tolSelect', 0.05, ...
    'selectionMetric', 'post_refit_train_censored_rmse');
defaults.lambdaPathConfig.alphaGrid = logspace(0, ...
    log10(defaults.lambdaPathConfig.lambdaMinRatio), ...
    defaults.lambdaPathConfig.nLambda).';
defaults.c = 0.01;
defaults.rho = 50;
defaults.delaySteps = 1;
defaults.excitationThreshold = 0.01;
defaults.sampleWeight = [];
defaults.robustInitIter = 2;   % Huber IRLS for initial OLS (0 disables)
defaults.robustIter = 2;       % Huber IRLS for post-selection refit
defaults.huberDeltaFactor = 1.345;
defaults.verbose = true;
defaults.groupLabels = {'stance', 'swing'};
defaults.groupPhases = {[0 1], [2 3 4]};
defaults.source = 'MocoInverse';
defaults.preprocessing = struct();

fn = fieldnames(defaults);
opts = options;
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}) || isempty(opts.(fn{i}))
        opts.(fn{i}) = defaults.(fn{i});
    end
end

% Nested lambdaPathConfig: fill missing sub-fields, keep user overrides.
lpc = defaults.lambdaPathConfig;
userHasAlphaGrid = false;
if isfield(options, 'lambdaPathConfig') && isstruct(options.lambdaPathConfig)
    u = options.lambdaPathConfig;
    if isfield(u, 'alphaGrid') && ~isempty(u.alphaGrid)
        userHasAlphaGrid = true;
    end
    fn2 = fieldnames(u);
    for i = 1:numel(fn2)
        lpc.(fn2{i}) = u.(fn2{i});
    end
end
if ~userHasAlphaGrid
    lpc.alphaGrid = logspace(0, log10(lpc.lambdaMinRatio), lpc.nLambda).';
end
opts.lambdaPathConfig = lpc;
end

%% ========================================================================
% 局部函数：加权 RMSE
% ========================================================================
% 用于 lambda 选择和报告。w=0 的样本（接缝/无效）不参与误差；
% 如果所有 w 都为 0（理论上不该发生），退化为普通 RMSE 以免 NaN。
function r = weighted_rmse(y, yhat, w)
w = max(w(:), 0);
if sum(w) <= 0
    r = sqrt(mean((y(:) - yhat(:)).^2));
else
    r = sqrt(sum(w .* (y(:) - yhat(:)).^2) / sum(w));
end
end

%% ========================================================================
% 局部函数：压缩 sgl_fit 的 stats
% ========================================================================
% sgl_fit 的 stats 含 objHist/stepHist 等长数组；对 18 肌肉×2 相位×60
% lambda，完整保存会有数百 MB。这里只保留诊断真正需要的标量/短字段：
%   exitType, iters, converged, nnz, lambda1/2, c, rho, U/C count,
%   finalObjective（objHist(end)）。
% 完整路径信息（beta path）仍在 LambdaPath 里，足够复现/画图。
function s = compact_stats(st)
%COMPACT_STATS  Keep only the small, useful part of an sgl_fit stats struct.
s = struct();
fields = {'exitType', 'iters', 'converged', 'nnz', 'lambda1', 'lambda2', ...
          'c', 'rho', 'Ucount', 'Ccount'};
for i = 1:numel(fields)
    if isfield(st, fields{i})
        s.(fields{i}) = st.(fields{i});
    end
end
if isfield(st, 'objHist') && ~isempty(st.objHist)
    s.finalObjective = st.objHist(end);
end
end

%% ========================================================================
% 局部函数：feature 名 -> 行号 的反向映射
% ========================================================================
% extractStoMuscleFeatures 返回 containers.Map(featureName -> column)。
% lasso.featureNames 需要按行号排列的 cellstr，便于人和调试工具核对
% beta 的每一行到底对应哪个特征。这里做严格的双射检查。
function indexToFeatureName = buildInverseFeatureMap(featureToIndexMap, numFeatures)
ks = keys(featureToIndexMap);
vs = cell2mat(values(featureToIndexMap));
if numel(ks) ~= numFeatures
    error('buildInverseFeatureMap:CountMismatch', ...
        'feature map size (%d) does not match numFeatures (%d).', ...
        numel(ks), numFeatures);
end
if any(vs < 1) || any(vs > numFeatures) || any(abs(vs - round(vs)) > 0)
    error('buildInverseFeatureMap:BadIndex', ...
        'feature indices must be integers in [1, %d].', numFeatures);
end
if numel(unique(vs)) ~= numel(vs)
    error('buildInverseFeatureMap:NonBijective', 'feature map is not one-to-one.');
end
indexToFeatureName = cell(numFeatures, 1);
for t = 1:numel(ks)
    indexToFeatureName{vs(t)} = ks{t};
end
end
