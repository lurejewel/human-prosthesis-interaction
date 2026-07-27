% Adaptive LASSO runner script for phase-dependent controlled muscles.
% 
% P.S. While the stance-swing split may work well for normal walking,
% dividing the gait cycle into additional phases (thus introducing more
% segmented controllers) may provide a better fit. Properly defining and
% segmenting gait phases is therefore an important reserach consideration.
% ========================================================================
% PURPOSE:
%   Fit a lower-censored adaptive LASSO model to predict muscle excitation
%   from biomechanical features using sgl_fit as a black-box backend.
%
% WORKFLOW:
%   1. Load and preprocess data
%   2. Compute OLS initial estimate on uncensored samples
%   3. Derive adaptive weights w_j = 1/max(|beta_OLS_j|, eps)^gamma
%   4. Column-scale A to convert adaptive LASSO to standard LASSO
%   5. Solve via sgl_fit (lambda1 = 0)
%   6. Recover coefficients in original space
%   7. Post-selection OLS on support set via sgl_fit (lambda = 0)
%   8. Report and visualize
%
% REFERENCE:
%   Zou, H. (2006). "The Adaptive Lasso and Its Oracle Properties."
%   Journal of the American Statistical Association, 101(476), 1418-1429.
% 
% TODO:
% - 所有肌肉成功后的行走正动力学验证
% - 基于多人行走实测数据的肌肉控制律预测和分析
% - 其他具备现有控制律的动作（跑步？）的验证
% - 其他动作（跑步、上下楼梯、蹲起、游泳等周期性动作，以及跳跃等非周期性动作，现有预测动力学无法做到的）的仿真及实验数据下的预测动力学验证
% - 病态动作的验证
% 
% IMPORTANT:
% - maskMatrix需要根据肌肉特性赋值（重读经典论文，导入先验证据）；是否还需要sparse group？
% - 是否需要为不同动作/相位设置不同的lambda、c、rho？
% - (正动力学完成后）代理模型研究：对应Annuals Review of Biomedical Engineering中的大脑预测模型）
% - gamma & lambda的二维网格搜索：cross-validation for grid selection (one-SE rule)
% - 消融实验：pure lasso vs adaptive lasso vs adaptive lasso with priori vs
% something else
% - c, rho，以及其他参数的灵敏度分析
% - 按照source对象分组：[TA], [SOL], ..., [hip], [knee], [ankle]，而非按照反馈类型分组
% - 后续考虑按照muscular function分组，但会有重叠的问题（因为有跨双关节的肌肉），要想想怎么解决
% - stability selection, 解决source的可信度问题（重采样、保留selection probability高于阈值的source）
% ========================================================================

close all; clear; clc;

scriptFullPath = mfilename('fullpath');
if isempty(scriptFullPath)
    logFilePath = fullfile(pwd, 'Adaptive_LASSO_StanceSwing.log');
else
    [scriptDir, scriptBase] = fileparts(scriptFullPath);
    logFilePath = fullfile(scriptDir, [scriptBase, '.log']);
end
if exist(logFilePath, 'file')
    delete(logFilePath);
end
diary(logFilePath);
diary on;
logCleanup = onCleanup(@() diary('off'));
fprintf('Logging to %s\n', logFilePath);

%% ===== STEP 1: DATA LOADING & PREPROCESSING =====
% Load target muscle excitations and one shared biomechanical feature matrix.
stoFilePath = 'UN.sto';

% ---- auto-detect right-leg muscle names from UN.sto header ----
fid = fopen(stoFilePath, 'r');
if fid == -1
    error('Cannot open STO file: %s', stoFilePath);
end
ln = '';
while ischar(ln) && ~contains(ln, 'endheader')
    ln = fgetl(fid);
end
if ~ischar(ln)
    fclose(fid);
    error('endheader not found in STO file: %s', stoFilePath);
end
headerLine = fgetl(fid);
fclose(fid);
colNames = strsplit(strtrim(headerLine), '\t');

% find columns whose names end with '_r.activation' (right-leg muscles)
isRightActivation = endsWith(colNames, '_r.activation');
muscleIdx = find(isRightActivation);
muscleNames = cell(1, numel(muscleIdx));
for i = 1:numel(muscleIdx)
    % e.g. 'hamstrings_r.activation' -> 'hamstrings_r'
    muscleNames{i} = extractBefore(colNames{muscleIdx(i)}, '.activation');
end
fprintf('Auto-detected %d right-leg muscles from %s:\n', numel(muscleNames), stoFilePath);
fprintf('  %s\n', muscleNames{:});
numTargets = numel(muscleNames);

% Tunable parameters (adaptive LASSO + censoring)
gamma = 1;
epsilon_ols = 1e-8;
nzTol = 1e-5;
lambda1 = 0;
lambdaPathConfig = struct();
lambdaPathConfig.nLambda = 60;
lambdaPathConfig.lambdaMinRatio = 1e-4;
lambdaPathConfig.alphaGrid = logspace(0, log10(lambdaPathConfig.lambdaMinRatio), lambdaPathConfig.nLambda).';
lambdaPathConfig.tolSelect = 0.05;
lambdaPathConfig.selectionMetric = 'post_refit_train_censored_rmse';
c = 0.01;
rho = 50;
delaySteps = 1;
excitationThreshold = 0.01;

assert(lambda1 == 0, ...
    'Adaptive LASSO column-scaling only adapts L1 penalty. lambda1 must be 0.');

fprintf('Lambda path config: nLambda=%d, alpha from %.3g to %.3g, tolSelect=%.3g, metric=%s\n', ...
    lambdaPathConfig.nLambda, lambdaPathConfig.alphaGrid(1), lambdaPathConfig.alphaGrid(end), ...
    lambdaPathConfig.tolSelect, lambdaPathConfig.selectionMetric);

% Extract base muscle to initialize shared X and M.
[yAllBase, XAll, featureToIndexMap, M, phase] = extractStoMuscleFeatures(stoFilePath, muscleNames{1});
[yBase, X] = preprocessExcitationAndFeatures(yAllBase, XAll, delaySteps, excitationThreshold);

Y = zeros(numel(yBase), numTargets);
Y(:, 1) = yBase;

% Load remaining muscles and verify shared feature/mask assumptions.
for m = 2:numTargets
    [yAllM, XAllM, featureToIndexMapM, MM, ~] = extractStoMuscleFeatures(stoFilePath, muscleNames{m});

    if ~isequal(size(XAllM), size(XAll)) || max(abs(XAllM(:) - XAll(:))) > 1e-12
        error('Adaptive_LASSO:FeatureMismatch', ...
            'Feature matrix for %s is not identical to base muscle %s.', muscleNames{m}, muscleNames{1});
    end

    if ~isequal(MM, M)
        error('Adaptive_LASSO:MaskMismatch', ...
            'Mask matrix for %s is not identical to base muscle %s.', muscleNames{m}, muscleNames{1});
    end

    if ~isequal(sort(keys(featureToIndexMapM)), sort(keys(featureToIndexMap)))
        error('Adaptive_LASSO:FeatureMapMismatch', ...
            'Feature map for %s differs from base muscle %s.', muscleNames{m}, muscleNames{1});
    end

    [yM, XM] = preprocessExcitationAndFeatures(yAllM, XAllM, delaySteps, excitationThreshold);
    if ~isequal(size(XM), size(X)) || max(abs(XM(:) - X(:))) > 1e-12
        error('Adaptive_LASSO:PreprocessMismatch', ...
            'Preprocessed feature matrix for %s is not identical to base muscle %s.', muscleNames{m}, muscleNames{1});
    end

    Y(:, m) = yM;
end

fprintf('Loaded %d muscles with shared feature matrix: X size = [%d, %d], Y size = [%d, %d].\n', ...
    numTargets, size(X, 1), size(X, 2), size(Y, 1), size(Y, 2));

% Build one shared A = X * M and solve each muscle independently.
A = X * M; % feature matrix X (augmented), masked with the sparse matrix M

% Center and scale features before OLS/optimization (population std).
A_orig = A;
muA = mean(A_orig, 1);
sA = std(A_orig, 1, 1);
sA(sA == 0) = 1;
A = (A_orig - muA) ./ sA;
p = size(A, 2);

nMus = (p - 8) / 2;
if nMus < 1 || abs(nMus - round(nMus)) > 1e-12
    error('Adaptive_LASSO:InvalidPredictorCount', ...
        'size(A,2) must satisfy p = 2*#mus + 8. Got p = %d.', p);
end

%% ===== STEP 2: OPTIMIZATION SETTINGS =====
optsCommon = struct();
optsCommon.maxIter = 20000;
optsCommon.tol = 1e-10;
optsCommon.verbose = false;
optsCommon.doWarmStart = false; % true;
optsCommon.backtrackBeta = 0.5;
optsCommon.L0 = 1000;
optsCommon.tolCensor = 1e-12;

%% ===== STEP 2B: MAP PHASE INDICES TO DELAYED DATA =====
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
    error('Adaptive_LASSO:PhaseOverlap', ...
        'Stance and swing indices overlap after delay mapping.');
end

if any(~(stanceMask | swingMask))
    warning('Adaptive_LASSO:PhaseGap', ...
        'Phase mapping leaves %d samples unlabeled after delay.', nnz(~(stanceMask | swingMask)));
end

phaseNames = {'stance', 'swing'};
phaseIdxMap = struct('stance', stanceIdx, 'swing', swingIdx);
phaseMaskMap = struct('stance', stanceMask, 'swing', swingMask);

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

%% ===== STEP 3-6: OLS -> ADAPTIVE LASSO -> POST-SELECTION OLS (PER PHASE) =====
for m = 1:numTargets
    fprintf('\n===== Fitting target muscle: %s (%d/%d) =====\n', muscleNames{m}, m, numTargets);

    for phaseIter = 1:numel(phaseNames)
        phaseName = phaseNames{phaseIter};
        idx = phaseIdxMap.(phaseName);
        A_phase = A(idx, :);
        A_phase_orig = A_orig(idx, :);
        Y_phase = Y(idx, m);

        % --- STEP 3: INITIAL OLS ESTIMATE ---
        c_m = c;
        tolC = optsCommon.tolCensor;
        Uidx = find(Y_phase > c_m + tolC);
        nU = numel(Uidx);

        if nU == 0
            warning('Adaptive_LASSO:NoUncensored', ...
                'No uncensored samples for %s (%s). Using zero OLS coefficients.', ...
                muscleNames{m}, phaseName);
            beta_OLS_std = zeros(p, 1);
            b_OLS_std = c_m;
        else
            A_U = A_phase(Uidx, :);
            y_U = Y_phase(Uidx);
            augA = [A_U, ones(nU, 1)];

            if rank(augA) < min(size(augA))
                warning('Adaptive_LASSO:RankDeficient', ...
                    'OLS for %s (%s) is rank-deficient (nU=%d, p=%d). Using pinv.', ...
                    muscleNames{m}, phaseName, nU, p);
                coeffs = pinv(augA) * y_U;
            else
                coeffs = augA \ y_U;
            end     

            beta_OLS_std = coeffs(1:p);
            b_OLS_std = coeffs(end);
        end

        beta_OLS_orig = beta_OLS_std ./ sA.';
        b_OLS_orig = b_OLS_std - muA * beta_OLS_orig;
        beta_OLS.(phaseName)(:, m) = beta_OLS_orig;
        b_OLS.(phaseName)(m) = b_OLS_orig;

        fprintf('[%s|%s] OLS: nU=%d, nnz(beta_OLS)=%d\n', ...
            muscleNames{m}, phaseName, nU, nnz(abs(beta_OLS_orig) > nzTol));

        % --- STEP 4: COMPUTE ADAPTIVE WEIGHTS ---
        w = 1 ./ max(abs(beta_OLS_std), epsilon_ols).^gamma;
        if any(~isfinite(w))
            error('Adaptive_LASSO:BadWeights', ...
                'Non-finite adaptive weights for %s (%s).', muscleNames{m}, phaseName);
        end
        W.(phaseName)(:, m) = w;

        % --- STEP 5: ADAPTIVE LASSO VIA sgl_fit (LAMBDA PATH) ---
        % Column-scale A to absorb adaptive weights into standard L1.
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

        for kk = 1:nPath
            lambda2_k = lambdaGridThis(kk);
            if lambda2Max <= 0
                betaStarStdK = zeros(p, 1);
                interceptAdaStdK = b0AtZero;
                statsAdaK = struct('exitType', 'degenerate_lambda2max_zero', ...
                    'lambda1', lambda1, 'lambda2', lambda2_k, 'converged', true, 'iters', 0);
            else
                optsM = optsCommon;
                optsM.muscleName = [muscleNames{m}, '_', phaseName, '_path_', num2str(kk)];
                optsM.c = c_m;
                optsM.rho = rho;
                optsM.beta0 = prevBetaStarStd;
                optsM.b0 = prevInterceptStd;
                optsM.doWarmStart = true;

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
                optsRefit.muscleName = [muscleNames{m}, '_', phaseName, '_refit_path_', num2str(kk)];
                optsRefit.c = c_m;
                optsRefit.rho = rho;
                optsRefit.doWarmStart = true;
                optsRefit.beta0 = beta0Refit;
                optsRefit.b0 = interceptAdaStdK;

                [betaRefitStdK, interceptRefitStdK, ~, statsRefitK] = ...
                    sgl_fit(A_refit, Y_phase, 0, 0, optsRefit);
                betaRefitStdK(nonS) = 0;
            end

            yPredAdaK = max(A_phase * betaAdaStdK + interceptAdaStdK, c_m);
            trainErrorAdaptive(kk) = sqrt(mean((Y_phase - yPredAdaK).^2));

            yPredRefitK = max(A_phase * betaRefitStdK + interceptRefitStdK, c_m);
            trainErrorRefit(kk) = sqrt(mean((Y_phase - yPredRefitK).^2));

            betaStarStdPath(:, kk) = betaStarStdK;
            betaAdaStdPath(:, kk) = betaAdaStdK;
            interceptAdaStdPath(kk) = interceptAdaStdK;
            betaRefitStdPath(:, kk) = betaRefitStdK;
            interceptRefitStdPath(kk) = interceptRefitStdK;
            supportPath{kk} = S_k;
            supportSize(kk) = numel(S_k);
            statsAdaPath{kk} = statsAdaK;
            statsRefitPath{kk} = statsRefitK;

            prevBetaStarStd = betaStarStdK;
            prevInterceptStd = interceptAdaStdK;
        end

        validErr = isfinite(trainErrorRefit);
        if ~any(validErr)
            error('Adaptive_LASSO:NoValidPathRMSE', ...
                'No finite refit errors for %s (%s).', muscleNames{m}, phaseName);
        end

        [minErr, localMinPos] = min(trainErrorRefit(validErr));
        validIdx = find(validErr);
        minIdx = validIdx(localMinPos);
        thresholdErr = (1 + lambdaPathConfig.tolSelect) * minErr;

        feasibleIdx = find(validErr & trainErrorRefit <= thresholdErr + 10 * eps(max(1, thresholdErr)));
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

        LambdaPath.(phaseName){m} = pathInfo;

        fprintf(['[%s|%s] lambda2Max=%.6g | selectedIdx=%d/%d | selectedAlpha=%.6g | ' ...
            'selectedLambda2=%.6g | minRefitRMSE=%.6g | selectedRefitRMSE=%.6g | selectedSupportSize=%d\n'], ...
            muscleNames{m}, phaseName, lambda2Max, selectedIdx, nPath, alphaGridThis(selectedIdx), ...
            lambdaGridThis(selectedIdx), minErr, trainErrorRefit(selectedIdx), numel(supportPath{selectedIdx}));

        betaAdaOrig = betaAdaStd ./ sA.';
        interceptAdaOrig = interceptAdaStd - muA * betaAdaOrig;
        betaRefitOrig = betaRefitStd ./ sA.';
        interceptRefitOrig = interceptRefitStd - muA * betaRefitOrig;

        Beta_adaptive.(phaseName)(:, m) = betaAdaOrig;
        Intercept_adaptive.(phaseName)(m) = interceptAdaOrig;
        Beta_refit.(phaseName)(:, m) = betaRefitOrig;
        Intercept_refit.(phaseName)(m) = interceptRefitOrig;

        yPredAda = max(A_phase_orig * betaAdaOrig + interceptAdaOrig, c_m);
        YPred_adaptive.(phaseName){m} = yPredAda;
        RMSE_adaptive.(phaseName)(m) = sqrt(mean((Y_phase - yPredAda).^2));

        yPredRefit = max(A_phase_orig * betaRefitOrig + interceptRefitOrig, c_m);
        YPred_refit.(phaseName){m} = yPredRefit;
        RMSE_refit.(phaseName)(m) = sqrt(mean((Y_phase - yPredRefit).^2));

        YPred_adaptive_full(idx, m) = yPredAda;
        YPred_refit_full(idx, m) = yPredRefit;
    end
end

if any(isnan(YPred_adaptive_full(:))) || any(isnan(YPred_refit_full(:)))
    warning('Adaptive_LASSO:PredictionGap', ...
        'Stitched predictions contain NaN values. Check phase index mapping.');
end

RMSE_adaptive_full = sqrt(mean((Y - YPred_adaptive_full).^2, 1, 'omitnan')).';
RMSE_refit_full = sqrt(mean((Y - YPred_refit_full).^2, 1, 'omitnan')).';

%% ===== STEP 7: SUMMARY METRICS =====
fprintf('\n===== MULTI-MUSCLE ADAPTIVE LASSO SUMMARY (PHASE-SPECIFIC) =====\n');
for m = 1:numTargets
    for phaseIter = 1:numel(phaseNames)
        phaseName = phaseNames{phaseIter};
        pathInfo = LambdaPath.(phaseName){m};
        fprintf(['[%s|%s] OLS_nnz=%d | Ada_nnz=%d | Refit_nnz=%d | ' ...
            'RMSE_ada=%.4g | RMSE_refit=%.4g | exit_ada=%s | exit_refit=%s | ' ...
            'lambda2Max=%.4g | lambda2Sel=%.4g\n'], ...
            muscleNames{m}, phaseName, ...
            nnz(abs(beta_OLS.(phaseName)(:, m)) > nzTol), ...
            nnz(abs(Beta_adaptive.(phaseName)(:, m)) > nzTol), ...
            nnz(abs(Beta_refit.(phaseName)(:, m)) > nzTol), ...
            RMSE_adaptive.(phaseName)(m), RMSE_refit.(phaseName)(m), ...
            stats_ada.(phaseName){m}.exitType, stats_refit.(phaseName){m}.exitType, ...
            pathInfo.lambda2Max, pathInfo.selectedLambda2);
    end
end

%% ===== STEP 8: ANALYZE SPARSE SOLUTION =====
indexToFeatureName = buildInverseFeatureMap(featureToIndexMap, size(X, 2));

for m = 1:numTargets
    for phaseIter = 1:numel(phaseNames)
        phaseName = phaseNames{phaseIter};

        nzIdx = find(abs(Beta_adaptive.(phaseName)(:, m)) > nzTol);
        % fprintf('\nAdaptive LASSO nonzero controls for %s (%s) (nnz = %.0f):\n', ...
            % muscleNames{m}, phaseName, numel(nzIdx));
        if isempty(nzIdx)
            fprintf('(none)\n');
        else
            for k = 1:numel(nzIdx)
                i = nzIdx(k);
                featureName = indexToFeatureName{i};
                % fprintf('  %s -> %s: %.5g  (index %d)\n', featureName, muscleNames{m}, ...
                    % Beta_adaptive.(phaseName)(i, m), i);
            end
        end
        % fprintf('prestimulation (%s, %s adaptive): %.5g\n', ...
            % muscleNames{m}, phaseName, Intercept_adaptive.(phaseName)(m));

        nzIdxRefit = find(abs(Beta_refit.(phaseName)(:, m)) > nzTol);
        fprintf('\nPost-selection OLS nonzero controls for %s (%s) (nnz = %.0f):\n', ...
            muscleNames{m}, phaseName, numel(nzIdxRefit));
        if isempty(nzIdxRefit)
            fprintf('(none)\n');
        else
            for k = 1:numel(nzIdxRefit)
                i = nzIdxRefit(k);
                featureName = indexToFeatureName{i};
                fprintf('  %s -> %s: %.5g  (index %d)\n', featureName, muscleNames{m}, ...
                    Beta_refit.(phaseName)(i, m), i);
            end
        end
        fprintf('prestimulation (%s, %s refit): %.5g\n', ...
            muscleNames{m}, phaseName, Intercept_refit.(phaseName)(m));
    end
end

%% ===== STEP 9: VISUALIZE FIT QUALITY =====
stanceRanges = mask_to_ranges(phaseMaskMap.stance);
swingRanges = mask_to_ranges(phaseMaskMap.swing);
stanceColor = [0.78, 0.86, 0.96]; % blue
swingColor = [0.92, 0.84, 0.84]; % light orange
shadeAlpha = 0.25;

for m = 1:numTargets
    figure('Name', sprintf('Adaptive LASSO Fitting Results: %s', muscleNames{m}), 'Color', 'w');
    vals = [Y(:, m); YPred_adaptive_full(:, m); YPred_refit_full(:, m)];
    vals = vals(isfinite(vals));
    yMin = min(vals);
    yMax = max(vals);
    pad = 0.05 * max(1e-6, yMax - yMin);
    yMin = yMin - pad;
    yMax = yMax + pad;

    hold on;
    for r = 1:size(stanceRanges, 1)
        x1 = stanceRanges(r, 1);
        x2 = stanceRanges(r, 2);
        patch([x1, x2, x2, x1], [yMin, yMin, yMax, yMax], stanceColor, ...
            'EdgeColor', 'none', 'FaceAlpha', shadeAlpha, 'HandleVisibility', 'off');
    end
    for r = 1:size(swingRanges, 1)
        x1 = swingRanges(r, 1);
        x2 = swingRanges(r, 2);
        patch([x1, x2, x2, x1], [yMin, yMin, yMax, yMax], swingColor, ...
            'EdgeColor', 'none', 'FaceAlpha', shadeAlpha, 'HandleVisibility', 'off');
    end

    plot(Y(:, m), 'LineWidth', 1.5);
    plot(YPred_adaptive_full(:, m), 'LineWidth', 1.0);
    plot(YPred_refit_full(:, m), 'LineWidth', 1.0);
    xlim([1 nSamples]);
    ylim([yMin yMax]);
    xlabel('#timesteps');
    ylabel('excitation');
    legend('SCONE muscle excitation', 'Adaptive LASSO', 'Post-selection OLS');
    lambda2Stance = LambdaPath.stance{m}.selectedLambda2;
    lambda2Swing = LambdaPath.swing{m}.selectedLambda2;
    title(sprintf(['%s: Adaptive LASSO (RMSE=%.4g) vs Refit (RMSE=%.4g)\n', ...
        'lambda2 stance=%.4g, swing=%.4g'], ...
        muscleNames{m}, RMSE_adaptive_full(m), RMSE_refit_full(m), ...
        lambda2Stance, lambda2Swing), 'Interpreter','none');
    grid on;
end

%% ===== COLLECT RESULTS =====
results = struct();
results.muscleNames = muscleNames;
results.A = A_orig;
results.A_std = A;
results.featureScaling = struct('muA', muA, 'sA', sA);
results.Y = Y;
results.phase = struct();
results.phase.stanceIdx = stanceIdx;
results.phase.swingIdx = swingIdx;
results.phase.stanceMask = stanceMask;
results.phase.swingMask = swingMask;
results.beta_OLS = beta_OLS;
results.b_OLS = b_OLS;
results.W = W;
results.Beta_adaptive = Beta_adaptive;
results.Intercept_adaptive = Intercept_adaptive;
results.YPred_adaptive = YPred_adaptive_full;
results.YPred_adaptive_phase = YPred_adaptive;
results.RMSE_adaptive = RMSE_adaptive;
results.RMSE_adaptive_full = RMSE_adaptive_full;
results.stats_ada = stats_ada;
results.S = S;
results.Beta_refit = Beta_refit;
results.Intercept_refit = Intercept_refit;
results.YPred_refit = YPred_refit_full;
results.YPred_refit_phase = YPred_refit;
results.RMSE_refit = RMSE_refit;
results.RMSE_refit_full = RMSE_refit_full;
results.stats_refit = stats_refit;
results.lambdaPathConfig = lambdaPathConfig;
results.LambdaPath = LambdaPath;

%% ===== EXPORT TO LASSO CONTROLLER FORMAT =====
% Map the 2-phase LASSO result (stance/swing) to the 5-phase controller
% expected by demo_predictiveForwardSimulation_humanModel.m.
%
% Phase mapping:  stance → phases 0,1    swing → phases 2,3,4
%
% The LASSO feature matrix has 26 rows (9 muscles × 2 + 8 kinematics).
% All 9 muscles and 26 features are exported directly.

% ---- build lasso struct (grouped format) ----
% Two distinct controllers: stance (phases 0,1) and swing (phases 2,3,4).
% The grouped format avoids duplicating parameters for phases that share
% the same controller, making the CMA-ES optimisation more efficient.

betaStance = results.Beta_refit.stance;  % 26 × 9
biasStance = results.Intercept_refit.stance(:)';  % 1 × 9
maskStance = abs(betaStance) > nzTol;

betaSwing = results.Beta_refit.swing;   % 26 × 9
biasSwing = results.Intercept_refit.swing(:)';    % 1 × 9
maskSwing = abs(betaSwing) > nzTol;

lasso = struct();
lasso.format = 'grouped';
lasso.nPhases = 5;
lasso.phaseIds = 0:4;
lasso.nGroups = 2;
lasso.groups = struct();
lasso.groups(1).label  = 'stance';
lasso.groups(1).phases = [0 1];
lasso.groups(1).beta   = betaStance;
lasso.groups(1).bias   = biasStance;
lasso.groups(1).mask   = maskStance;
lasso.groups(2).label  = 'swing';
lasso.groups(2).phases = [2 3 4];
lasso.groups(2).beta   = betaSwing;
lasso.groups(2).bias   = biasSwing;
lasso.groups(2).mask   = maskSwing;

% ---- save ----
outDir = fullfile(fileparts(scriptFullPath), '..', 'results');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
save(fullfile(outDir, 'lasso_controller_result.mat'), 'lasso');

fprintf('\n===== LASSO controller exported (grouped format) =====\n');
fprintf('  File : %s\n', fullfile(outDir, 'lasso_controller_result.mat'));
fprintf('  nPhases = 5, phaseIds = [0 1 2 3 4]\n');
fprintf('  nGroups = 2\n');
fprintf('  Group 1 (stance): phases [0,1],  nnz(beta) = %d\n', nnz(maskStance));
fprintf('  Group 2 (swing):  phases [2,3,4], nnz(beta) = %d\n', nnz(maskSwing));
fprintf('  Feature rows: 26 (9 muscles × 2 + 8 kinematics)\n');
fprintf('  Total optimizable params: %d\n', ...
    nnz(maskStance) + 9 + nnz(maskSwing) + 9);
fprintf('=======================================\n');

%% functions
function indexToFeatureName = buildInverseFeatureMap(featureToIndexMap, numFeatures)
ks = keys(featureToIndexMap);
vs = cell2mat(values(featureToIndexMap));

if numel(ks) ~= numFeatures
    error('buildInverseFeatureMap:CountMismatch', ...
        'feature map size (%d) does not match numFeatures (%d).', numel(ks), numFeatures);
end

if any(vs < 1) || any(vs > numFeatures) || any(abs(vs - round(vs)) > 0)
    error('buildInverseFeatureMap:BadIndex', 'feature indices must be integers in [1, %d].', numFeatures);
end

if numel(unique(vs)) ~= numel(vs)
    error('buildInverseFeatureMap:NonBijective', 'feature map is not one-to-one.');
end

indexToFeatureName = cell(numFeatures, 1);
for t = 1:numel(ks)
    indexToFeatureName{vs(t)} = ks{t};
end
end

function ranges = mask_to_ranges(mask)
mask = mask(:);
if isempty(mask)
    ranges = zeros(0, 2);
    return;
end

d = diff([false; mask; false]);
startIdx = find(d == 1);
endIdx = find(d == -1) - 1;
ranges = [startIdx, endIdx];
end

% TESTING CHECKLIST:
% [ ] With lambda2 very large: all Beta_adaptive.(phase) should be zero
% [ ] With lambda2 = 0: adaptive LASSO should approximate OLS per phase
% [ ] Post-selection OLS nnz should equal adaptive LASSO nnz per phase
% [ ] Post-selection OLS RMSE should be <= adaptive LASSO RMSE per phase
%     (debiasing reduces shrinkage bias)
% [ ] max(abs(A_refit(:,nonS) * Beta_refit(nonS,m))) should be 0
% [ ] YPred >= c for all predictions (censoring floor)
% [ ] stats_ada.(phase){m}.objHist should be non-increasing
