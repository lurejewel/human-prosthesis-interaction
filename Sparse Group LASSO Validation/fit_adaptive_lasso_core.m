function [fit, lasso] = fit_adaptive_lasso_core(D, options)
%FIT_ADAPTIVE_LASSO_CORE  Per-phase adaptive LASSO on a prepared design system.
%
%   [fit, lasso] = fit_adaptive_lasso_core(D, options)
%
%   This is the fitting half of the original fit_adaptive_lasso_phase_model.m.
%   It takes an already loaded/preprocessed design system D (see
%   load_lasso_design_matrices.m) and runs, for every (muscle, phase group):
%
%       y_m(t) = X(t) * beta_{g,m} + b_{g,m} + noise
%
%   with the censored adaptive-LASSO pipeline:
%     1. weighted OLS (optionally Huber IRLS) on uncensored samples;
%     2. adaptive weights w_j = 1/max(|beta_OLS_j|, eps)^gamma;
%     3. column-scaled LASSO along a lambda path (sgl_fit);
%     4. post-selection refit at every lambda, pick by refit RMSE;
%     5. Huber IRLS on the selected support (de-biasing).
%
%   Two modes
%   ---------
%   mode = 'full'          (default) the complete pipeline above.
%   mode = 'mask_override' when options.maskOverride is a {1 x nGroups} cell
%                          of logical masks.  The lambda path is skipped and
%                          the coefficients are re-estimated on the FULL data
%                          with the support frozen to the supplied masks.  This
%                          is how the consensus structure of the multi-seed
%                          identification is turned into a controller.
%
%   Input D (as produced by load_lasso_design_matrices.m)
%   ----------------------------------------------------
%   .X, .X_std, .muX, .sX, .Y, .w, .phaseIdxMap, .muscleNames, .featureNames
%   .stoFile (optional)
%
%   Options: see lasso_fit_default_options.m.  Relevant here:
%     verbose, storeLambdaPath, maskOverride, suppressWarnings,
%     robustInitIter, robustIter, huberDeltaFactor, gamma, epsilonOLS, nzTol,
%     c, rho, lambdaPathConfig, groupLabels, groupPhases, source, preprocessing
%
%   Output fit
%   ----------
%   .Beta_refit / .Intercept_refit      final coefficients (original scale)
%   .YPred_refit / .RMSE_refit          per-muscle-phase prediction / error
%   .RMSE_refit_full                    weighted RMSE over all samples
%   .S                                  selected support per phase/muscle
%   .LambdaPath                         per-lambda diagnostics (mode 'full')
%   .mode                               'full' | 'mask_override'
%   plus the adaptive/OLS snapshots (empty in 'mask_override' mode).
%
%   Output lasso
%   ------------
%   grouped v2 controller: .groups(g).beta/bias/mask/phases/label plus
%   featureNames/muscleNames/source/preprocessing metadata.
%
%   See also load_lasso_design_matrices, fit_adaptive_lasso_phase_model.

if nargin < 2 || isempty(options)
    options = struct();
end
opts = lasso_fit_default_options(options);

%% ===== VALIDATE THE PREPARED SYSTEM ======================================
req = {'X', 'X_std', 'muX', 'sX', 'Y', 'w', 'phaseIdxMap', 'muscleNames'};
for i = 1:numel(req)
    if ~isfield(D, req{i})
        error('fit_adaptive_lasso_core:MissingField', ...
            'D.%s is required (see load_lasso_design_matrices).', req{i});
    end
end

X_orig    = D.X;
X_std     = D.X_std;
muX       = D.muX(:).';
sX        = D.sX(:).';
Y         = D.Y;
wAll      = D.w(:);
muscleNames = D.muscleNames;
nSamples  = size(X_std, 1);
p         = size(X_std, 2);
numTargets = size(Y, 2);

if numel(muscleNames) ~= numTargets
    error('fit_adaptive_lasso_core:MuscleCountMismatch', ...
        'D.muscleNames has %d entries but Y has %d columns.', ...
        numel(muscleNames), numTargets);
end
nMus = (p - 8) / 2;
if nMus < 1 || abs(nMus - round(nMus)) > 1e-12 || nMus ~= numTargets
    error('fit_adaptive_lasso_core:InvalidPredictorCount', ...
        'size(X,2) must satisfy p = 2*#mus + 8 with #mus = %d. Got p = %d.', ...
        numTargets, p);
end
if numel(wAll) ~= nSamples
    error('fit_adaptive_lasso_core:BadWeights', ...
        'D.w has %d entries, expected %d.', numel(wAll), nSamples);
end
if any(~isfinite(wAll)) || any(wAll < 0)
    error('fit_adaptive_lasso_core:BadWeights', ...
        'D.w must be finite and nonnegative.');
end

phaseNames = {'stance', 'swing'};
for i = 1:numel(phaseNames)
    if ~isfield(D.phaseIdxMap, phaseNames{i})
        error('fit_adaptive_lasso_core:MissingPhase', ...
            'D.phaseIdxMap.%s is required.', phaseNames{i});
    end
end

nGroups = numel(opts.groupPhases);
maskOverride = opts.maskOverride;
useOverride = ~isempty(maskOverride);
if useOverride
    if numel(maskOverride) ~= nGroups
        error('fit_adaptive_lasso_core:BadMaskOverride', ...
            'options.maskOverride must have %d entries, got %d.', ...
            nGroups, numel(maskOverride));
    end
    for g = 1:nGroups
        if ~isequal(size(maskOverride{g}), [p, numTargets])
            error('fit_adaptive_lasso_core:BadMaskOverride', ...
                'options.maskOverride{%d} must be %d x %d, got %d x %d.', ...
                g, p, numTargets, size(maskOverride{g}, 1), size(maskOverride{g}, 2));
        end
        maskOverride{g} = logical(maskOverride{g});
    end
end
if opts.suppressWarnings
    prevWarn = warning('off', 'fit_adaptive_lasso_phase_model:RankDeficient');
    warnCleanup = onCleanup(@() warning(prevWarn));
end

gamma = opts.gamma;
epsilon_ols = opts.epsilonOLS;
nzTol = opts.nzTol;
lambda1 = opts.lambda1;
lambdaPathConfig = opts.lambdaPathConfig;
c = opts.c;
rho = opts.rho;

assert(lambda1 == 0, ...
    'Adaptive LASSO column-scaling only adapts L1 penalty. lambda1 must be 0.');

%% ===== STEP 2: OPTIMIZATION SETTINGS ====================================
optsCommon = struct();
optsCommon.maxIter = 20000;
optsCommon.tol = 1e-10;
optsCommon.verbose = false;
optsCommon.doWarmStart = false;
optsCommon.backtrackBeta = 0.5;
optsCommon.L0 = 1000;
optsCommon.tolCensor = 1e-12;

%% ===== STEP 2B: PHASE INDICES (ALREADY DELAY-ALIGNED) ===================
% load_lasso_design_matrices already shifted the phase indices by delaySteps
% and dropped out-of-range samples.  Here we only rebuild the boolean masks
% and check that the two groups stay disjoint.
stanceIdx = D.phaseIdxMap.stance(:);
swingIdx  = D.phaseIdxMap.swing(:);
if any(stanceIdx < 1 | stanceIdx > nSamples) || ...
        any(swingIdx < 1 | swingIdx > nSamples)
    error('fit_adaptive_lasso_core:BadPhaseIndex', ...
        'Phase indices must lie in [1, %d].', nSamples);
end

stanceMask = false(nSamples, 1);
stanceMask(stanceIdx) = true;
swingMask = false(nSamples, 1);
swingMask(swingIdx) = true;

if any(stanceMask & swingMask)
    error('fit_adaptive_lasso_core:PhaseOverlap', ...
        'Stance and swing indices overlap.');
end

phaseIdxMap = struct('stance', stanceIdx, 'swing', swingIdx);

%% ===== PRE-ALLOCATIONS ==================================================
beta_OLS = struct();
b_OLS = struct();
W = struct();
Beta_adaptive = struct();
Intercept_adaptive = struct();
YPred_adaptive = struct();
RMSE_adaptive = struct();
stats_ada = struct();
Support = struct();          % selected support (per phase, per muscle)
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
    Support.(phaseName) = cell(numTargets, 1);
    LambdaPath.(phaseName) = cell(numTargets, 1);

    Beta_refit.(phaseName) = zeros(p, numTargets);
    Intercept_refit.(phaseName) = zeros(numTargets, 1);
    YPred_refit.(phaseName) = cell(numTargets, 1);
    RMSE_refit.(phaseName) = zeros(numTargets, 1);
    stats_refit.(phaseName) = cell(numTargets, 1);
end

YPred_adaptive_full = nan(nSamples, numTargets);
YPred_refit_full = nan(nSamples, numTargets);

%% ===== STEP 3-6: OLS -> ADAPTIVE LASSO -> ROBUST REFIT ==================
% Outer loop: one fit per muscle (shared design matrix, different target).
% Inner loop: one fit per phase group (stance and swing obey different laws).
% Muscles cannot be pooled because the censored/uncensored split depends on y.
for m = 1:numTargets
    if opts.verbose
        fprintf('\n===== Fitting target muscle: %s (%d/%d) =====\n', ...
            muscleNames{m}, m, numTargets);
    end

    for phaseIter = 1:numel(phaseNames)
        phaseName = phaseNames{phaseIter};
        idx = phaseIdxMap.(phaseName);
        X_phase = X_std(idx, :);
        X_phase_orig = X_orig(idx, :);
        Y_phase = Y(idx, m);
        wPhase = wAll(idx);
        if sum(wPhase) <= 0
            error('fit_adaptive_lasso_core:NoValidSamples', ...
                'All samples are excluded for %s (%s).', ...
                muscleNames{m}, phaseName);
        end

        c_m = c;
        tolC = optsCommon.tolCensor;

        if useOverride
            %% ---- mask-constrained refit on the frozen support -----------
            S_k = find(maskOverride{phaseIter}(:, m));
            [betaRefitStdK, interceptRefitStdK, statsRefitFinal] = ...
                ols_support_init(X_phase, Y_phase, wPhase, S_k, c_m, tolC, ...
                muscleNames{m}, phaseName);
            % The "adaptive" fields are not produced in this mode; mirror the
            % refit values so no NaN propagates into the (then discarded)
            % adaptive snapshots.
            betaAdaStd = betaRefitStdK;
            interceptAdaStd = interceptRefitStdK;
            betaRefitStd = betaRefitStdK;
            interceptRefitStd = interceptRefitStdK;
            Support.(phaseName){m} = S_k(:);
            stats_ada.(phaseName){m} = struct( ...
                'exitType', 'skipped_mask_override', 'lambda1', 0, ...
                'lambda2', 0, 'converged', true, 'iters', 0);
            pathInfo = struct();
            pathInfo.mode = 'mask_override';
            pathInfo.selectedIdx = 0;
            pathInfo.selectedAlpha = 1;
            pathInfo.selectedLambda2 = 0;
            pathInfo.selectedSupport = S_k(:);
            pathInfo.supportSize = numel(S_k);
            pathInfo.minTrainErrorRefit = NaN;
            pathInfo.thresholdTrainErrorRefit = NaN;
            pathInfo.selectedTrainErrorRefit = NaN;
            pathInfo.selectedTrainErrorAdaptive = NaN;
            pathInfo.tolSelect = lambdaPathConfig.tolSelect;
            pathInfo.selectionMetric = lambdaPathConfig.selectionMetric;
        else
            %% ---- STEP 3: INITIAL WEIGHTED OLS ESTIMATE ------------------
            % Only uncensored samples (y > c+tol) are true observations; the
            % adaptive weights are built from this OLS estimate.
            Uidx = find(Y_phase > c_m + tolC);
            nU = numel(Uidx);

            if nU == 0
                warning('fit_adaptive_lasso_phase_model:NoUncensored', ...
                    'No uncensored samples for %s (%s). Using zero OLS.', ...
                    muscleNames{m}, phaseName);
                beta_OLS_std = zeros(p, 1);
                b_OLS_std = c_m;
            else
                X_U = X_phase(Uidx, :);
                y_U = Y_phase(Uidx);
                w_U = max(wPhase(Uidx), 0);
                augX = [X_U, ones(nU, 1)];
                rankDef = rank(augX) < min(size(augX));
                if rankDef
                    warning('fit_adaptive_lasso_phase_model:RankDeficient', ...
                        'OLS for %s (%s) is rank-deficient (nU=%d, p=%d). Using pinv.', ...
                        muscleNames{m}, phaseName, nU, p);
                end

                coeffs = solve_weighted_ls(augX, y_U, w_U, rankDef);

                if opts.robustInitIter > 0 && nU >= 5 && sum(w_U) > 0
                    for it = 1:opts.robustInitIter
                        res = y_U - augX * coeffs;
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
                        coeffs = solve_weighted_ls(augX, y_U, w_U .* wH, rankDef);
                    end
                end

                beta_OLS_std = coeffs(1:p);
                b_OLS_std = coeffs(end);
            end

            beta_OLS_orig = beta_OLS_std ./ sX.';
            b_OLS_orig = b_OLS_std - muX * beta_OLS_orig;
            beta_OLS.(phaseName)(:, m) = beta_OLS_orig;
            b_OLS.(phaseName)(m) = b_OLS_orig;

            %% ---- STEP 4: ADAPTIVE WEIGHTS ------------------------------
            % Zou (2006): large OLS coefficients get a small penalty.  The
            % epsilon_ols floor avoids division by zero.
            w = 1 ./ max(abs(beta_OLS_std), epsilon_ols).^gamma;
            if any(~isfinite(w))
                error('fit_adaptive_lasso_core:BadWeights', ...
                    'Non-finite adaptive weights for %s (%s).', ...
                    muscleNames{m}, phaseName);
            end
            W.(phaseName)(:, m) = w;

            %% ---- STEP 5: ADAPTIVE LASSO VIA sgl_fit (LAMBDA PATH) ------
            % Column scaling X_star = X ./ w turns the weighted L1 penalty
            % into a plain one.  The path starts at the censored lambda2Max
            % and is walked from strong to weak penalty with warm starts.
            X_star = X_phase .* (1 ./ w)';
            [lambda2Max, b0AtZero, lambdaMaxInfo] = ...
                compute_lambda2_max_censored(X_star, Y_phase, c_m, rho, tolC);

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
                        sgl_fit(X_star, Y_phase, lambda1, lambda2_k, optsM);
                end

                betaAdaStdK = betaStarStdK ./ w;
                S_k = find(abs(betaAdaStdK) > nzTol);

                if isempty(S_k)
                    betaRefitStdK = zeros(p, 1);
                    interceptRefitStdK = b0AtZero;
                    statsRefitK = struct('exitType', 'intercept_only_empty_support', ...
                        'lambda1', 0, 'lambda2', 0, 'converged', true, 'iters', 0);
                else
                    X_refit = X_phase;
                    nonS = setdiff(1:p, S_k);
                    X_refit(:, nonS) = 0;
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
                        sgl_fit(X_refit, Y_phase, 0, 0, optsRefit);
                    betaRefitStdK(nonS) = 0;
                end

                yPredAdaK = max(X_phase * betaAdaStdK + interceptAdaStdK, c_m);
                trainErrorAdaptive(kk) = weighted_rmse(Y_phase, yPredAdaK, wPhase);

                yPredRefitK = max(X_phase * betaRefitStdK + interceptRefitStdK, c_m);
                trainErrorRefit(kk) = weighted_rmse(Y_phase, yPredRefitK, wPhase);

                betaStarStdPath(:, kk) = betaStarStdK;
                betaAdaStdPath(:, kk) = betaAdaStdK;
                interceptAdaStdPath(kk) = interceptAdaStdK;
                betaRefitStdPath(:, kk) = betaRefitStdK;
                interceptRefitStdPath(kk) = interceptRefitStdK;
                supportPath{kk} = S_k;
                supportSize(kk) = numel(S_k);
                statsAdaPath{kk} = compact_stats(statsAdaK);
                statsRefitPath{kk} = compact_stats(statsRefitK);

                prevBetaStarStd = betaStarStdK;
                prevInterceptStd = interceptAdaStdK;
            end

            % Selection rule (unchanged): within (1+tolSelect) of the minimum
            % refit RMSE, take the sparsest lambda (first feasible on the path,
            % which is walked from strong to weak penalty).
            validErr = isfinite(trainErrorRefit);
            if ~any(validErr)
                error('fit_adaptive_lasso_core:NoValidPathRMSE', ...
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

            interceptAdaStd = interceptAdaStdPath(selectedIdx);
            betaAdaStd = betaAdaStdPath(:, selectedIdx);
            betaRefitStd = betaRefitStdPath(:, selectedIdx);
            interceptRefitStd = interceptRefitStdPath(selectedIdx);
            S_k = supportPath{selectedIdx};
            Support.(phaseName){m} = supportPath{selectedIdx};
            stats_ada.(phaseName){m} = statsAdaPath{selectedIdx};
            statsRefitFinal = statsRefitPath{selectedIdx};

            pathInfo = struct();
            pathInfo.mode = 'full';
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
            if opts.storeLambdaPath
                pathInfo.betaStarStdPath = betaStarStdPath;
                pathInfo.betaAdaStdPath = betaAdaStdPath;
                pathInfo.interceptAdaStdPath = interceptAdaStdPath;
                pathInfo.betaRefitStdPath = betaRefitStdPath;
                pathInfo.interceptRefitStdPath = interceptRefitStdPath;
            end
            pathInfo.statsAdaPath = statsAdaPath;
            pathInfo.statsRefitPath = statsRefitPath;
            pathInfo.robustIter = opts.robustIter;
        end

        %% ---- STEP 6: ROBUST HUBER IRLS ON THE SELECTED SUPPORT ---------
        % The path coefficients are only a warm start; the final values come
        % from Huber IRLS on the fixed support.  In mask_override mode the
        % support is the frozen consensus mask and the same loop is used.
        [betaRefitStdK, interceptRefitStdK, statsRefitFinal, wCur] = ...
            robust_support_refit(X_phase, Y_phase, wPhase, S_k, c_m, rho, ...
            optsCommon, opts, betaRefitStd, interceptRefitStd, ...
            statsRefitFinal, muscleNames{m}, phaseName);
        stats_refit.(phaseName){m} = compact_stats(statsRefitFinal);
        pathInfo.sampleWeight = wPhase;
        pathInfo.robustFinalWeights = wCur;
        LambdaPath.(phaseName){m} = pathInfo;

        if opts.verbose
            selLambda2 = 0;
            if isfield(pathInfo, 'selectedLambda2')
                selLambda2 = pathInfo.selectedLambda2;
            end
            fprintf(['[%s|%s] selectedIdx=%d | selectedLambda2=%.6g | ' ...
                'RMSE=%.6g | nnz=%d\n'], ...
                muscleNames{m}, phaseName, pathInfo.selectedIdx, selLambda2, ...
                weighted_rmse(Y_phase, max(X_phase * betaRefitStdK + ...
                interceptRefitStdK, c_m), wPhase), numel(S_k));
        end

        betaAdaOrig = betaAdaStd ./ sX.';
        interceptAdaOrig = interceptAdaStd - muX * betaAdaOrig;
        betaRefitOrig = betaRefitStdK ./ sX.';
        interceptRefitOrig = interceptRefitStdK - muX * betaRefitOrig;

        Beta_adaptive.(phaseName)(:, m) = betaAdaOrig;
        Intercept_adaptive.(phaseName)(m) = interceptAdaOrig;
        Beta_refit.(phaseName)(:, m) = betaRefitOrig;
        Intercept_refit.(phaseName)(m) = interceptRefitOrig;

        yPredAda = max(X_phase_orig * betaAdaOrig + interceptAdaOrig, c_m);
        YPred_adaptive.(phaseName){m} = yPredAda;
        RMSE_adaptive.(phaseName)(m) = weighted_rmse(Y_phase, yPredAda, wPhase);

        yPredRefit = max(X_phase_orig * betaRefitOrig + interceptRefitOrig, c_m);
        YPred_refit.(phaseName){m} = yPredRefit;
        RMSE_refit.(phaseName)(m) = weighted_rmse(Y_phase, yPredRefit, wPhase);

        YPred_adaptive_full(idx, m) = yPredAda;
        YPred_refit_full(idx, m) = yPredRefit;
    end
end

if ~useOverride && ...
        (any(isnan(YPred_adaptive_full(:))) || any(isnan(YPred_refit_full(:))))
    warning('fit_adaptive_lasso_phase_model:PredictionGap', ...
        'Stitched predictions contain NaN values. Check phase mapping.');
end

if useOverride
    RMSE_adaptive_full = nan(1, numTargets);
else
    RMSE_adaptive_full = sqrt(sum(wAll .* (Y - YPred_adaptive_full).^2, 1) ./ ...
        max(sum(wAll), eps)).';
end
RMSE_refit_full = sqrt(sum(wAll .* (Y - YPred_refit_full).^2, 1) ./ ...
    max(sum(wAll), eps)).';

%% ===== SUMMARY ==========================================================
if opts.verbose && ~useOverride
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
elseif opts.verbose
    fprintf('\n===== MASK-CONSTRAINED REFIT SUMMARY =====\n');
    for m = 1:numTargets
        fprintf('[%s] frozen nnz: stance=%d swing=%d | RMSE_refit stance=%.4g swing=%.4g\n', ...
            muscleNames{m}, numel(Support.stance{m}), numel(Support.swing{m}), ...
            RMSE_refit.stance(m), RMSE_refit.swing(m));
    end
end

%% ===== ASSEMBLE OUTPUTS =================================================
stoFile = '';
if isfield(D, 'stoFile')
    stoFile = D.stoFile;
end
if isfield(D, 'featureNames')
    featureNames = D.featureNames;
else
    featureNames = {};
end

fit = struct();
fit.source       = opts.source;
fit.stoFile      = stoFile;
fit.mode         = ternary_str(useOverride, 'mask_override', 'full');
fit.refitOnly    = useOverride;
fit.muscleNames  = muscleNames;
fit.X            = X_orig;
fit.X_std        = X_std;
fit.featureScaling = struct('muX', muX, 'sX', sX);
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
fit.S            = Support;
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
fit.featureNames = featureNames;

% In mask_override mode the adaptive snapshots were never produced; remove
% them instead of leaving misleading numbers behind.
if useOverride
    fit.beta_OLS = [];
    fit.b_OLS = [];
    fit.W = [];
    fit.Beta_adaptive = [];
    fit.Intercept_adaptive = [];
    fit.YPred_adaptive = [];
    fit.YPred_adaptive_phase = [];
    fit.RMSE_adaptive = [];
    fit.RMSE_adaptive_full = [];
    fit.stats_ada = [];
    fit.LambdaPath = [];
end

% ---- grouped controller struct ----
% beta entries at/below nzTol are forced to exact structural zeros because
% load_lasso_reflex_controller.m rejects any non-zero where mask is false.
% In mask_override mode the mask comes from the consensus and the coefficients
% are left untouched (a small in-mask coefficient is still a real parameter).
if useOverride
    betaStance = fit.Beta_refit.stance;
    biasStance = fit.Intercept_refit.stance(:)';
    maskStance = maskOverride{1};

    betaSwing  = fit.Beta_refit.swing;
    biasSwing  = fit.Intercept_refit.swing(:)';
    maskSwing  = maskOverride{2};
else
    betaStance = fit.Beta_refit.stance;
    betaStance(abs(betaStance) <= nzTol) = 0;
    biasStance = fit.Intercept_refit.stance(:)';
    maskStance = abs(betaStance) > nzTol;

    betaSwing  = fit.Beta_refit.swing;
    betaSwing(abs(betaSwing) <= nzTol) = 0;
    biasSwing  = fit.Intercept_refit.swing(:)';
    maskSwing  = abs(betaSwing) > nzTol;
end

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

lasso.featureNames = featureNames;
lasso.muscleNames  = muscleNames;
lasso.source       = opts.source;
lasso.datasetFile  = stoFile;
lasso.created      = char(datetime);
lasso.preprocessing = opts.preprocessing;

if opts.verbose
    fprintf('\n===== LASSO controller assembled (grouped format) =====\n');
    fprintf('  stance nnz(beta) = %d, swing nnz(beta) = %d\n', ...
        nnz(maskStance), nnz(maskSwing));
    fprintf('  total optimizable parameters = %d\n', ...
        nnz(maskStance) + 9 + nnz(maskSwing) + 9);
    fprintf('=======================================================\n');
end
end

%% ========================================================================
% Local: robust Huber IRLS on a fixed support.
% ========================================================================
% Shared by the normal pipeline (warm start = the selected lambda solution)
% and by the mask-constrained refit (warm start = a support-restricted OLS).
% Only uncensored samples get residual-scale Huber weights; censored samples
% keep their base weight so rho still enforces the hinge constraint.
function [betaS, bS, statsFinal, wCur] = robust_support_refit( ...
    X_phase, Y_phase, wPhase, S_k, c_m, rho, optsCommon, opts, ...
    beta0, b0, stats0, muscleName, phaseName)
p = size(X_phase, 2);
nonS = setdiff(1:p, S_k);
X_refit = X_phase;
X_refit(:, nonS) = 0;
wCur = wPhase;
betaS = beta0(:);
bS = b0;
statsFinal = stats0;

for it = 1:opts.robustIter
    optsRefit = optsCommon;
    optsRefit.muscleName = [muscleName, '_', phaseName, '_robust_', num2str(it)];
    optsRefit.c = c_m;
    optsRefit.rho = rho;
    optsRefit.doWarmStart = true;
    optsRefit.beta0 = betaS;
    optsRefit.b0 = bS;
    optsRefit.sampleWeight = max(wCur, 0);

    [betaS, bS, ~, statsFinal] = sgl_fit(X_refit, Y_phase, 0, 0, optsRefit);
    betaS(nonS) = 0;

    latent = X_phase * betaS + bS;
    isCensored = Y_phase <= c_m + optsCommon.tolCensor;
    fitU = ~isCensored & wPhase > 0;

    if any(fitU)
        resU = Y_phase(fitU) - latent(fitU);
        s = 1.4826 * median(abs(resU - median(resU)));
        if ~isfinite(s) || s < 1e-8
            s = 1e-8;
        end
        delta = opts.huberDeltaFactor * s;

        wH = ones(size(Y_phase));
        wH(fitU) = min(1, delta ./ max(abs(resU), eps));
        wCur = wPhase .* wH;
    else
        wCur = wPhase;
    end
end
end

%% ========================================================================
% Local: support-restricted weighted OLS, used to initialise the frozen-mask
% refit (mirrors STEP 3 but keeps only the frozen columns).
% ========================================================================
function [betaS, bS, stats0] = ols_support_init( ...
    X_phase, Y_phase, wPhase, S_k, c_m, tolC, muscleName, phaseName)
p = size(X_phase, 2);
nonS = setdiff(1:p, S_k);
X_refit = X_phase;
X_refit(:, nonS) = 0;

Uidx = find(Y_phase > c_m + tolC);
nU = numel(Uidx);
if nU == 0
    betaS = zeros(p, 1);
    bS = c_m;
else
    augX = [X_refit(Uidx, :), ones(nU, 1)];
    rankDef = rank(augX) < min(size(augX));
    coeffs = solve_weighted_ls(augX, Y_phase(Uidx), ...
        max(wPhase(Uidx), 0), rankDef);
    betaS = coeffs(1:p);
    betaS(nonS) = 0;
    bS = coeffs(end);
end
stats0 = struct('exitType', 'ols_support_init', 'lambda1', 0, 'lambda2', 0, ...
    'converged', true, 'iters', 0, 'muscle', muscleName, 'phase', phaseName);
end

%% ========================================================================
% Local: weighted least squares via row scaling.
% ========================================================================
function coeffs = solve_weighted_ls(augX, y, w, usePinv)
w = max(w(:), 0);
sw = sqrt(w);
if sum(w) <= 0
    sw = ones(size(w));
end
if usePinv
    coeffs = pinv(augX .* sw) * (y .* sw);
else
    coeffs = (augX .* sw) \ (y .* sw);
end
end

%% ========================================================================
% Local: weighted RMSE used for lambda selection and reporting.
% ========================================================================
function r = weighted_rmse(y, yhat, w)
w = max(w(:), 0);
if sum(w) <= 0
    r = sqrt(mean((y(:) - yhat(:)).^2));
else
    r = sqrt(sum(w .* (y(:) - yhat(:)).^2) / sum(w));
end
end

%% ========================================================================
% Local: keep only the small, useful part of an sgl_fit stats struct.
% ========================================================================
function s = compact_stats(st)
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
% Local: tiny ternary for strings (avoids an if/else in a struct literal).
% ========================================================================
function out = ternary_str(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
