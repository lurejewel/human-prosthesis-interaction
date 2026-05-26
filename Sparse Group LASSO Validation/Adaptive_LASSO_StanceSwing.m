% Adaptive LASSO runner script.
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
% - gamma & lambda的二维网格搜索
% - 消融实验：pure lasso vs adaptive lasso vs adaptive lasso with priori vs
% something else
% - c, rho，以及其他参数的灵敏度分析
% ========================================================================

close all; clear; clc;

%% ===== STEP 1: DATA LOADING & PREPROCESSING =====
% Load target muscle excitations and one shared biomechanical feature matrix.
stoFilePath = 'UN.sto';
muscleNames = {'hamstrings_r', 'bifemsh_r', 'glut_max_r', 'iliopsoas_r', ...
    'rect_fem_r', 'vasti_r'};
numTargets = numel(muscleNames);

% Tunable parameters (adaptive LASSO + censoring)
gamma = 1;
epsilon_ols = 1e-8;
nzTol = 1e-5;
lambda1 = 0;
lambda2 = 2;
c = 0.01;
rho = 50;
delaySteps = 1;
excitationThreshold = 0.01;

assert(lambda1 == 0, ...
    'Adaptive LASSO column-scaling only adapts L1 penalty. lambda1 must be 0.');

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
optsCommon.doCenter = false; % true;
optsCommon.doScale = false; % true;
optsCommon.scaleType = 'std';
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
            beta_OLS.(phaseName)(:, m) = zeros(p, 1);
            b_OLS.(phaseName)(m) = c_m;
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

            beta_OLS.(phaseName)(:, m) = coeffs(1:p);
            b_OLS.(phaseName)(m) = coeffs(end);
        end

        fprintf('[%s|%s] OLS: nU=%d, nnz(beta_OLS)=%d\n', ...
            muscleNames{m}, phaseName, nU, nnz(abs(beta_OLS.(phaseName)(:, m)) > nzTol));

        % --- STEP 4: COMPUTE ADAPTIVE WEIGHTS ---
        w = 1 ./ max(abs(beta_OLS.(phaseName)(:, m)), epsilon_ols).^gamma;
        if any(~isfinite(w))
            error('Adaptive_LASSO:BadWeights', ...
                'Non-finite adaptive weights for %s (%s).', muscleNames{m}, phaseName);
        end
        W.(phaseName)(:, m) = w;

        % --- STEP 5: ADAPTIVE LASSO VIA sgl_fit ---
        % Column-scale A to absorb adaptive weights into standard L1.
        A_star = A_phase .* (1 ./ w)';

        optsM = optsCommon;
        optsM.muscleName = [muscleNames{m}, '_', phaseName];
        optsM.c = c;
        optsM.rho = rho;
        optsM.beta0 = zeros(p, 1);
        optsM.b0 = 0;
        optsM.doWarmStart = true;

        % Optional warm start from OLS (in A* space):
        % optsM.beta0 = beta_OLS.(phaseName)(:, m) .* w;
        % optsM.b0 = b_OLS.(phaseName)(m);

        [betaOrig_star, intercept_star, ~, stats_ada.(phaseName){m}] = ...
            sgl_fit(A_star, Y_phase, lambda1, lambda2, optsM);

        Beta_adaptive.(phaseName)(:, m) = betaOrig_star ./ w;
        Intercept_adaptive.(phaseName)(m) = intercept_star;

        yPredAda = max(A_phase * Beta_adaptive.(phaseName)(:, m) + ...
            Intercept_adaptive.(phaseName)(m), c_m);
        YPred_adaptive.(phaseName){m} = yPredAda;
        RMSE_adaptive.(phaseName)(m) = sqrt(mean((Y_phase - yPredAda).^2));

        S.(phaseName){m} = find(abs(Beta_adaptive.(phaseName)(:, m)) > nzTol);

        % --- STEP 6: POST-SELECTION OLS VIA sgl_fit ---
        if isempty(S.(phaseName){m})
            warning('Adaptive_LASSO:EmptySupport', ...
                'Support set for %s (%s) is empty. Skipping post-selection OLS.', ...
                muscleNames{m}, phaseName);
            Beta_refit.(phaseName)(:, m) = zeros(p, 1);
            Intercept_refit.(phaseName)(m) = 0;
            yPredRefit = c_m * ones(numel(idx), 1);
            YPred_refit.(phaseName){m} = yPredRefit;
            RMSE_refit.(phaseName)(m) = sqrt(mean((Y_phase - yPredRefit).^2));
            stats_refit.(phaseName){m} = struct('exitType', 'skipped_empty_support');
        else
            A_refit = A_phase;
            nonS = setdiff(1:p, S.(phaseName){m});
            A_refit(:, nonS) = 0;

            optsRefit = optsCommon;
            optsRefit.muscleName = [muscleNames{m}, '_', phaseName, '_refit'];
            optsRefit.c = c;
            optsRefit.rho = rho;
            optsRefit.doWarmStart = true;
            optsRefit.beta0 = Beta_adaptive.(phaseName)(:, m);
            optsRefit.b0 = Intercept_adaptive.(phaseName)(m);

            [Beta_refit.(phaseName)(:, m), Intercept_refit.(phaseName)(m), ~, ...
                stats_refit.(phaseName){m}] = sgl_fit(A_refit, Y_phase, 0, 0, optsRefit);

            yPredRefit = max(A_phase * Beta_refit.(phaseName)(:, m) + ...
                Intercept_refit.(phaseName)(m), c_m);
            YPred_refit.(phaseName){m} = yPredRefit;
            RMSE_refit.(phaseName)(m) = sqrt(mean((Y_phase - yPredRefit).^2));
        end

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
        fprintf(['[%s|%s] OLS_nnz=%d | Ada_nnz=%d | Refit_nnz=%d | ' ...
            'RMSE_ada=%.4g | RMSE_refit=%.4g | exit_ada=%s | exit_refit=%s\n'], ...
            muscleNames{m}, phaseName, ...
            nnz(abs(beta_OLS.(phaseName)(:, m)) > nzTol), ...
            nnz(abs(Beta_adaptive.(phaseName)(:, m)) > nzTol), ...
            nnz(abs(Beta_refit.(phaseName)(:, m)) > nzTol), ...
            RMSE_adaptive.(phaseName)(m), RMSE_refit.(phaseName)(m), ...
            stats_ada.(phaseName){m}.exitType, stats_refit.(phaseName){m}.exitType);
    end
end

%% ===== STEP 8: ANALYZE SPARSE SOLUTION =====
indexToFeatureName = buildInverseFeatureMap(featureToIndexMap, size(X, 2));

for m = 1:numTargets
    for phaseIter = 1:numel(phaseNames)
        phaseName = phaseNames{phaseIter};

        nzIdx = find(abs(Beta_adaptive.(phaseName)(:, m)) > nzTol);
        fprintf('\nAdaptive LASSO nonzero controls for %s (%s) (nnz = %.0f):\n', ...
            muscleNames{m}, phaseName, numel(nzIdx));
        if isempty(nzIdx)
            fprintf('(none)\n');
        else
            for k = 1:numel(nzIdx)
                i = nzIdx(k);
                featureName = indexToFeatureName{i};
                fprintf('  %s -> %s: %.5g  (index %d)\n', featureName, muscleNames{m}, ...
                    Beta_adaptive.(phaseName)(i, m), i);
            end
        end
        fprintf('prestimulation (%s, %s adaptive): %.5g\n', ...
            muscleNames{m}, phaseName, Intercept_adaptive.(phaseName)(m));

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
    title(sprintf('%s: Adaptive LASSO (RMSE=%.4g) vs Refit (RMSE=%.4g)', ...
        muscleNames{m}, RMSE_adaptive_full(m), RMSE_refit_full(m)), 'Interpreter','none');
    grid on;
end

%% ===== COLLECT RESULTS =====
results = struct();
results.muscleNames = muscleNames;
results.A = A;
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
