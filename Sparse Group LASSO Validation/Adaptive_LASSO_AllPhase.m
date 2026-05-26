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
% - 按照支撑相和摆动相切分数据，对其他几个肌肉的控制律进行验证
% - 所有肌肉成功后的行走正动力学验证
% - 基于多人行走实测数据的肌肉控制律预测和分析
% - 其他具备现有控制律的动作（跑步？）的验证
% - 其他动作（跑步、上下楼梯、蹲起、游泳等周期性动作，以及跳跃等非周期性动作，现有预测动力学无法做到的）的仿真及实验数据下的预测动力学验证
% - 病态动作的验证
% 
% IMPORTANT:
% - maskMatrix需要根据肌肉特性赋值（重读经典论文，导入先验证据）；是否还需要sparse group？
% - (正动力学完成后）代理模型研究：对应Annuals Review of Biomedical Engineering中的大脑预测模型）
% - gamma & lambda的二维网格搜索
% - 消融实验：pure lasso vs adaptive lasso vs adaptive lasso with priori vs
% something else
% - c, rho，以及其他参数的灵敏度分析
% ========================================================================

close all; clear; clc;

%% ===== STEP 1: DATA LOADING & PREPROCESSING =====
% Load three muscle excitations and one shared biomechanical feature matrix.
stoFilePath = 'UN.sto';
muscleNames = {'soleus_r', 'tib_ant_r', 'gastroc_r'};
numTargets = numel(muscleNames);

% Tunable parameters (adaptive LASSO + censoring)
gamma = 1;
epsilon_ols = 1e-8;
nzTol = 1e-5;
lambda1Vec = [0, 0, 0];
lambda2Vec = [2, 2, 2];
cVec = [0.01, 0.01, 0.01];
rhoVec = [50, 50, 50];
delaySteps = 1;
excitationThreshold = 0.01;

assert(all(lambda1Vec == 0), ...
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

%% ===== PRE-ALLOCATIONS =====
nSamples = size(A, 1);

beta_OLS = zeros(p, numTargets);
b_OLS = zeros(numTargets, 1);
W = zeros(p, numTargets);

Beta_adaptive = zeros(p, numTargets);
Intercept_adaptive = zeros(numTargets, 1);
YPred_adaptive = zeros(nSamples, numTargets);
RMSE_adaptive = zeros(numTargets, 1);
stats_ada = cell(numTargets, 1);
S = cell(numTargets, 1);

Beta_refit = zeros(p, numTargets);
Intercept_refit = zeros(numTargets, 1);
YPred_refit = zeros(nSamples, numTargets);
RMSE_refit = zeros(numTargets, 1);
stats_refit = cell(numTargets, 1);

%% ===== STEP 3-6: OLS -> ADAPTIVE LASSO -> POST-SELECTION OLS =====
for m = 1:numTargets
    fprintf('\n===== Fitting target muscle: %s (%d/%d) =====\n', muscleNames{m}, m, numTargets);

    % --- STEP 3: INITIAL OLS ESTIMATE ---
    c_m = cVec(m);
    tolC = optsCommon.tolCensor;
    Uidx = find(Y(:, m) > c_m + tolC);
    nU = numel(Uidx);

    if nU == 0
        warning('Adaptive_LASSO:NoUncensored', ...
            'No uncensored samples for %s. Using zero OLS coefficients.', muscleNames{m});
        beta_OLS(:, m) = zeros(p, 1);
        b_OLS(m) = c_m;
    else
        A_U = A(Uidx, :);
        y_U = Y(Uidx, m);
        augA = [A_U, ones(nU, 1)];

        if rank(augA) < min(size(augA))
            warning('Adaptive_LASSO:RankDeficient', ...
                'OLS for %s is rank-deficient (nU=%d, p=%d). Using pinv.', ...
                muscleNames{m}, nU, p);
            coeffs = pinv(augA) * y_U;
        else
            coeffs = augA \ y_U;
        end

        beta_OLS(:, m) = coeffs(1:p);
        b_OLS(m) = coeffs(end);
    end

    fprintf('[%s] OLS: nU=%d, nnz(beta_OLS)=%d\n', ...
        muscleNames{m}, nU, nnz(abs(beta_OLS(:, m)) > nzTol));

    % --- STEP 4: COMPUTE ADAPTIVE WEIGHTS ---
    w = 1 ./ max(abs(beta_OLS(:, m)), epsilon_ols).^gamma;
    if any(~isfinite(w))
        error('Adaptive_LASSO:BadWeights', ...
            'Non-finite adaptive weights for %s.', muscleNames{m});
    end
    W(:, m) = w;

    % --- STEP 5: ADAPTIVE LASSO VIA sgl_fit ---
    % Column-scale A to absorb adaptive weights into standard L1.
    A_star = A .* (1 ./ w)';

    optsM = optsCommon;
    optsM.muscleName = muscleNames{m};
    optsM.c = cVec(m);
    optsM.rho = rhoVec(m);
    optsM.beta0 = zeros(p, 1);
    optsM.b0 = 0;
    optsM.doWarmStart = true;

    % Optional warm start from OLS (in A* space):
    % optsM.beta0 = beta_OLS(:, m) .* w;
    % optsM.b0 = b_OLS(m);

    [betaOrig_star, intercept_star, ~, stats_ada{m}] = ...
        sgl_fit(A_star, Y(:, m), lambda1Vec(m), lambda2Vec(m), optsM);

    Beta_adaptive(:, m) = betaOrig_star ./ w;
    Intercept_adaptive(m) = intercept_star;

    YPred_adaptive(:, m) = max(A * Beta_adaptive(:, m) + Intercept_adaptive(m), cVec(m));
    RMSE_adaptive(m) = sqrt(mean((Y(:, m) - YPred_adaptive(:, m)).^2));

    S{m} = find(abs(Beta_adaptive(:, m)) > nzTol);

    % --- STEP 6: POST-SELECTION OLS VIA sgl_fit ---
    if isempty(S{m})
        warning('Adaptive_LASSO:EmptySupport', ...
            'Support set for %s is empty. Skipping post-selection OLS.', ...
            muscleNames{m});
        Beta_refit(:, m) = zeros(p, 1);
        Intercept_refit(m) = 0;
        YPred_refit(:, m) = cVec(m) * ones(nSamples, 1);
        RMSE_refit(m) = sqrt(mean((Y(:, m) - YPred_refit(:, m)).^2));
        stats_refit{m} = struct('exitType', 'skipped_empty_support');
        continue;
    end

    A_refit = A;
    nonS = setdiff(1:p, S{m});
    A_refit(:, nonS) = 0;

    optsRefit = optsCommon;
    optsRefit.muscleName = [muscleNames{m}, '_refit'];
    optsRefit.c = cVec(m);
    optsRefit.rho = rhoVec(m);
    optsRefit.doWarmStart = true;
    optsRefit.beta0 = Beta_adaptive(:, m);
    optsRefit.b0 = Intercept_adaptive(m);

    [Beta_refit(:, m), Intercept_refit(m), ~, stats_refit{m}] = ...
        sgl_fit(A_refit, Y(:, m), 0, 0, optsRefit);

    YPred_refit(:, m) = max(A * Beta_refit(:, m) + Intercept_refit(m), cVec(m));
    RMSE_refit(m) = sqrt(mean((Y(:, m) - YPred_refit(:, m)).^2));
end

%% ===== STEP 7: SUMMARY METRICS =====
fprintf('\n===== MULTI-MUSCLE ADAPTIVE LASSO SUMMARY =====\n');
for m = 1:numTargets
    fprintf(['[%s] OLS_nnz=%d | Ada_nnz=%d | Refit_nnz=%d | ' ...
        'RMSE_ada=%.4g | RMSE_refit=%.4g | exit=%s\n'], ...
        muscleNames{m}, ...
        nnz(abs(beta_OLS(:, m)) > nzTol), ...
        nnz(abs(Beta_adaptive(:, m)) > nzTol), ...
        nnz(abs(Beta_refit(:, m)) > nzTol), ...
        RMSE_adaptive(m), RMSE_refit(m), stats_ada{m}.exitType);
end

%% ===== STEP 8: ANALYZE SPARSE SOLUTION =====
indexToFeatureName = buildInverseFeatureMap(featureToIndexMap, size(X, 2));

for m = 1:numTargets
    nzIdx = find(abs(Beta_adaptive(:, m)) > nzTol);
    fprintf('\nAdaptive LASSO nonzero controls for %s (nnz = %.0f):\n', muscleNames{m}, numel(nzIdx));
    if isempty(nzIdx)
        fprintf('(none)\n');
    else
        for k = 1:numel(nzIdx)
            i = nzIdx(k);
            featureName = indexToFeatureName{i};
            fprintf('  %s -> %s: %.5g  (index %d)\n', featureName, muscleNames{m}, Beta_adaptive(i, m), i);
        end
    end
    fprintf('prestimulation (%s, adaptive): %.5g\n', muscleNames{m}, Intercept_adaptive(m));

    nzIdxRefit = find(abs(Beta_refit(:, m)) > nzTol);
    fprintf('\nPost-selection OLS nonzero controls for %s (nnz = %.0f):\n', muscleNames{m}, numel(nzIdxRefit));
    if isempty(nzIdxRefit)
        fprintf('(none)\n');
    else
        for k = 1:numel(nzIdxRefit)
            i = nzIdxRefit(k);
            featureName = indexToFeatureName{i};
            fprintf('  %s -> %s: %.5g  (index %d)\n', featureName, muscleNames{m}, Beta_refit(i, m), i);
        end
    end
    fprintf('prestimulation (%s, refit): %.5g\n', muscleNames{m}, Intercept_refit(m));
end

%% ===== STEP 9: VISUALIZE FIT QUALITY =====
for m = 1:numTargets
    figure('Name', sprintf('Adaptive LASSO Fitting Results: %s', muscleNames{m}), 'Color', 'w');
    plot(Y(:, m), 'LineWidth', 1.0); hold on;
    plot(YPred_adaptive(:, m), 'LineWidth', 1.0);
    plot(YPred_refit(:, m), 'LineWidth', 1.0);
    xlim([1 size(Y, 1)]);
    xlabel('#timesteps');
    ylabel('excitation');
    legend('SCONE muscle excitation', 'Adaptive LASSO', 'Post-selection OLS');
    title(sprintf('%s: Adaptive LASSO (RMSE=%.4g) vs Refit (RMSE=%.4g)', ...
        muscleNames{m}, RMSE_adaptive(m), RMSE_refit(m)), 'Interpreter','none');
    grid on;
end

%% ===== COLLECT RESULTS =====
results = struct();
results.muscleNames = muscleNames;
results.A = A;
results.Y = Y;
results.beta_OLS = beta_OLS;
results.b_OLS = b_OLS;
results.W = W;
results.Beta_adaptive = Beta_adaptive;
results.Intercept_adaptive = Intercept_adaptive;
results.YPred_adaptive = YPred_adaptive;
results.RMSE_adaptive = RMSE_adaptive;
results.stats_ada = stats_ada;
results.S = S;
results.Beta_refit = Beta_refit;
results.Intercept_refit = Intercept_refit;
results.YPred_refit = YPred_refit;
results.RMSE_refit = RMSE_refit;
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

% TESTING CHECKLIST:
% [ ] With lambda2 very large: all beta_adaptive should be zero
% [ ] With lambda2 = 0: adaptive LASSO should approximate OLS
% [ ] Post-selection OLS nnz should equal adaptive LASSO nnz
% [ ] Post-selection OLS RMSE should be <= adaptive LASSO RMSE
%     (debiasing reduces shrinkage bias)
% [ ] max(abs(A_refit(:,nonS) * Beta_refit(nonS,m))) should be 0
% [ ] YPred >= c for all predictions (censoring floor)
% [ ] stats_ada{m}.objHist should be non-increasing
