% Sparse Group LASSO runner script, now deprecated. Run "Adaptive_LASSO.m"
% instead.
% ========================================================================
% PURPOSE:
%   Fit a lower-censored sparse-group-lasso model to predict muscle
%   excitation from biomechanical features using FISTA optimization.
%
% WORKFLOW:
%   1. Load muscle excitation (y) and biomechanical features (X) from SCONE
%   2. Preprocess: apply 1-step delay, remove small observations
%   3. Define regularization: lambda1 (group L2) and lambda2 (L1)
%   4. Configure optimization: FISTA, centering, scaling, censor threshold
%   5. Call sgl_fit() to solve: min ||max(y-Aβ-b,0)||²/2 + λ₁∑||β_g||₂
%   6. Report nonzero connections in sparse solution
%   7. Visualize fit quality: predictions vs SCONE outputs
% ========================================================================
%
% TODO:
% - maskMatrix需要根据肌肉特性，导入先验（重读经典论文，寻找先验依据）
% - 研究SGL代码，看是否需要针对性改造（权重？惩罚？）√（改造后代码为Adaptive_LASSO.m）
% - （soleus_r成功后）切分数据（stance-swing），对其他肌肉的控制律进行验证
% - （所有肌肉成功后）进行正动力学验证
% - （正动力学成功后）进行代理模型研究（对应Annuals Review of Biomedical Engineering中的大脑预测模型）
% - （用于论文）lambda1 lambda2的灵敏度分析

% close all; 
clear; clc;

%% ===== STEP 1: DATA LOADING & PREPROCESSING =====
% Load three muscle excitations and one shared biomechanical feature matrix.
stoFilePath = 'UN.sto';
muscleNames = {'soleus_r', 'tib_ant_r', 'gastroc_r'};
numTargets = numel(muscleNames);

delaySteps = 1;
excitationThreshold = 0.01;

% Per-muscle penalties (independent tuning).
lambda1Vec = [0, 0, 0];
lambda2Vec = [20, 20, 20];

% Extract base muscle to initialize shared X and M.
[yAllBase, XAll, featureToIndexMap, M, phase] = extractStoMuscleFeatures(stoFilePath, muscleNames{1});
[yBase, X] = preprocessExcitationAndFeatures(yAllBase, XAll, delaySteps, excitationThreshold);

Y = zeros(numel(yBase), numTargets); 
Y(:, 1) = yBase;

% Load remaining muscles and verify shared feature/mask assumptions.
for m = 2:numTargets
    [yAllM, XAllM, featureToIndexMapM, MM, ~] = extractStoMuscleFeatures(stoFilePath, muscleNames{m});

    if ~isequal(size(XAllM), size(XAll)) || max(abs(XAllM(:) - XAll(:))) > 1e-12
        error('Sparse_Group_LASSO:FeatureMismatch', ...
            'Feature matrix for %s is not identical to base muscle %s.', muscleNames{m}, muscleNames{1});
    end

    if ~isequal(MM, M)
        error('Sparse_Group_LASSO:MaskMismatch', ...
            'Mask matrix for %s is not identical to base muscle %s.', muscleNames{m}, muscleNames{1});
    end

    if ~isequal(sort(keys(featureToIndexMapM)), sort(keys(featureToIndexMap)))
        error('Sparse_Group_LASSO:FeatureMapMismatch', ...
            'Feature map for %s differs from base muscle %s.', muscleNames{m}, muscleNames{1});
    end

    [yM, XM] = preprocessExcitationAndFeatures(yAllM, XAllM, delaySteps, excitationThreshold);
    if ~isequal(size(XM), size(X)) || max(abs(XM(:) - X(:))) > 1e-12
        error('Sparse_Group_LASSO:PreprocessMismatch', ...
            'Preprocessed feature matrix for %s is not identical to base muscle %s.', muscleNames{m}, muscleNames{1});
    end

    Y(:, m) = yM;
end

fprintf('Loaded %d muscles with shared feature matrix: X size = [%d, %d], Y size = [%d, %d].\n', ...
    numTargets, size(X, 1), size(X, 2), size(Y, 1), size(Y, 2));

%% ===== STEP 2: WARM START (optional, from prior run or domain knowledge) =====
% If you have a good initial guess from a previous run, warm-start can speed convergence.
% Specify in ORIGINAL space (not standardized). sgl_fit() auto-converts internally.

% Build per-muscle options with independent censoring and warm starts.
optsCommon = struct();
optsCommon.maxIter = 20000;
optsCommon.tol = 1e-10;
optsCommon.verbose = true;
optsCommon.doWarmStart = true;
optsCommon.doCenter = true;
optsCommon.doScale = true;
optsCommon.scaleType = 'std';
optsCommon.backtrackBeta = 0.5;
optsCommon.L0 = 1000;
optsCommon.tolCensor = 1e-12;

% Independent per-muscle censoring and penalty configuration.
cVec = [0.01, 0.01, 0.01];
rhoVec = [50, 50, 50];

% Optional warm starts in ORIGINAL space.
warmStartBeta = cell(numTargets, 1);
warmStartB = zeros(numTargets, 1);
for m = 1:numTargets
    warmStartBeta{m} = zeros(size(M, 2), 1);
end

% 这里需要修改：
% 通过文件实现warm start读取
% 分步态阶段读取（站立相和摆动相不同）
% soleus_r: 
warmStartBeta{1}(featureToIndexMap('soleus_r.fiber_length_norm')) = 0.681;
warmStartBeta{1}(featureToIndexMap('soleus_r.mtu_force_norm')) = 0.920;
warmStartBeta{1}(featureToIndexMap('tib_ant_r.mtu_force_norm')) = -0.704;
warmStartB(1) = -0.610;

% tib_ant_r
warmStartBeta{2}(featureToIndexMap('tib_ant_r.fiber_length_norm')) = 0.685;
warmStartBeta{2}(featureToIndexMap('tib_ant_r.mtu_force_norm')) = 0.657;
warmStartBeta{2}(featureToIndexMap('soleus_r.fiber_length_norm')) = -0.050;
warmStartBeta{2}(featureToIndexMap('soleus_r.mtu_force_norm')) = -0.986;
warmStartBeta{2}(featureToIndexMap('gastroc_r.fiber_length_norm')) = -0.029;
warmStartBeta{2}(featureToIndexMap('gastroc_r.mtu_force_norm')) = -0.771;
warmStartB(2) = -0.314;

% gastroc_r
warmStartBeta{3}(featureToIndexMap('gastroc_r.fiber_length_norm')) = 0.452;
warmStartBeta{3}(featureToIndexMap('gastroc_r.mtu_force_norm')) = 0.686;
warmStartBeta{3}(featureToIndexMap('tib_ant_r.mtu_force_norm')) = -0.819;
warmStartB(3) = -0.308;


%% ===== STEP 4: BUILD FEATURE MATRIX & RUN OPTIMIZATION =====
% Build one shared A = X * M and solve each muscle independently.
A = X * M; % feature matrix X (augmented), masked with the sparse matrix M

numCoeffs = size(A, 2);
Beta = zeros(numCoeffs, numTargets);
Intercept = zeros(numTargets, 1);
BetaStd = zeros(numCoeffs, numTargets);
Stats = cell(numTargets, 1);
YPred = zeros(size(Y));
RMSE = zeros(numTargets, 1);

for m = 1:numTargets
    optsM = optsCommon;
    optsM.muscleName = muscleNames{m};
    optsM.c = cVec(m);
    optsM.rho = rhoVec(m);
    if optsCommon.doWarmStart
        optsM.beta0 = warmStartBeta{m};
        optsM.b0 = warmStartB(m);
    else
        optsM.beta0 = [];
        optsM.b0 = 0;
    end

    fprintf('\n===== Fitting target muscle: %s (%d/%d) =====\n', muscleNames{m}, m, numTargets);
    [Beta(:, m), Intercept(m), BetaStd(:, m), Stats{m}] = sgl_fit(A, Y(:, m), lambda1Vec(m), lambda2Vec(m), optsM);
    YPred(:, m) = max(A * Beta(:, m) + Intercept(m), optsM.c);
    RMSE(m) = sqrt(mean((Y(:, m) - YPred(:, m)).^2));
end

% test
% yPred = max(0.01, X * M * beta + intercept);
% figure, plot(y), hold on, plot(yPred)
% sqrt(mean((y - yPred).^2))
% beta(:) = 0; beta(26) = 0.68247; beta(27) = -0.58974; beta(37) = 0.44245; intercept = 0.097071;
% yPredAllSwing = XAllSwing * M * beta + intercept;
% figure, plot(yAllSwing), hold on, plot(yPredAllSwing)
% sqrt(mean((yAllSwing - yPredAllSwing).^2))


%% ===== STEP 5: SUMMARY METRICS =====
fprintf('\n===== MULTI-MUSCLE SGL SUMMARY =====\n');
for m = 1:numTargets
    fprintf('[%s] objective = %.12g | RMSE = %.6g | iters = %d | exit = %s\n', ...
        muscleNames{m}, Stats{m}.objHist(end), RMSE(m), Stats{m}.iters, Stats{m}.exitType);
end

%% ===== STEP 6: ANALYZE SPARSE SOLUTION =====
% Extract and report nonzero connections for each muscle target.
indexToFeatureName = buildInverseFeatureMap(featureToIndexMap, size(X, 2));
nzTol = 1e-5;

for m = 1:numTargets
    nzIdx = find(abs(Beta(:, m)) > nzTol);
    fprintf('\nNonzero muscle reflex controls for %s (nnz = %.0f):\n', muscleNames{m}, numel(nzIdx));
    if isempty(nzIdx)
        fprintf('(none)\n');
    else
        for k = 1:numel(nzIdx)
            i = nzIdx(k);
            featureName = indexToFeatureName{i};
            fprintf('  %s -> %s: %.5g  (index %d)\n', featureName, muscleNames{m}, Beta(i, m), i);
        end
    end

    fprintf('prestimulation (%s): %.5g\n', muscleNames{m}, Intercept(m));
end

%% ===== STEP 7: VISUALIZE FIT QUALITY =====
% Plot actual vs predicted excitation for each muscle.
figure('Name', 'SGL Fitting Results for Multiple Muscles', 'Color', 'w');
for m = 1:numTargets
    subplot(numTargets, 1, m);
    plot(Y(:, m), 'LineWidth', 1.0); hold on;
    plot(YPred(:, m), 'LineWidth', 1.0);
    xlim([0 size(Y, 1)]);
    xlabel('#timesteps');
    ylabel('excitation');
    legend('SCONE muscle excitation', 'SGL fitting results');
    title(sprintf('%s: SGL fitting vs SCONE (RMSE = %.4g)', muscleNames{m}, RMSE(m)), Interpreter="none");
    grid on;
end

% Optional exported variables for downstream analysis.
results = struct();
results.muscleNames = muscleNames;
results.A = A;
results.Y = Y;
results.Beta = Beta;
results.Intercept = Intercept;
results.BetaStd = BetaStd;
results.YPred = YPred;
results.RMSE = RMSE;
results.Stats = Stats;

% INTERPRETATION GUIDE:
%   Good fit (red ≈ blue):  model captures activation dynamics well
%   Peaks missed:           reduce λ₁/λ₂ for less sparsity, or increase ρ
%   High-frequency noise:   increase λ₁/λ₂ for more regularization

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