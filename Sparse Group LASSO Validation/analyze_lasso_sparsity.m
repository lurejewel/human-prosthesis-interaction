function A = analyze_lasso_sparsity(stoFile, options)
%ANALYZE_LASSO_SPARSITY  How much accuracy does each extra nonzero buy?
%
%   A = analyze_lasso_sparsity(stoFile, options)
%
%   The reflex controller should be CONCISE (the 0-1 structure is the point,
%   not a precise fit of the measured excitation).  This tool quantifies the
%   sparsity/accuracy trade-off of the current identification so that
%   `tolSelect` and the consensus pi-threshold can be chosen on evidence
%   rather than by taste.
%
%   Three analyses
%   -------------
%   1. lambda-path marginal return: for every (muscle, phase group) the path
%      already contains 60 (support size, refit RMSE) pairs.  We ask: how many
%      nonzeros are needed to stay within +2/+5/+10/+20/+50% of the best
%      achievable RMSE?  If the curve is flat, the 5% rule is over-selecting.
%   2. pi sweep: threshold the consensus selection frequency at several levels,
%      refit the coefficients on the FULL data and (for a generalization check)
%      on two run-folds, and report (parameters, train RMSE, held-out RMSE).
%      NOTE: the mask comes from the full ensemble, so the held-out number is
%      a coefficient-generalization check, not an unbiased structure-selection
%      CV; it is meant to be compared ACROSS thresholds.
%   3. tolSelect sweep: how strongly the base selector's 5% tolerance drives
%      the support size (full-data fits only).
%
%   Output A
%   --------
%   .path        .targets, .nnzAtTarget (nTarget x nGroups), .nnzSelected,
%                .medianNnz{g}/.medianRelErr{g} (60-point median curves),
%                .knee{g} = .nnz/.relErr/.tolSelect
%   .piSweep     .thr, .params, .nnz,
%                .trainRMSE  refit on ALL data, evaluated on ALL data
%                .trainCV    in-fold train (fit on a fold, evaluated on it)
%                .cvRMSE     out-of-fold test (fit on the other fold)
%                .gap        cvRMSE - trainCV (overfitting gap at equal sizes)
%                .knee       .params/.cvRMSE/.pi
%   .tolSweep    .tolSelect, .params, .nnz, .trainRMSE
%   .reference   unfrozen full-data fit (params, RMSE)
%   .files       where the .mat / figure were written
%
%   See also fit_adaptive_lasso_core, lasso_stability_from_design.

if nargin < 1 || isempty(stoFile)
    here = fileparts(mfilename('fullpath'));
    stoFile = fullfile(here, 'moco_dataset', 'moco_lasso_dataset.sto');
end
if nargin < 2 || isempty(options)
    options = struct();
end
here = fileparts(mfilename('fullpath'));
repoRoot = fileparts(here);
resultsDir = fullfile(repoRoot, 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

fitOpts = base_fit_options(options);
fitOpts.verbose = false;
fitOpts.storeLambdaPath = true;

% Match the production pipeline exactly: sample weights (seam/invalid frames
% get weight 0) and the preprocessing metadata come from the dataset .mat.
matFile = fullfile(fileparts(stoFile), 'moco_lasso_dataset.mat');
if isfile(matFile)
    Sm = load(matFile, 'ds');
    fitOpts.sampleWeight = Sm.ds.sampleWeight(:);
    fitOpts.preprocessing = Sm.ds.meta;
end

fprintf('[sparsity] loading %s ...\n', stoFile);
D = load_lasso_design_matrices(stoFile, fitOpts);
groupLabels = fitOpts.groupLabels;
nMus = numel(D.muscleNames);

%% ===== 1. LAMBDA-PATH MARGINAL RETURN ==================================
fprintf('[sparsity] full fit with lambda paths (tolSelect=%.2f) ...\n', ...
    fitOpts.lambdaPathConfig.tolSelect);
[fitPath, lassoPath] = fit_adaptive_lasso_core(D, fitOpts);

targets = [1.02 1.05 1.10 1.20 1.50 2.00];
nnzAtTarget = nan(numel(targets), 2, nMus);
nnzSelected = nan(2, nMus);
medianNnz = cell(1, 2);
medianRelErr = cell(1, 2);
for g = 1:2
    ph = groupLabels{g};
    for m = 1:nMus
        p = fitPath.LambdaPath.(ph){m};
        errRel = p.trainErrorRefit(:) / max(p.minTrainErrorRefit, eps);
        nnzPath = p.supportSize(:);
        nnzSelected(g, m) = numel(fitPath.S.(ph){m});
        for t = 1:numel(targets)
            k = find(errRel <= targets(t), 1, 'first');
            if ~isempty(k)
                nnzAtTarget(t, g, m) = nnzPath(k);
            end
        end
    end
    % median curve across the 9 muscles (same 60 lambda indices)
    nP = numel(fitPath.LambdaPath.(ph){1}.supportSize);
    medN = zeros(nP, 1);
    medE = zeros(nP, 1);
    for k = 1:nP
        vn = zeros(1, nMus);
        ve = zeros(1, nMus);
        for m = 1:nMus
            p = fitPath.LambdaPath.(ph){m};
            vn(m) = p.supportSize(k);
            ve(m) = p.trainErrorRefit(k) / max(p.minTrainErrorRefit, eps);
        end
        medN(k) = median(vn);
        medE(k) = median(ve);
    end
    medianNnz{g} = medN;
    medianRelErr{g} = medE;
end
A.path.targets = targets;
A.path.nnzAtTarget = nnzAtTarget;      % target x group x muscle
A.path.nnzSelected = nnzSelected;      % group x muscle
A.path.medianNnz = medianNnz;
A.path.medianRelErr = medianRelErr;

% ---- knee point of the sparsity/accuracy curve (per muscle-group curve) ---
A.path.knee = cell(1, 2);
fprintf('\n  knee of the median path curve (Pareto knee, Kneedle criterion):\n');
for g = 1:2
    [kx, ky, kk] = knee_point(medianNnz{g}, medianRelErr{g});
    A.path.knee{g} = struct('nnz', kx, 'relErr', ky, 'tolSelect', ky - 1);
    fprintf('    %-6s : nnz=%.1f, RMSE/RMSEbest=%.3f -> implied tolSelect=%.2f\n', ...
        groupLabels{g}, kx, ky, ky - 1);
end

fprintf('\n===== 1. LAMBDA-PATH MARGINAL RETURN (mean over 9 muscles) =====\n');
fprintf('  tolerance to best RMSE |  stance nnz |  swing nnz\n');
for t = 1:numel(targets)
    fprintf('          +%5.1f%%          | %10.1f | %10.1f\n', ...
        100 * (targets(t) - 1), mean(nnzAtTarget(t, 1, :)), mean(nnzAtTarget(t, 2, :)));
end
fprintf('  (current tolSelect=%.2f selects: stance %.1f, swing %.1f)\n', ...
    fitOpts.lambdaPathConfig.tolSelect, mean(nnzSelected(1, :)), mean(nnzSelected(2, :)));

%% ===== 2. PI SWEEP (with held-out check) ===============================
diagFile = fullfile(resultsDir, 'lasso_stability_diagnostics.mat');
if ~isfile(diagFile)
    error('analyze_lasso_sparsity:NoDiagnostics', ...
        'Run Adaptive_LASSO_MocoInverse_Stability first (%s missing).', diagFile);
end
R = load(diagFile, 'R');
C = R.R.consensus;

% two run-folds (odd / even strides within each phase group)
sA = fold_idx(D.phaseIdxMap.stance, 1);
wA = fold_idx(D.phaseIdxMap.swing, 1);
sB = fold_idx(D.phaseIdxMap.stance, 2);
wB = fold_idx(D.phaseIdxMap.swing, 2);
DA = pack_design(D, [sA; wA], numel(sA));
DB = pack_design(D, [sB; wB], numel(sB));

thrGrid = [0.30 0.40 0.50 0.60 0.70 0.80 0.85 0.90 0.95 1.00];
piThr = struct('thr', thrGrid, 'params', nan(size(thrGrid)), ...
    'nnz', nan(numel(thrGrid), 2), 'trainRMSE', nan(size(thrGrid)), ...
    'trainCV', nan(size(thrGrid)), 'cvRMSE', nan(size(thrGrid)), ...
    'gap', nan(size(thrGrid)));

fprintf('\n===== 2. PI THRESHOLD SWEEP =====\n');
fprintf('  pi>= | params | stance | swing | train(full) | train(fold) | held-out | gap\n');
o = fitOpts;
o.storeLambdaPath = false;
for i = 1:numel(thrGrid)
    mask = cell(1, 2);
    for g = 1:2
        mask{g} = (C.Pi{g} >= thrGrid(i) - 1e-12);
    end
    o.maskOverride = mask;
    [fitF, lassoF] = fit_adaptive_lasso_core(D, o);
    piThr.nnz(i, :) = [nnz(mask{1}), nnz(mask{2})];
    piThr.params(i) = sum(piThr.nnz(i, :)) + 18;
    piThr.trainRMSE(i) = mean(fitF.RMSE_refit_full);

    % Proper 2-fold CV: within each fold compare in-fold train with
    % out-of-fold test, so the train/test numbers share the same sample size.
    [~, lassoA] = fit_adaptive_lasso_core(DA, o);   % fit on fold A
    [~, lassoB] = fit_adaptive_lasso_core(DB, o);   % fit on fold B
    trainA = eval_lasso(lassoA, DA, fitOpts.c);
    trainB = eval_lasso(lassoB, DB, fitOpts.c);
    testA = eval_lasso(lassoB, DA, fitOpts.c);      % fit B -> test on A
    testB = eval_lasso(lassoA, DB, fitOpts.c);      % fit A -> test on B
    piThr.trainCV(i) = mean([trainA, trainB]);
    piThr.cvRMSE(i) = mean([testA, testB]);
    piThr.gap(i) = piThr.cvRMSE(i) - piThr.trainCV(i);

    fprintf('  %.2f | %6d | %6d | %5d | %11.5f | %11.5f | %9.5f | %+.1f%%\n', ...
        thrGrid(i), piThr.params(i), piThr.nnz(i, 1), piThr.nnz(i, 2), ...
        piThr.trainRMSE(i), piThr.trainCV(i), piThr.cvRMSE(i), ...
        100 * piThr.gap(i) / max(piThr.trainCV(i), eps));
end
A.piSweep = piThr;

% knee of the (parameters, out-of-fold RMSE) trade-off
[kx, ky] = knee_point(piThr.params(:), piThr.cvRMSE(:));
row = find(piThr.params == kx, 1, 'first');   % knee_point sorts/dedups x
A.piSweep.knee = struct('params', kx, 'cvRMSE', ky, 'pi', piThr.thr(row));
fprintf('\n  knee of the pi sweep: params=%.0f (pi>=%.2f), out-of-fold RMSE=%.5f\n', ...
    kx, A.piSweep.knee.pi, ky);

%% ===== 3. tolSelect SWEEP (full-data fits) ============================
tolGrid = [0.05 0.10 0.20 0.30];
tolSweep = struct('tolSelect', tolGrid, 'params', nan(size(tolGrid)), ...
    'nnz', nan(numel(tolGrid), 2), 'trainRMSE', nan(size(tolGrid)));

fprintf('\n===== 3. tolSelect SWEEP =====\n');
fprintf('  tolSelect | params | stance | swing | trainRMSE\n');
for i = 1:numel(tolGrid)
    if abs(tolGrid(i) - fitOpts.lambdaPathConfig.tolSelect) < 1e-12
        [fitT, lassoT] = deal(fitPath, lassoPath);
    else
        o2 = fitOpts;
        o2.lambdaPathConfig.tolSelect = tolGrid(i);
        o2.storeLambdaPath = false;
        [fitT, lassoT] = fit_adaptive_lasso_core(D, o2);
    end
    tolSweep.nnz(i, :) = [nnz(lassoT.groups(1).mask), nnz(lassoT.groups(2).mask)];
    tolSweep.params(i) = sum(tolSweep.nnz(i, :)) + 18;
    tolSweep.trainRMSE(i) = mean(fitT.RMSE_refit_full);
    fprintf('  %9.2f | %6d | %6d | %5d | %9.5f\n', tolGrid(i), ...
        tolSweep.params(i), tolSweep.nnz(i, 1), tolSweep.nnz(i, 2), ...
        tolSweep.trainRMSE(i));
end
A.tolSweep = tolSweep;

%% ===== reference: unfrozen full-data fit ==============================
A.reference.params = nnz(lassoPath.groups(1).mask) + 9 + ...
    nnz(lassoPath.groups(2).mask) + 9;
A.reference.RMSE = mean(fitPath.RMSE_refit_full);

%% ===== figure ==========================================================
fig = figure('Name', 'LASSO sparsity trade-off', 'Color', 'w', ...
    'Position', [80 80 1300 420]);

subplot(1, 3, 1); hold on; grid on;
hMed = gobjects(1, 2);
for g = 1:2
    ph = groupLabels{g};
    for m = 1:nMus
        p = fitPath.LambdaPath.(ph){m};
        plot(p.supportSize, p.trainErrorRefit / max(p.minTrainErrorRefit, eps), ...
            '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5, ...
            'HandleVisibility', 'off');
    end
    hMed(g) = plot(medianNnz{g}, medianRelErr{g}, '-o', 'LineWidth', 1.6, 'MarkerSize', 3);
    kn = A.path.knee{g};
    plot(kn.nnz, kn.relErr, 'p', 'MarkerSize', 14, 'LineWidth', 1.2, ...
        'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
end
set(gca, 'YScale', 'log');
xlabel('support size (nonzeros, per muscle-group)');
ylabel('RMSE / best RMSE');
title('lambda path: accuracy vs sparsity (knee = *)');
legend(hMed, groupLabels, 'Location', 'northeast');
ylim([0.99 3]);

subplot(1, 3, 2); hold on; grid on;
plot(piThr.params, piThr.trainRMSE, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0);
h1 = plot(piThr.params, piThr.trainCV, '-o', 'LineWidth', 1.6);
h2 = plot(piThr.params, piThr.cvRMSE, '-s', 'LineWidth', 1.6);
xlabel('CMA-ES parameters'); ylabel('mean RMSE');
title('consensus pi threshold (2-fold CV)');
legend([h1 h2], {'in-fold train', 'out-of-fold test'}, 'Location', 'northeast');
kn = A.piSweep.knee;
plot(kn.params, kn.cvRMSE, 'p', 'MarkerSize', 14, 'LineWidth', 1.2, ...
    'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');

subplot(1, 3, 3); hold on; grid on;
plot(tolSweep.params, tolSweep.trainRMSE, '-o', 'LineWidth', 1.6);
xlabel('CMA-ES parameters'); ylabel('mean RMSE');
title('tolSelect (5% = current)');

figFile = fullfile(resultsDir, 'lasso_sparsity_tradeoff.png');
try
    exportgraphics(fig, figFile, 'Resolution', 150);
catch
    print(fig, figFile, '-dpng', '-r150');
end

%% ===== save ============================================================
matFile = fullfile(resultsDir, 'lasso_sparsity_analysis.mat');
A.files = struct('mat', matFile, 'figure', figFile);
save(matFile, 'A', '-v7.3');
fprintf('\n[sparsity] saved %s and %s\n', matFile, figFile);
end

%% ========================================================================
function [kx, ky, idx] = knee_point(x, y)
%KNEE_POINT  Kneedle knee of a monotone trade-off curve.
%   Normalises both axes to [0, 1] and returns the point that lies farthest
%   "below" the chord joining the two extremes -- the classic elbow of an
%   L-shaped (Pareto) front.  Works for a decreasing curve y(x).
x = x(:); y = y(:);
ok = isfinite(x) & isfinite(y);
x = x(ok); y = y(ok);
if numel(x) < 3
    kx = x(end); ky = y(end); idx = numel(x);
    return;
end
% sort by x and keep the best (lowest) y for duplicated x
[xu, ~, ic] = unique(x);
yu = accumarray(ic, y, [], @min);
x = xu; y = yu;
xn = (x - min(x)) / max(max(x) - min(x), eps);
yn = (y - min(y)) / max(max(y) - min(y), eps);
d = (1 - xn) - yn;          % distance below the chord, decreasing curve
[~, idx] = max(d);
kx = x(idx);
ky = y(idx);
end

%% ========================================================================
function o = base_fit_options(user)
o = struct();
o.gamma = 1; o.epsilonOLS = 1e-8; o.nzTol = 1e-5; o.lambda1 = 0;
o.c = 0.01; o.rho = 50; o.delaySteps = 1; o.excitationThreshold = 0.01;
o.robustInitIter = 2; o.robustIter = 2; o.huberDeltaFactor = 1.345;
o.source = 'MocoInverse';
o.groupLabels = {'stance', 'swing'};
o.groupPhases = {[0 1], [2 3 4]};
o.lambdaPathConfig = struct('nLambda', 60, 'lambdaMinRatio', 1e-4, ...
    'tolSelect', 0.05, 'selectionMetric', 'post_refit_train_censored_rmse');
if nargin >= 1 && ~isempty(user)
    fn = fieldnames(user);
    for i = 1:numel(fn)
        o.(fn{i}) = user.(fn{i});
    end
end
end

%% ========================================================================
function runs = split_runs(idxIn)
idxIn = sort(idxIn(:));
d = diff(idxIn);
starts = [1; find(d > 1) + 1];
ends = [find(d > 1); numel(idxIn)];
runs = cell(numel(starts), 1);
for r = 1:numel(starts)
    runs{r} = idxIn(starts(r):ends(r));
end
end

function idxOut = fold_idx(idxIn, parity)
runs = split_runs(idxIn);
if parity == 1
    sel = find(mod(1:numel(runs), 2) == 1);
else
    sel = find(mod(1:numel(runs), 2) == 0);
end
idxOut = sort(vertcat(runs{sel}));
end

function Ds = pack_design(D, allIdx, nS)
Ds = struct();
Ds.X = D.X(allIdx, :);
Ds.X_std = D.X_std(allIdx, :);
Ds.muX = D.muX; Ds.sX = D.sX;
Ds.Y = D.Y(allIdx, :);
Ds.w = D.w(allIdx);
Ds.phaseIdxMap = struct('stance', (1:nS).', 'swing', ((nS + 1):numel(allIdx)).');
Ds.muscleNames = D.muscleNames;
Ds.featureNames = D.featureNames;
Ds.stoFile = D.stoFile;
end

%% ========================================================================
function r = eval_lasso(lasso, D, c)
% Weighted RMSE of the controller's prediction on design D, averaged over
% muscles and phase groups.  lasso holds beta/bias in the ORIGINAL scale.
nGroups = numel(lasso.groups);
rG = nan(1, nGroups);
names = {'stance', 'swing'};
for g = 1:nGroups
    idx = D.phaseIdxMap.(names{g});
    X = D.X(idx, :);
    Y = D.Y(idx, :);
    w = D.w(idx);
    yhat = max(X * lasso.groups(g).beta + lasso.groups(g).bias, c);
    wsum = max(sum(w), eps);
    rmseMus = sqrt(sum(w .* (Y - yhat).^2, 1) ./ wsum);
    rG(g) = mean(rmseMus);
end
r = mean(rG);
end
