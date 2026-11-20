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
%                .medianNnz{g}/.medianRelErr{g} (per-muscle median curves),
%                .nnzMat/.errMat (nGroups x nMus x nPath),
%                .totalParamsK / .allErrK / .aggRelErrK (aggregate curve:
%                one global path index k applied to all muscles),
%                .knee{g} (per-muscle) and .kneeAgg (total-parameter curve),
%                .nPoints*.nDominated* (Pareto dominance counts)
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
    % Median curve across the 9 muscles.  Paths are aligned on the alpha grid
    % but a muscle with lambda2Max == 0 gets a degenerate 1-point path, so do
    % NOT assume muscle 1's length applies to all muscles: use the longest
    % path and let shorter paths contribute NaN at the missing indices.
    % Median then uses the muscles that still have a point at that index.
    pathLens = zeros(1, nMus);
    for m = 1:nMus
        pathLens(m) = numel(fitPath.LambdaPath.(ph){m}.supportSize);
    end
    nP = max(pathLens);
    medN = nan(nP, 1);
    medE = nan(nP, 1);
    for k = 1:nP
        vn = nan(1, nMus);
        ve = nan(1, nMus);
        for m = 1:nMus
            if k > pathLens(m)
                continue;   % shorter/truncated path: NaN for this index
            end
            p = fitPath.LambdaPath.(ph){m};
            vn(m) = p.supportSize(k);
            ve(m) = p.trainErrorRefit(k) / max(p.minTrainErrorRefit, eps);
        end
        ok = isfinite(vn) & isfinite(ve);
        if any(ok)
            medN(k) = median(vn(ok));
            medE(k) = median(ve(ok));
        end
    end
    medianNnz{g} = medN;
    medianRelErr{g} = medE;
end

% ---- per-muscle path matrices (NaN where a degenerate path is shorter) ----
nP = max(cellfun(@numel, medianNnz));
nnzMat = nan(2, nMus, nP);     % support size of muscle m, group g, path index k
errMat = nan(2, nMus, nP);     % RMSE / best RMSE, same indexing
for g = 1:2
    ph = groupLabels{g};
    for m = 1:nMus
        p = fitPath.LambdaPath.(ph){m};
        n = numel(p.supportSize);
        nnzMat(g, m, 1:n) = p.supportSize(:).';
        errMat(g, m, 1:n) = (p.trainErrorRefit(:) / ...
            max(p.minTrainErrorRefit, eps)).';
    end
end

% ---- aggregate trade-off: ONE path index k applied to all muscles --------
% x = total CMA-ES parameters if the path is stopped at index k everywhere
%     (same k for every muscle-group; this is what a global tolerance does);
% y = median (and mean) relative RMSE over the 18 muscle-group curves.
totalParamsK = squeeze(sum(sum(nnzMat, 2, 'omitnan'), 1, 'omitnan')) + 18;
allErrK = reshape(permute(errMat, [2 1 3]), nMus * 2, nP);   % 18 x nP
aggRelErrK = median(allErrK, 1, 'omitnan').';
aggRelErrMeanK = mean(allErrK, 1, 'omitnan').';

A.path.targets = targets;
A.path.nnzAtTarget = nnzAtTarget;      % target x group x muscle
A.path.nnzSelected = nnzSelected;      % group x muscle
A.path.medianNnz = medianNnz;
A.path.medianRelErr = medianRelErr;
A.path.nnzMat = nnzMat;
A.path.errMat = errMat;
A.path.totalParamsK = totalParamsK;
A.path.allErrK = allErrK;
A.path.aggRelErrK = aggRelErrK;
A.path.aggRelErrMeanK = aggRelErrMeanK;

% ---- Pareto structure: which path points are dominated? ------------------
% A point (nnz, RMSE) is dominated if another lambda gives both no more
% nonzeros AND no worse RMSE.  The lambda path is a *trajectory*, not a
% Pareto front, so dominated points do exist (support-size plateaus, or a
% later lambda that is both denser and worse).
fprintf('\n  Pareto structure of the lambda paths (minimise nnz and RMSE):\n');
A.path.nDominated = zeros(1, 2);
A.path.nPoints = zeros(1, 2);
for g = 1:2
    nd = pareto_flags(medianNnz{g}, medianRelErr{g});
    A.path.nPoints(g) = sum(isfinite(medianNnz{g}) & isfinite(medianRelErr{g}));
    A.path.nDominated(g) = A.path.nPoints(g) - nnz(nd);
    fprintf('    %-6s median curve: %d/%d points non-dominated\n', ...
        groupLabels{g}, nnz(nd), A.path.nPoints(g));
end

nDomMuscle = 0; nPtMuscle = 0;
for g = 1:2
    for m = 1:nMus
        x = squeeze(nnzMat(g, m, :));
        y = squeeze(errMat(g, m, :));
        ok = isfinite(x) & isfinite(y);
        nd = pareto_flags(x(ok), y(ok));
        nPtMuscle = nPtMuscle + nnz(ok);
        nDomMuscle = nDomMuscle + (nnz(ok) - nnz(nd));
    end
end
A.path.nPointsMuscle = nPtMuscle;
A.path.nDominatedMuscle = nDomMuscle;
fprintf('    all 18 per-muscle paths: %d/%d points non-dominated\n', ...
    nPtMuscle - nDomMuscle, nPtMuscle);

% ---- knees --------------------------------------------------------------
% (a) per phase group, on the per-muscle median curve (x = nnz per muscle)
A.path.knee = cell(1, 2);
fprintf('\n  knee of the per-muscle median curve (Kneedle):\n');
for g = 1:2
    [kx, ky] = knee_point(medianNnz{g}, medianRelErr{g});
    A.path.knee{g} = struct('nnz', kx, 'relErr', ky, 'tolSelect', ky - 1);
    fprintf('    %-6s : nnz=%.1f, RMSE/RMSEbest=%.3f -> implied tolSelect=%.2f\n', ...
        groupLabels{g}, kx, ky, ky - 1);
end

% (b) aggregate curve (x = total CMA-ES parameters, ONE global tolerance)
ndAgg = pareto_flags(totalParamsK, aggRelErrK);
A.path.nPointsAgg = numel(totalParamsK);
A.path.nDominatedAgg = numel(totalParamsK) - nnz(ndAgg);
[kx, ky] = knee_point(totalParamsK, aggRelErrK);
A.path.kneeAgg = struct('params', kx, 'relErr', ky, 'tolSelect', ky - 1, ...
    'nDominated', A.path.nDominatedAgg);
fprintf(['    aggregate (total params): %d/%d non-dominated, knee at ' ...
    '%d params (+%.1f%% RMSE, equivalent global tolSelect~%.2f)\n'], ...
    nnz(ndAgg), numel(totalParamsK), kx, 100 * (ky - 1), ky - 1);

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
    [fitF, ~] = fit_adaptive_lasso_core(D, o);
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
% Panel 1 uses the per-muscle-group support size on x: the lambda is selected
% independently for every muscle-group, so a "total parameters" axis would mix
% 18 separate decisions into a hypothetical global setting.
figFile = fullfile(resultsDir, 'lasso_sparsity_tradeoff.png');
plot_lasso_sparsity_tradeoff(A, R.R, figFile);

%% ===== save ============================================================
matFile = fullfile(resultsDir, 'lasso_sparsity_analysis.mat');
A.files = struct('mat', matFile, 'figure', figFile);
save(matFile, 'A', '-v7.3');
fprintf('\n[sparsity] saved %s and %s\n', matFile, figFile);
end

%% ========================================================================
function nd = pareto_flags(x, y)
%PARETO_FLAGS  Non-dominated flags for a 2-objective MINIMISATION set.
%   A point is dominated when another point has x <= and y <= with at least
%   one strict inequality.  O(n^2); n here is at most a few hundred.
x = x(:); y = y(:);
n = numel(x);
nd = true(n, 1);
for i = 1:n
    if ~isfinite(x(i)) || ~isfinite(y(i))
        nd(i) = false;
        continue;
    end
    for j = 1:n
        if j == i || ~isfinite(x(j)) || ~isfinite(y(j))
            continue;
        end
        if x(j) <= x(i) && y(j) <= y(i) && (x(j) < x(i) || y(j) < y(i))
            nd(i) = false;
            break;
        end
    end
end
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
