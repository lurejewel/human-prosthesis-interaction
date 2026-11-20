function A = make_knee_selection_controller(stoFile, outDir)
%MAKE_KNEE_SELECTION_CONTROLLER  Build the controller that uses each
%   muscle-group's own lambda-path knee instead of a tolSelect rule.
%
%   A = make_knee_selection_controller(stoFile, outDir)
%
%   For every (muscle, phase group) this takes the Kneedle knee of its own
%   (support size, RMSE/best RMSE) path -- computed on the Pareto-filtered
%   points -- and uses the support at that path index as the structure.  The
%   coefficients are then re-estimated on the FULL data with that frozen mask
%   (mask_override), exactly like the consensus controllers.
%
%   This is an ALTERNATIVE selection rule to `tolSelect`; it is provided for
%   comparison / inspection, not as the default pipeline.
%
%   Outputs
%   -------
%   results/lasso_controller_result_moco_consensus_knee.mat   (the controller)
%   results/lasso_fit_quality_knee.png                        (30 s traces)
%   results/lasso_fit_quality_zoom_knee.png                   (3 s traces)
%   results/lasso_fit_rmse_comparison_knee.png                (RMSE bars)
%   A  : summary struct (mask, total parameters, per-muscle RMSE, knee nnz)
%
%   See also plot_lasso_fit_quality, plot_lasso_sparsity_tradeoff.

if nargin < 1 || isempty(stoFile)
    here = fileparts(mfilename('fullpath'));
    stoFile = fullfile(here, 'moco_dataset', 'moco_lasso_dataset.sto');
end
here = fileparts(mfilename('fullpath'));
repoRoot = fileparts(here);
resultsDir = fullfile(repoRoot, 'results');
if nargin < 2 || isempty(outDir)
    outDir = resultsDir;
end

S = load(fullfile(resultsDir, 'lasso_stability_diagnostics.mat'), 'R');
R = S.R;
ref = R.fitReference;
phs = {'stance', 'swing'};
nMus = numel(ref.muscleNames);
p = 26;

mask = cell(1, 2);
kneeNnz = zeros(2, nMus);
for g = 1:2
    ph = phs{g};
    mask{g} = false(p, nMus);
    for m = 1:nMus
        path = ref.LambdaPath.(ph){m};
        x = path.supportSize(:);
        y = path.trainErrorRefit(:) / max(path.minTrainErrorRefit, eps);
        k = knee_index(x, y);
        mask{g}(path.supportPath{k}, m) = true;
        kneeNnz(g, m) = numel(path.supportPath{k});
    end
end

o = ref.options;
o.maskOverride = mask;
o.verbose = false;
o.storeLambdaPath = false;
D = load_lasso_design_matrices(stoFile, o);
[fitK, lassoK] = fit_adaptive_lasso_core(D, o);

totalLen = nnz(mask{1}) + 9 + nnz(mask{2}) + 9;

% ---- save the controller -------------------------------------------------
lasso = lassoK; %#ok<NASGU>
ctrlFile = fullfile(outDir, 'lasso_controller_result_moco_consensus_knee.mat');
save(ctrlFile, 'lasso');

% ---- diagnostics struct for the plotting function ------------------------
Rk = struct();
Rk.fitFinal = fitK;
Rk.totalLen = totalLen;
Rk.fitReference = ref;
Rk.referenceTotalLen = R.referenceTotalLen;
Rk.mask = mask;
Rk.kneeNnz = kneeNnz;
R = Rk; %#ok<NASGU>
diagFile = fullfile(outDir, '_knee_diag.mat');
save(diagFile, 'R');

% ---- figures (same layout as the consensus fit-quality figures) ----------
dsMat = fullfile(here, 'moco_dataset', 'moco_lasso_dataset.mat');
plot_lasso_fit_quality(diagFile, dsMat, outDir, [], '_knee');

A = struct();
A.controllerFile = ctrlFile;
A.diagFile = diagFile;
A.totalLen = totalLen;
A.kneeNnz = kneeNnz;
A.meanRMSE = mean(fitK.RMSE_refit_full);
A.perMuscleRMSE = fitK.RMSE_refit_full(:).';
A.refMeanRMSE = mean(ref.RMSE_refit_full);

fprintf('\n[knee] stance nnz=%d swing nnz=%d | total params=%d (ref %d)\n', ...
    nnz(mask{1}), nnz(mask{2}), totalLen, Rk.referenceTotalLen);
fprintf('[knee] mean RMSE = %.5f  (unfrozen reference %.5f, %+.1f%%)\n', ...
    A.meanRMSE, A.refMeanRMSE, 100 * (A.meanRMSE - A.refMeanRMSE) / A.refMeanRMSE);
fprintf('[knee] per-muscle RMSE = %s\n', mat2str(A.perMuscleRMSE, 4));
fprintf('[knee] wrote %s\n', ctrlFile);
end

%% ========================================================================
function k = knee_index(x, y)
%KNEE_INDEX  Kneedle knee on one muscle-group path; returns the PATH INDEX.
ok = isfinite(x) & isfinite(y) & x >= 1;   % exclude the zero-support solution
x = x(ok); y = y(ok);
ii = find(ok);
n = numel(x);
if n < 3
    k = ii(end);
    return;
end
% Pareto filter (minimise both), then sort by x
nd = true(n, 1);
for i = 1:n
    for j = 1:n
        if j ~= i && x(j) <= x(i) && y(j) <= y(i) && (x(j) < x(i) || y(j) < y(i))
            nd(i) = false;
            break;
        end
    end
end
x = x(nd); y = y(nd); ii = ii(nd);
[x, is] = sort(x); y = y(is); ii = ii(is);
xn = (x - min(x)) / max(max(x) - min(x), eps);
yn = (y - min(y)) / max(max(y) - min(y), eps);
d = (1 - xn) - yn;
[~, j] = max(d);
k = ii(j);
end
