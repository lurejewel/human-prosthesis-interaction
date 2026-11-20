function A = plot_controller_fit_quality(controllerFile, stoFile, outDir, zoomWin, tag)
%PLOT_CONTROLLER_FIT_QUALITY  Measured-vs-predicted excitation for ANY
%   LASSO controller .mat (grouped v2), without needing a diagnostics file.
%
%   A = plot_controller_fit_quality(controllerFile, stoFile, outDir, zoomWin, tag)
%
%   The controller already contains beta / bias / mask per group, so the
%   prediction is simply  yhat = max(X*beta + bias, c)  on the phase samples
%   of that group.  This works for the majority consensus, the pi>=0.7
%   consensus, the knee-selected controller, or any hand-made controller that
%   `load_lasso_reflex_controller.m` accepts.
%
%   Outputs (PNG, in outDir; <tag> is appended to the name)
%   -------------------------------------------------------
%   lasso_fit_quality<tag>.png        full 30 s, 9 muscles
%   lasso_fit_quality_zoom<tag>.png   3 s zoom
%   lasso_fit_rmse_comparison<tag>.png per-muscle RMSE bars
%
%   Defaults
%   --------
%   controllerFile = results/lasso_controller_result_moco_consensus_knee.mat
%   stoFile        = Sparse Group LASSO Validation/moco_dataset/moco_lasso_dataset.sto
%   outDir         = results/
%   zoomWin        = 3 s window, 1/3 into the record
%   tag            = '' (use e.g. '_knee' or '_pi70' to keep sets apart)
%
%   Example
%   -------
%   plot_controller_fit_quality( ...
%       'results/lasso_controller_result_moco_consensus_knee.mat', ...
%       'Sparse Group LASSO Validation/moco_dataset/moco_lasso_dataset.sto', ...
%       'results', [], '_knee');
%
%   See also load_lasso_reflex_controller, plot_lasso_fit_quality,
%            make_knee_selection_controller.

here = fileparts(mfilename('fullpath'));
repoRoot = fileparts(here);
if nargin < 1 || isempty(controllerFile)
    controllerFile = fullfile(repoRoot, 'results', ...
        'lasso_controller_result_moco_consensus_knee.mat');
end
if nargin < 2 || isempty(stoFile)
    stoFile = fullfile(here, 'moco_dataset', 'moco_lasso_dataset.sto');
end
if nargin < 3 || isempty(outDir)
    outDir = fullfile(repoRoot, 'results');
end
if nargin < 5 || isempty(tag)
    tag = '';
end

S = load(controllerFile);
if isfield(S, 'lasso')
    lasso = S.lasso;
else
    fn = fieldnames(S);
    if numel(fn) == 1 && isstruct(S.(fn{1}))
        lasso = S.(fn{1});
    else
        error('plot_controller_fit_quality:BadFile', ...
            'No "lasso" variable and not a single struct in %s.', controllerFile);
    end
end
if ~isfield(lasso, 'groups') || isempty(lasso.groups)
    error('plot_controller_fit_quality:NotGrouped', ...
        'Only grouped-v2 controllers are supported (got no .groups).');
end

% design matrix + targets (same preprocessing/weights as the fit)
matFile = fullfile(fileparts(stoFile), 'moco_lasso_dataset.mat');
o = lasso_fit_default_options(struct());
if isfile(matFile)
    Sm = load(matFile, 'ds');
    o.sampleWeight = Sm.ds.sampleWeight(:);
    o.preprocessing = Sm.ds.meta;
end
o.verbose = false;
o.storeLambdaPath = false;
D = load_lasso_design_matrices(stoFile, o);

Y = D.Y;
nMus = numel(D.muscleNames);
Yp = nan(size(Y));
phaseOfGroup = {'stance', 'swing'};
nGroups = numel(lasso.groups);
for g = 1:nGroups
    ph = phaseOfGroup{g};
    idx = D.phaseIdxMap.(ph);
    beta = lasso.groups(g).beta;
    bias = lasso.groups(g).bias(:).';
    if size(beta, 1) ~= size(D.X, 2)
        error('plot_controller_fit_quality:DimMismatch', ...
            'Group %d beta has %d rows but the design has %d features.', ...
            g, size(beta, 1), size(D.X, 2));
    end
    Yp(idx, :) = max(D.X(idx, :) * beta + bias, o.c);
end

wsum = max(sum(D.w), eps);
rmse = sqrt(sum(D.w .* (Y - Yp).^2, 1) ./ wsum);
totalLen = 0;
for g = 1:nGroups
    totalLen = totalLen + nnz(lasso.groups(g).mask) + 9;
end

% package a diagnostics-like struct and reuse the standard plotting routine
fitK = struct();
fitK.Y = Y;
fitK.YPred_refit = Yp;
fitK.muscleNames = D.muscleNames;
fitK.phase = struct('stanceMask', D.phaseMask.stance, ...
    'swingMask', D.phaseMask.swing);
fitK.options = struct('delaySteps', o.delaySteps);
fitK.RMSE_refit_full = rmse(:);
Rk = struct();
Rk.fitFinal = fitK;
Rk.totalLen = totalLen;
Rk.fitReference = [];
Rk.referenceTotalLen = NaN;

tmpFile = [tempname, '.mat'];
R = Rk; %#ok<NASGU>
save(tmpFile, 'R');
cleanup = onCleanup(@() delete_if_exists(tmpFile));
plot_lasso_fit_quality(tmpFile, matFile, outDir, zoomWin, tag);

A = struct();
A.controllerFile = controllerFile;
A.totalLen = totalLen;
A.meanRMSE = mean(rmse);
A.perMuscleRMSE = rmse(:).';

fprintf('[controller-fit] %s\n', controllerFile);
fprintf('[controller-fit] total params = %d | mean RMSE = %.5f\n', ...
    totalLen, A.meanRMSE);
fprintf('[controller-fit] per-muscle RMSE = %s\n', mat2str(rmse, 4));
end

%% ========================================================================
function delete_if_exists(f)
if isfile(f)
    delete(f);
end
end
