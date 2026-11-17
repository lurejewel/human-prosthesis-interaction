function C = compare_bootstrap_resampling(stoFile, nSeeds)
%COMPARE_BOOTSTRAP_RESAMPLING  Run-level vs i.i.d. resampling, head to head.
%
%   C = compare_bootstrap_resampling(stoFile, nSeeds)
%
%   The LASSO solver (`sgl_fit`) treats every row as an independent sample and
%   never uses time, so one may expect run-level and i.i.d. resampling to give
%   the same stability.  The difference is in what the resample REPRESENTS:
%   i.i.d. draws break the within-stride correlation and scatter duplicated
%   frames across all strides, which makes each replicate look like it saw a
%   richer, more independent sample -> stability is usually OVERSTATED.
%   This function measures that directly.
%
%   It re-runs the whole identification with resampleMode='iid_bootstrap'
%   (K = nSeeds) and compares it with the saved run_bootstrap diagnostics in
%   results/lasso_stability_diagnostics.mat.
%
%   Output C: .run, .iid (consensus structs), .table (per-group summary),
%             .settings
%
%   See also lasso_stability_identification, bootstrap_phase_runs.

if nargin < 1 || isempty(stoFile)
    here = fileparts(mfilename('fullpath'));
    stoFile = fullfile(here, 'moco_dataset', 'moco_lasso_dataset.sto');
end
if nargin < 2 || isempty(nSeeds)
    nSeeds = 20;
end
here = fileparts(mfilename('fullpath'));
repoRoot = fileparts(here);
resultsDir = fullfile(repoRoot, 'results');

diagFile = fullfile(resultsDir, 'lasso_stability_diagnostics.mat');
if ~isfile(diagFile)
    error('compare_bootstrap_resampling:NoRunBootstrap', ...
        'Run Adaptive_LASSO_MocoInverse_Stability first (%s missing).', diagFile);
end
Rrun = load(diagFile, 'R');
Rrun = Rrun.R;

S = load(fullfile(here, 'moco_dataset', 'moco_lasso_dataset.mat'), 'ds');
fitOptions = struct();
fitOptions.gamma = 1; fitOptions.epsilonOLS = 1e-8; fitOptions.nzTol = 1e-5;
fitOptions.lambda1 = 0; fitOptions.c = 0.01; fitOptions.rho = 50;
fitOptions.delaySteps = 1; fitOptions.excitationThreshold = 0.01;
fitOptions.robustInitIter = 2; fitOptions.robustIter = 2;
fitOptions.huberDeltaFactor = 1.345;
fitOptions.source = 'MocoInverse';
fitOptions.groupLabels = {'stance', 'swing'};
fitOptions.groupPhases = {[0 1], [2 3 4]};
fitOptions.lambdaPathConfig = struct('nLambda', 60, 'lambdaMinRatio', 1e-4, ...
    'tolSelect', 0.05, 'selectionMetric', 'post_refit_train_censored_rmse');
fitOptions.sampleWeight = S.ds.sampleWeight(:);
fitOptions.preprocessing = S.ds.meta;
fitOptions.verbose = false;

opt = struct();
opt.nSeeds = nSeeds;
opt.baseSeed = 2026;
opt.useParallel = true;
opt.saveResults = false;      % never overwrite the run-bootstrap artifacts
opt.verbose = true;
opt.fitOptions = fitOptions;
opt.resampleMode = 'iid_bootstrap';
opt.computeReferenceFit = false;   % the reference fit is mode-independent

fprintf('[compare] running %d replicates with resampleMode=iid_bootstrap ...\n', nSeeds);
Riid = lasso_stability_identification(stoFile, opt);

Cr = Rrun.consensus;
Ci = Riid.consensus;
groupLabels = fitOptions.groupLabels;

fprintf('\n===== RUN-LEVEL vs I.I.D. BOOTSTRAP (K=%d) =====\n', nSeeds);
fprintf('%-34s | %-17s | %-17s\n', 'metric', groupLabels{1}, groupLabels{2});
rows = {};
metricNames = {'mean pairwise Jaccard', 'stable-on (pi>=0.9)', ...
    'stable-off (pi<=0.1)', 'ambiguous', 'majority mask nnz', ...
    'modal mask fraction', 'mean coverage'};
valsRun = [Cr.jaccardMean; Cr.nStableOn; Cr.nStableOff; Cr.nAmbiguous; ...
    cellfun(@nnz, Cr.mask); Cr.modeFraction; Cr.meanCoverage];
valsIid = [Ci.jaccardMean; Ci.nStableOn; Ci.nStableOff; Ci.nAmbiguous; ...
    cellfun(@nnz, Ci.mask); Ci.modeFraction; Ci.meanCoverage];
fprintf('%-34s | %8s %8s | %8s %8s\n', '', 'run', 'iid', 'run', 'iid');
for i = 1:numel(metricNames)
    fprintf('%-34s | %8.3f %8.3f | %8.3f %8.3f\n', metricNames{i}, ...
        valsRun(i, 1), valsIid(i, 1), valsRun(i, 2), valsIid(i, 2));
end

C = struct();
C.run = Cr;
C.iid = Ci;
C.metricNames = metricNames;
C.valsRun = valsRun;
C.valsIid = valsIid;
C.settings = struct('nSeeds', nSeeds, 'baseSeed', 2026);
C.files = fullfile(resultsDir, 'lasso_bootstrap_mode_comparison.mat');
save(C.files, 'C', '-v7.3');
fprintf('\n[compare] saved %s\n', C.files);
end
