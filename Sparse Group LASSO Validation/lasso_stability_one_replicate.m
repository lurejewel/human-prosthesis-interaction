function r = lasso_stability_one_replicate(D, k, opts)
%LASSO_STABILITY_ONE_REPLICATE  One bootstrap replicate of the LASSO structure.
%
%   r = lasso_stability_one_replicate(D, k, opts)
%
%   Draws a resample of the gait-phase runs for replicate k (seed = baseSeed+k),
%   re-runs the unchanged adaptive-LASSO pipeline on it, and returns a
%   LIGHTWEIGHT result: only the 0-1 masks, per-muscle RMSE and a few scalars.
%   The full fit objects are intentionally not kept -- K of them would be
%   hundreds of megabytes and nothing downstream needs them.
%
%   D      : design system from load_lasso_design_matrices
%   k      : replicate index (1..nSeeds)
%   opts   : resolved struct from lasso_stability_default_options, with
%            opts.fitOptions already merged through lasso_fit_default_options
%
%   r fields: k, seed, ok, error, attempts, masks{1 x nGroups},
%             rmse{1 x nGroups} (1 x nMus), supportSize{1 x nGroups},
%             nSamples, nStance, nSwing, tElapsed
%
%   A failing replicate is retried with seed + 1000*attempt; if all attempts
%   fail, r.ok stays false and the replicate is excluded from the consensus
%   for every group (see consensus_reflex_structure's opts.valid).

nGroups = numel(opts.fitOptions.groupPhases);

% Inside a parfor worker, pin BLAS to a few threads.  These are 26x26
% least-squares problems; with 16 workers the default per-worker threading
% oversubscribes the CPU and makes the ensemble much slower.
if isfield(opts, 'workerThreads') && ~isempty(opts.workerThreads)
    try
        if ~isempty(getCurrentTask())
            maxNumCompThreads(opts.workerThreads);
        end
    catch
        % getCurrentTask/maxNumCompThreads unavailable -> leave as is
    end
end

r = struct();
r.k = k;
r.seed = NaN;
r.ok = false;
r.error = '';
r.attempts = 0;
r.masks = cell(1, nGroups);
r.rmse = cell(1, nGroups);
r.supportSize = cell(1, nGroups);
r.nSamples = NaN;
r.nStance = NaN;
r.nSwing = NaN;
r.tElapsed = NaN;

t0 = tic;
seed0 = opts.baseSeed + k;

for attempt = 1:opts.nRetry
    seed = seed0 + 1000 * (attempt - 1);
    try
        rs = RandStream('Threefry', 'Seed', seed);
        S = build_resampled_system(D, opts, rs);

        fitOpts = opts.fitOptions;
        fitOpts.verbose = false;
        fitOpts.storeLambdaPath = false;
        fitOpts.suppressWarnings = true;

        [fitK, lassoK] = fit_adaptive_lasso_core(S, fitOpts);

        for g = 1:nGroups
            r.masks{g} = lassoK.groups(g).mask;
            r.rmse{g} = fitK.RMSE_refit_full(:).';
            r.supportSize{g} = nnz(lassoK.groups(g).mask);
        end
        r.seed = seed;
        r.ok = true;
        r.attempts = attempt;
        r.error = '';
        r.nSamples = size(S.X, 1);
        r.nStance = numel(S.phaseIdxMap.stance);
        r.nSwing = numel(S.phaseIdxMap.swing);
        r.tElapsed = toc(t0);
        return;
    catch ME
        r.error = ME.message;
    end
end

r.attempts = opts.nRetry;
r.tElapsed = toc(t0);
end

%% ========================================================================
% Local: resample the two phase groups and build a local design system.
% ========================================================================
% The resampled rows are stacked as [stance block; swing block]; the local
% phase index map is therefore trivial (1:nStance and nStance+1:end).  The
% standardisation constants muX/sX are kept from the FULL data on purpose:
% fixing them removes one nuisance source of variation and keeps coefficients
% comparable across replicates.
function S = build_resampled_system(D, opts, rs)
bsOpts = struct();
bsOpts.randStream = rs;
bsOpts.maxBlockLen = opts.maxBlockLen;
bsOpts.subsampleFraction = opts.subsampleFraction;

idxStance = bootstrap_phase_runs(D.phaseIdxMap.stance, opts.resampleMode, bsOpts);
idxSwing  = bootstrap_phase_runs(D.phaseIdxMap.swing,  opts.resampleMode, bsOpts);

allIdx = [idxStance; idxSwing];

S = struct();
S.X = D.X(allIdx, :);
S.X_std = D.X_std(allIdx, :);
S.muX = D.muX;
S.sX = D.sX;
S.Y = D.Y(allIdx, :);
S.w = D.w(allIdx);
nS = numel(idxStance);
nW = numel(idxSwing);
S.phaseIdxMap = struct('stance', (1:nS).', 'swing', ((nS + 1):(nS + nW)).');
S.muscleNames = D.muscleNames;
S.featureNames = D.featureNames;
if isfield(D, 'stoFile')
    S.stoFile = D.stoFile;
end
end
