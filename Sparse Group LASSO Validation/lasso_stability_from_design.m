function R = lasso_stability_from_design(D, options)
%LASSO_STABILITY_FROM_DESIGN  Multi-seed consensus structure identification.
%
%   R = lasso_stability_from_design(D, options)
%
%   Runs K replicates of the (unchanged) adaptive-LASSO identification on
%   run-level bootstrap resamples of the same experimental data, converts each
%   replicate into a 0-1 mask per phase group, and turns the ensemble into one
%   consensus controller:
%
%     1. draw K resamples (seed = baseSeed + k);
%     2. fit every resample with the SAME identification procedure;
%     3. collect the 26x9 masks and their selection frequencies Pi;
%     4. consensus mask (default: per-entry majority vote, Pi >= 0.5);
%     5. refit the coefficients on the FULL data with the consensus mask frozen
%        -> the parameters CMA-ES starts from.
%
%   Nothing here depends on OpenSim, so it can be exercised with a synthetic D
%   (see tests/test_lasso_stability_identification.m).  Use
%   lasso_stability_identification.m to go from a .sto file to the controller.
%
%   D       : design system from load_lasso_design_matrices (FULL data)
%   options : see lasso_stability_default_options
%
%   Output R
%   --------
%   .lasso        consensus controller (grouped v2 + .identification/.stability)
%   .fitFinal     full-data refit with the frozen mask
%   .consensus    everything from consensus_reflex_structure
%   .replicates   1 x K lightweight replicate results
%   .valid        K x nGroups logical
%   .options      resolved options
%   .rmseDelta    relative RMSE change of the frozen-mask refit vs replicates
%   .warnings     cellstr of issued warnings
%   .timing       .totalSeconds / .replicateSeconds
%
%   See also lasso_stability_identification, lasso_stability_one_replicate,
%            consensus_reflex_structure.

if nargin < 2 || isempty(options)
    options = struct();
end
opts = lasso_stability_default_options(options);
fitOpts = lasso_fit_default_options(opts.fitOptions);
opts.fitOptions = fitOpts;

K = opts.nSeeds;
nGroups = numel(fitOpts.groupPhases);
muscleNames = D.muscleNames;
t0 = tic;

if opts.verbose
    fprintf('\n========== LASSO stability identification ==========\n');
    fprintf('  seeds        : %d (baseSeed=%d, %s)\n', K, opts.baseSeed, opts.resampleMode);
    fprintf('  consensus    : %s (piThreshold=%.2f, tieBreak=%s)\n', ...
        opts.consensusRule, opts.piThreshold, opts.tieBreak);
    fprintf('  data         : %d samples (stance %d / swing %d), %d muscles\n', ...
        D.nSamples, numel(D.phaseIdxMap.stance), numel(D.phaseIdxMap.swing), numel(muscleNames));
end

%% ===== EXISTING REPLICATES (RESUME) =====================================
rep = empty_reps(K, nGroups);
completed = false(K, 1);
if opts.resumeFromCheckpoint && isfile(opts.checkpointFile)
    ck = load(opts.checkpointFile);
    okCompatible = isfield(ck, 'rep') && numel(ck.rep) == K && ...
        isfield(ck, 'nSeeds') && isequal(ck.nSeeds, K) && ...
        isfield(ck, 'baseSeed') && isequal(ck.baseSeed, opts.baseSeed);
    if okCompatible
        rep = ck.rep(:).';
        if isfield(ck, 'completed') && numel(ck.completed) == K
            completed = logical(ck.completed(:));
        end
        % A replicate that was attempted but failed must be retried on resume.
        for kk = 1:K
            completed(kk) = completed(kk) && isfield(rep(kk), 'ok') && rep(kk).ok;
        end
        if opts.verbose
            fprintf('[stability] resumed %d/%d replicates from %s\n', ...
                nnz(completed), K, opts.checkpointFile);
        end
    else
        warning('lasso_stability:CheckpointIncompatible', ...
            'Checkpoint exists but does not match nSeeds/baseSeed; ignoring it.');
    end
end

todo = find(~completed);
todo = todo(:).';

%% ===== REPLICATE LOOP ===================================================
if ~isempty(todo)
    res = cell(numel(todo), 1);

    usePar = opts.useParallel && numel(todo) > 1;
    if usePar
        pool = gcp('nocreate');
        if isempty(pool)
            try
                nWorkers = min(feature('numcores'), numel(todo));
                parpool('local', nWorkers);
            catch ME
                warning('lasso_stability:NoPool', ...
                    'Could not start a parallel pool (%s); falling back to serial.', ME.message);
                usePar = false;
            end
        end
    end

    if usePar
        if opts.verbose
            fprintf('[stability] running %d replicates on a parallel pool ...\n', numel(todo));
        end
        parfor ii = 1:numel(todo)
            res{ii} = lasso_stability_one_replicate(D, todo(ii), opts);
        end
    else
        for ii = 1:numel(todo)
            res{ii} = lasso_stability_one_replicate(D, todo(ii), opts);
            if opts.verbose
                rr = res{ii};
                fprintf('[stability] %2d/%2d (seed %d) %s | stance nnz=%d swing nnz=%d | %.1f s\n', ...
                    ii, numel(todo), rr.seed, ternary(rr.ok, 'ok  ', 'FAIL'), ...
                    local_nnz(rr, 1), local_nnz(rr, 2), rr.tElapsed);
            end
        end
    end

    for ii = 1:numel(todo)
        rep(todo(ii)) = res{ii};
    end
    completed = true(K, 1);

    if opts.saveResults
        ensure_dir(opts.checkpointFile);
        nSeeds = K; baseSeed = opts.baseSeed; timestamp = char(datetime);
        save(opts.checkpointFile, 'rep', 'completed', 'nSeeds', 'baseSeed', 'timestamp');
    end
end

%% ===== ASSEMBLE MASKS / VALIDITY ========================================
masks = cell(K, nGroups);
valid = false(K, nGroups);
failed = {};
for k = 1:K
    if isfield(rep(k), 'ok') && rep(k).ok
        for g = 1:nGroups
            masks{k, g} = rep(k).masks{g};
            valid(k, g) = true;
        end
    else
        failed{end + 1} = sprintf('#%d: %s', k, rep(k).error); %#ok<AGROW>
    end
end
if ~any(valid(:))
    error('lasso_stability:AllReplicatesFailed', ...
        'Every replicate failed. First error: %s', strjoin(failed, ' | '));
end
if opts.verbose && ~isempty(failed)
    fprintf('[stability] %d failed replicate(s):\n', numel(failed));
    fprintf('    %s\n', failed{:});
end

%% ===== CONSENSUS ========================================================
cOpts = struct();
cOpts.rule = opts.consensusRule;
cOpts.piThreshold = opts.piThreshold;
cOpts.tieBreak = opts.tieBreak;
cOpts.stableOnPi = opts.stableOnPi;
cOpts.stableOffPi = opts.stableOffPi;
cOpts.clusterJaccard = opts.clusterJaccard;
cOpts.conformJaccard = opts.conformJaccard;
cOpts.valid = valid;
C = consensus_reflex_structure(masks, cOpts);

%% ===== WARNINGS / DIAGNOSTIC MESSAGES ===================================
warnList = {};
for g = 1:nGroups
    gLabel = fitOpts.groupLabels{g};
    % Stability gate: how similar are the replicate masks to EACH OTHER?
    % This is the honest "is the 0-1 structure reproducible" measure; the
    % consensus itself is a superset, so its Jaccard to a single seed is low by
    % construction and must not be used as a gate.
    if C.jaccardMean(g) < opts.minStabilityJaccard
        msg = sprintf(['group %s: mean pairwise Jaccard between replicate masks ' ...
            'is %.2f (< minStabilityJaccard %.2f). Only %d of %d entries are ' ...
            'stable-on (pi>=%.1f) and %d are ambiguous; the 0-1 structure is ' ...
            'weakly reproducible, so treat individual entries with care.'], ...
            gLabel, C.jaccardMean(g), opts.minStabilityJaccard, ...
            C.nStableOn(g), C.nElem, opts.stableOnPi, C.nAmbiguous(g));
        warnList{end + 1} = msg; %#ok<AGROW>
        warning('lasso_stability:WeakStability', '%s', msg);
    end
    if nnz(C.mask{g}) == 0
        error('lasso_stability:EmptyConsensus', ...
            'Consensus mask for group %s is empty.', gLabel);
    end
end

%% ===== FINAL REFIT ON FULL DATA WITH THE FROZEN MASK ====================
fitFinalOpts = fitOpts;
fitFinalOpts.maskOverride = C.mask;
fitFinalOpts.verbose = opts.verbose;
fitFinalOpts.storeLambdaPath = false;
[fitFinal, lassoFinal] = fit_adaptive_lasso_core(D, fitFinalOpts);

% ---- diagnostic: frozen refit vs the replicate distribution -------------
% Replicate RMSE is IN-BAG on a bootstrap resample, hence optimistically
% biased; this number is informational only and must not gate anything.
rmseVsReplicates = nan(1, nGroups);
for g = 1:nGroups
    vk = find(valid(:, g));
    Rmat = nan(numel(vk), numel(muscleNames));
    for i = 1:numel(vk)
        Rmat(i, :) = rep(vk(i)).rmse{g};
    end
    medRep = median(Rmat, 1);
    rmseVsReplicates(g) = mean( ...
        (fitFinal.RMSE_refit_full(:).' - medRep) ./ max(medRep, eps));
end

% ---- apples-to-apples gate: frozen refit vs an UNfrozen full-data fit ---
% Both are evaluated on the full data, so the only difference is the
% structure: this is exactly "what does freezing the consensus mask cost?".
rmseDelta = nan(1, numel(muscleNames));
fitReference = [];
lassoReference = [];
referenceTotalLen = NaN;
if opts.computeReferenceFit
    refOpts = fitOpts;
    refOpts.verbose = false;
    refOpts.storeLambdaPath = false;
    [fitReference, lassoReference] = fit_adaptive_lasso_core(D, refOpts);
    rmseDelta = (fitFinal.RMSE_refit_full(:).' - fitReference.RMSE_refit_full(:).') ...
        ./ max(fitReference.RMSE_refit_full(:).', eps);
    referenceTotalLen = 0;
    for g = 1:nGroups
        referenceTotalLen = referenceTotalLen + nnz(lassoReference.groups(g).mask) + 9;
    end
    if mean(rmseDelta) > opts.maxRMSEIncrease
        msg = sprintf(['frozen-mask refit RMSE is %.1f%% above the full-data fit ' ...
            '(threshold %.1f%%); the consensus structure may be too sparse or ' ...
            'missing a feature that the full-data selection needs.'], ...
            100 * mean(rmseDelta), 100 * opts.maxRMSEIncrease);
        warnList{end + 1} = msg; %#ok<AGROW>
        warning('lasso_stability:RMSEIncrease', '%s', msg);
    end
end

%% ===== METADATA ON THE CONTROLLER =======================================
maskSignature = zeros(1, nGroups);
totalLen = 0;
for g = 1:nGroups
    v = double(C.mask{g}(:));
    maskSignature(g) = mod(sum(v .* (1:numel(v)).'), 4294967296);
    totalLen = totalLen + nnz(C.mask{g}) + 9;
end

lassoFinal = attach_identification(lassoFinal, opts, C, maskSignature, ...
    totalLen, rmseDelta, rmseVsReplicates, rep, valid, D, fitOpts, warnList, ...
    referenceTotalLen);

%% ===== PACKAGE RESULT ===================================================
R = struct();
R.lasso = lassoFinal;
R.fitFinal = fitFinal;
R.consensus = C;
R.replicates = rep;
R.valid = valid;
R.options = opts;
R.rmseDelta = rmseDelta;
R.rmseVsReplicates = rmseVsReplicates;
R.fitReference = fitReference;
R.lassoReference = lassoReference;
R.referenceTotalLen = referenceTotalLen;
R.warnings = warnList;
R.maskSignature = maskSignature;
R.totalLen = totalLen;
R.timing = struct();
R.timing.totalSeconds = toc(t0);
R.timing.replicateSeconds = [rep.tElapsed];

if opts.verbose
    fprintf('\n===== CONSENSUS STRUCTURE =====\n');
    for g = 1:nGroups
        fprintf(['  %-6s : nnz=%3d / %d | modal %d/%d (%.2f) | coverage %.2f' ...
            ' | stable-on %d, ambiguous %d | mean Jaccard %.3f\n'], ...
            fitOpts.groupLabels{g}, nnz(C.mask{g}), C.nElem, C.modeCount(g), ...
            C.nValid(g), C.modeFraction(g), C.meanCoverage(g), ...
            C.nStableOn(g), C.nAmbiguous(g), C.jaccardMean(g));
    end
    fprintf('  total CMA-ES parameters = %d (frozen 0-1 structure)\n', totalLen);
    if ~isnan(referenceTotalLen)
        fprintf('  full-data single fit     = %d parameters (unfrozen reference)\n', ...
            referenceTotalLen);
    end
    fprintf('  RMSE vs full-data fit (mean over muscles) = %+.1f%%\n', ...
        100 * mean(rmseDelta));
    fprintf('  RMSE vs replicate median (in-bag, informational) = %+.1f%%\n', ...
        100 * mean(rmseVsReplicates));
    fprintf('  total time = %.1f s\n', R.timing.totalSeconds);
    fprintf('================================\n');
end

%% ===== SAVE =============================================================
if opts.saveResults
    ensure_dir(opts.outFile);
    ensure_dir(opts.diagFile);
    lasso = lassoFinal;
    save(opts.outFile, 'lasso');
    save(opts.diagFile, 'R', '-v7.3');
    if opts.verbose
        fprintf('[stability] controller  : %s\n', opts.outFile);
        fprintf('[stability] diagnostics : %s\n', opts.diagFile);
    end
end
end

%% ========================================================================
% Local: attach identification / stability metadata to the controller struct.
% ========================================================================
function lasso = attach_identification(lasso, opts, C, maskSignature, ...
    totalLen, rmseDelta, rmseVsReplicates, rep, valid, D, fitOpts, warnList, ...
    referenceTotalLen)
lasso.identification = struct();
lasso.identification.method = 'run_level_bootstrap_stability_selection';
lasso.identification.nSeeds = opts.nSeeds;
lasso.identification.baseSeed = opts.baseSeed;
lasso.identification.resampleMode = opts.resampleMode;
lasso.identification.subsampleFraction = opts.subsampleFraction;
lasso.identification.maxBlockLen = opts.maxBlockLen;
lasso.identification.consensusRule = opts.consensusRule;
lasso.identification.piThreshold = opts.piThreshold;
lasso.identification.tieBreak = opts.tieBreak;
lasso.identification.stableOnPi = opts.stableOnPi;
lasso.identification.stableOffPi = opts.stableOffPi;
lasso.identification.clusterJaccard = opts.clusterJaccard;
lasso.identification.minStabilityJaccard = opts.minStabilityJaccard;
lasso.identification.minDominantFraction = opts.minDominantFraction;
lasso.identification.nValid = C.nValid;
lasso.identification.modeCount = C.modeCount;
lasso.identification.modeFraction = C.modeFraction;
lasso.identification.conformJaccard = opts.conformJaccard;
lasso.identification.nConforming = C.nConforming;
lasso.identification.conformingFraction = C.conformingFraction;
lasso.identification.meanCoverage = C.meanCoverage;
lasso.identification.minStabilityJaccard = opts.minStabilityJaccard;
lasso.identification.dominantClusterFraction = C.dominantClusterFraction;
lasso.identification.jaccardMean = C.jaccardMean;
lasso.identification.nStableOn = C.nStableOn;
lasso.identification.nStableOff = C.nStableOff;
lasso.identification.nAmbiguous = C.nAmbiguous;
lasso.identification.maskSignature = maskSignature;
lasso.identification.totalParameters = totalLen;
lasso.identification.referenceTotalParameters = referenceTotalLen;
lasso.identification.rmseDelta = rmseDelta;            % vs unfrozen full-data fit
lasso.identification.rmseVsReplicates = rmseVsReplicates;  % in-bag, informational
lasso.identification.seeds = [rep.seed];
lasso.identification.validMask = valid;
lasso.identification.warnings = warnList;
lasso.identification.groupLabels = fitOpts.groupLabels;
lasso.identification.groupPhases = fitOpts.groupPhases;
lasso.identification.lambdaSelectionMode = fitOpts.lambdaPathConfig.selectionMode;
lasso.identification.minSupport = fitOpts.lambdaPathConfig.minSupport;
lasso.identification.tolSelect = fitOpts.lambdaPathConfig.tolSelect;
lasso.identification.featureNames = D.featureNames;
lasso.identification.created = char(datetime);

% Compact stability payload (the full jaccard matrices stay in the diagnostics
% file; the controller only carries what a consumer may want to inspect).
lasso.stability = struct();
lasso.stability.Pi = C.Pi;
lasso.stability.mask = C.mask;
lasso.stability.modeMask = C.modeMask;
lasso.stability.medoidMask = C.medoidMask;
lasso.stability.clusterSizes = C.clusterSizes;
lasso.stability.supportSizes = C.supportSizes;
end

%% ========================================================================
% Local helpers.
% ========================================================================
function rep = empty_reps(K, nGroups)
rep = repmat(struct('k', 0, 'seed', NaN, 'ok', false, 'error', '', ...
    'attempts', 0, 'masks', {cell(1, nGroups)}, 'rmse', {cell(1, nGroups)}, ...
    'supportSize', {cell(1, nGroups)}, 'nSamples', NaN, 'nStance', NaN, ...
    'nSwing', NaN, 'tElapsed', NaN), 1, K);
for k = 1:K
    rep(k).k = k;
end
end

function n = local_nnz(r, g)
if isfield(r, 'supportSize') && numel(r.supportSize) >= g && ~isempty(r.supportSize{g})
    n = r.supportSize{g};
else
    n = NaN;
end
end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end

function ensure_dir(f)
d = fileparts(f);
if ~isempty(d) && ~exist(d, 'dir')
    mkdir(d);
end
end
