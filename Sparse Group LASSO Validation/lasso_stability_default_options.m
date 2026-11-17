function opts = lasso_stability_default_options(options)
%LASSO_STABILITY_DEFAULT_OPTIONS  Defaults for multi-seed structure identification.
%
%   opts = lasso_stability_default_options(options)
%
%   Fields
%   ------
%   nSeeds                20      number of random seeds (= replicates)
%   baseSeed              2026    replicate k uses baseSeed + k
%   resampleMode          'run_bootstrap'  'run_bootstrap'|'run_subsample'|'iid_bootstrap'
%   subsampleFraction     0.8     fraction of runs kept by 'run_subsample'
%   maxBlockLen           0       split runs into blocks of at most this length
%   consensusRule         'majority'  'majority'|'mode'|'medoid'|'stability'
%   piThreshold           0.7     threshold of the 'stability' rule
%   tieBreak              'include'   Pi == 0.5 resolution for majority
%   stableOnPi            0.9     diagnostic threshold
%   stableOffPi           0.1     diagnostic threshold
%   clusterJaccard        0.9     mask-similarity cluster threshold
%   conformJaccard        0.8     a replicate supports the consensus when its
%                                 Jaccard to the consensus mask is >= this
%   minDominantFraction   0.6     reported only (strict "same framework" share)
%   minStabilityJaccard   0.35    WARN if the mean pairwise Jaccard between
%                                 replicate masks is below this -- the honest
%                                 "is the 0-1 structure reproducible" gate
%   maxRMSEIncrease       0.15    warn if the frozen-mask refit RMSE grows more
%                                  than this versus an UNFROZEN full-data fit
%                                  (replicate RMSE is in-bag, so it is not used
%                                  as the reference)
%   computeReferenceFit   true    run that unfrozen full-data fit once for the gate
%   useParallel           true
%   workerThreads         1       BLAS threads per parfor worker; these are
%                                 26x26 least-squares problems, so multithreaded
%                                 BLAS only causes oversubscription (16 workers
%                                 x N threads).  [] disables the setting.
%   nRetry                3       retries per replicate (seed + 1000*attempt)
%   saveResults           true
%   outFile               'results/lasso_controller_result_moco_consensus.mat'
%   diagFile              'results/lasso_stability_diagnostics.mat'
%   checkpointFile        'results/lasso_stability_checkpoint.mat'
%   resumeFromCheckpoint  false
%   verbose               true
%   fitOptions            struct()  merged into lasso_fit_default_options
%
%   NOTE: the fitting options inside fitOptions are deliberately shared by all
%   replicates -- the random seed changes only the resampled data, never the
%   identification procedure.

if nargin < 1 || isempty(options)
    options = struct();
end

defaults = struct();
defaults.nSeeds = 20;
defaults.baseSeed = 2026;
defaults.resampleMode = 'run_bootstrap';
defaults.subsampleFraction = 0.8;
defaults.maxBlockLen = 0;
defaults.consensusRule = 'majority';
defaults.piThreshold = 0.7;
defaults.tieBreak = 'include';
defaults.stableOnPi = 0.9;
defaults.stableOffPi = 0.1;
defaults.clusterJaccard = 0.9;
defaults.conformJaccard = 0.8;
defaults.minDominantFraction = 0.6;
defaults.minStabilityJaccard = 0.35;
defaults.maxRMSEIncrease = 0.15;
defaults.computeReferenceFit = true;   % unfrozen full-data fit for a fair RMSE gate
defaults.useParallel = true;
defaults.workerThreads = 1;    % BLAS threads per parfor worker (tiny matrices)
defaults.nRetry = 3;
defaults.saveResults = true;
defaults.outFile = fullfile('results', 'lasso_controller_result_moco_consensus.mat');
defaults.diagFile = fullfile('results', 'lasso_stability_diagnostics.mat');
defaults.checkpointFile = fullfile('results', 'lasso_stability_checkpoint.mat');
defaults.resumeFromCheckpoint = false;
defaults.verbose = true;
defaults.fitOptions = struct();

fn = fieldnames(defaults);
opts = options;
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}) || isempty(opts.(fn{i}))
        opts.(fn{i}) = defaults.(fn{i});
    end
end

% ---- validation ---------------------------------------------------------
if ~(isscalar(opts.nSeeds) && opts.nSeeds >= 2 && floor(opts.nSeeds) == opts.nSeeds)
    error('lasso_stability_default_options:BadNSeeds', ...
        'nSeeds must be an integer >= 2 (got %s).', num2str(opts.nSeeds));
end
if ~(isscalar(opts.baseSeed) && floor(opts.baseSeed) == opts.baseSeed)
    error('lasso_stability_default_options:BadBaseSeed', ...
        'baseSeed must be an integer scalar.');
end
if ~any(strcmpi(opts.resampleMode, {'run_bootstrap', 'run_subsample', 'iid_bootstrap'}))
    error('lasso_stability_default_options:BadMode', ...
        'resampleMode must be run_bootstrap, run_subsample or iid_bootstrap.');
end
if ~any(strcmpi(opts.consensusRule, {'majority', 'mode', 'medoid', 'stability'}))
    error('lasso_stability_default_options:BadRule', ...
        'consensusRule must be majority, mode, medoid or stability.');
end
if ~(isscalar(opts.minDominantFraction) && opts.minDominantFraction > 0 && ...
        opts.minDominantFraction <= 1)
    error('lasso_stability_default_options:BadFraction', ...
        'minDominantFraction must be in (0, 1].');
end
if ~isstruct(opts.fitOptions)
    error('lasso_stability_default_options:BadFitOptions', ...
        'fitOptions must be a struct.');
end
end
