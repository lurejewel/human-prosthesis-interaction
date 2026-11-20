function opts = lasso_fit_default_options(options)
%LASSO_FIT_DEFAULT_OPTIONS  Shared option defaults for LASSO structure fitting.
%
%   opts = lasso_fit_default_options(options) fills every missing or empty
%   field of `options` with the defaults used by fit_adaptive_lasso_core.m.
%   A field explicitly set to [] is treated as missing (legacy behaviour of
%   the old local apply_defaults in fit_adaptive_lasso_phase_model.m).
%
%   Fields (defaults)
%   ----------------
%   gamma 1, epsilonOLS 1e-8, nzTol 1e-5, lambda1 0
%   lambdaPathConfig.nLambda 60
%   lambdaPathConfig.lambdaMinRatio 1e-4
%   lambdaPathConfig.tolSelect 0.05  scalar, or 1 x nGroups (e.g. [0.05 0.20]
%   lambdaPathConfig.selectionMode 'tolSelect'  'tolSelect' | 'knee'
%   lambdaPathConfig.minSupport 1     knee mode: require support >= this
%   lambdaPathConfig.selectionMetric 'post_refit_train_censored_rmse'
%   c 0.01, rho 50, delaySteps 1, excitationThreshold 0.01
%   sampleWeight []            (all ones; length n or n+delaySteps)
%   robustInitIter 2, robustIter 2, huberDeltaFactor 1.345
%   verbose true
%   groupLabels {'stance','swing'}, groupPhases {[0 1],[2 3 4]}
%   source 'MocoInverse', preprocessing struct()
%
%   Structure-identification extensions
%   -----------------------------------
%   storeLambdaPath   true    keep the full per-lambda beta paths inside fit
%   maskOverride      {}      {1 x nGroups} logical masks -> refit-only mode
%   suppressWarnings  false   silence the rank-deficiency warning
%
%   The nested lambdaPathConfig is merged field-by-field: user overrides win,
%   missing sub-fields fall back to the defaults.  If the caller supplies
%   nLambda (or lambdaMinRatio) without an explicit alphaGrid, the grid is
%   regenerated from the new values, exactly as the legacy implementation did.

if nargin < 1 || isempty(options)
    options = struct();
end

defaults = struct();
defaults.gamma = 1;
defaults.epsilonOLS = 1e-8;
defaults.nzTol = 1e-5;
defaults.lambda1 = 0;
defaults.lambdaPathConfig = struct( ...
    'nLambda', 60, ...
    'lambdaMinRatio', 1e-4, ...
    'tolSelect', 0.05, ...
    'selectionMode', 'tolSelect', ...   % 'tolSelect' | 'knee'
    'minSupport', 1, ...                % knee mode: smallest admissible support
    'selectionMetric', 'post_refit_train_censored_rmse');
defaults.c = 0.01;
defaults.rho = 50;
defaults.delaySteps = 1;
defaults.excitationThreshold = 0.01;
defaults.sampleWeight = [];
defaults.robustInitIter = 2;
defaults.robustIter = 2;
defaults.huberDeltaFactor = 1.345;
defaults.verbose = true;
defaults.groupLabels = {'stance', 'swing'};
defaults.groupPhases = {[0 1], [2 3 4]};
defaults.source = 'MocoInverse';
defaults.preprocessing = struct();
defaults.storeLambdaPath = true;
defaults.maskOverride = {};
defaults.suppressWarnings = false;

fn = fieldnames(defaults);
opts = options;
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}) || isempty(opts.(fn{i}))
        opts.(fn{i}) = defaults.(fn{i});
    end
end

% Nested lambdaPathConfig: fill missing sub-fields, keep user overrides.
lpc = defaults.lambdaPathConfig;
userHasAlphaGrid = false;
if isfield(options, 'lambdaPathConfig') && isstruct(options.lambdaPathConfig)
    u = options.lambdaPathConfig;
    if isfield(u, 'alphaGrid') && ~isempty(u.alphaGrid)
        userHasAlphaGrid = true;
    end
    fn2 = fieldnames(u);
    for i = 1:numel(fn2)
        lpc.(fn2{i}) = u.(fn2{i});
    end
end
if ~userHasAlphaGrid
    lpc.alphaGrid = logspace(0, log10(lpc.lambdaMinRatio), lpc.nLambda).';
end
opts.lambdaPathConfig = lpc;

% maskOverride must be a (possibly empty) cell array of logical masks.
if ~iscell(opts.maskOverride)
    error('lasso_fit_default_options:BadMaskOverride', ...
        'options.maskOverride must be a cell array of logical masks.');
end
end
