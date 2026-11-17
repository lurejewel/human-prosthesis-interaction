function ok = check_lasso_refactor_equivalence(mode, baselineFile)
%CHECK_LASSO_REFACTOR_EQUIVALENCE  Prove a LASSO refactor is numerically inert.
%
%   check_lasso_refactor_equivalence('record')            % before refactoring
%   ok = check_lasso_refactor_equivalence('check')        % after refactoring
%
%   The adaptive-LASSO fit was split into
%       load_lasso_design_matrices.m  (data loading, once)
%       fit_adaptive_lasso_core.m     (fitting, per resample)
%   with fit_adaptive_lasso_phase_model.m left as a thin wrapper.  Such a
%   refactor is only safe if it changes no number, so this script records a
%   baseline BEFORE the refactor and re-runs the identical call AFTER it.
%
%   Usage
%   -----
%   1) on the pre-refactor checkout / git stash:
%          check_lasso_refactor_equivalence('record')
%   2) after the refactor:
%          check_lasso_refactor_equivalence('check')
%      -> asserts mask / beta / bias / RMSE / support / selectedIdx are equal
%         to 1e-12 and errors otherwise.
%
%   The baseline and the data paths are configurable via the second argument
%   (default: results/baseline_fit_pre_refactor.mat).
%
%   NOTE: this needs MATLAB + OpenSim (it reads the real .sto), unlike
%   tests/test_lasso_stability_identification.m which is pure MATLAB.

if nargin < 1 || isempty(mode)
    mode = 'check';
end
thisDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(thisDir);
lassoDir = fullfile(repoRoot, 'Sparse Group LASSO Validation');
addpath(lassoDir);

if nargin < 2 || isempty(baselineFile)
    baselineFile = fullfile(repoRoot, 'results', 'baseline_fit_pre_refactor.mat');
end

stoFile = fullfile(lassoDir, 'moco_dataset', 'moco_lasso_dataset.sto');
matFile = fullfile(lassoDir, 'moco_dataset', 'moco_lasso_dataset.mat');
if ~isfile(stoFile)
    error('check_lasso_refactor_equivalence:NoData', ...
        'Dataset not found: %s (run moco_inverse_to_lasso_dataset first).', stoFile);
end

options = base_options(matFile);

switch lower(mode)
    case 'record'
        t0 = tic;
        [fit, lasso] = fit_adaptive_lasso_phase_model(stoFile, options);
        elapsed = toc(t0); %#ok<NASGU>
        [d, ~, ~] = fileparts(baselineFile);
        if ~isempty(d) && ~exist(d, 'dir')
            mkdir(d);
        end
        save(baselineFile, 'fit', 'lasso', 'elapsed', '-v7.3');
        fprintf('[record] baseline written to %s (%.1f s)\n', baselineFile, elapsed);
        ok = true;
        return;

    case 'check'
        if ~isfile(baselineFile)
            error('check_lasso_refactor_equivalence:NoBaseline', ...
                ['Baseline not found: %s\n' ...
                 'Run check_lasso_refactor_equivalence(''record'') on the ' ...
                 'pre-refactor code first.'], baselineFile);
        end
        B = load(baselineFile, 'fit', 'lasso');
        t0 = tic;
        [fit2, lasso2] = fit_adaptive_lasso_phase_model(stoFile, options);
        elapsed = toc(t0);

        fails = compare_fits(B.fit, B.lasso, fit2, lasso2);
        fprintf('\n===== REFACTOR EQUIVALENCE =====\n');
        fprintf('new run: %.1f s\n', elapsed);
        if isempty(fails)
            fprintf('RESULT: PASS - masks/beta/bias/RMSE/support/selectedIdx identical (tol 1e-12)\n');
            fprintf('masks: stance nnz=%d, swing nnz=%d, total params=%d\n', ...
                nnz(lasso2.groups(1).mask), nnz(lasso2.groups(2).mask), ...
                nnz(lasso2.groups(1).mask) + 9 + nnz(lasso2.groups(2).mask) + 9);
            ok = true;
        else
            fprintf('RESULT: FAIL\n  %s\n', fails{:});
            ok = false;
            error('check_lasso_refactor_equivalence:Failed', ...
                'Refactor changed numerics (%d difference(s)).', numel(fails));
        end

    otherwise
        error('check_lasso_refactor_equivalence:BadMode', ...
            'mode must be ''record'' or ''check''.');
end
end

%% ========================================================================
function options = base_options(matFile)
% Same identification settings as Adaptive_LASSO_MocoInverse(_Stability).m.
S = load(matFile, 'ds');
options = struct();
options.gamma            = 1;
options.epsilonOLS       = 1e-8;
options.nzTol            = 1e-5;
options.lambda1          = 0;
options.c                = 0.01;
options.rho              = 50;
options.delaySteps       = 1;
options.excitationThreshold = 0.01;
options.robustInitIter   = 2;
options.robustIter       = 2;
options.huberDeltaFactor = 1.345;
options.source           = 'MocoInverse';
options.groupLabels      = {'stance', 'swing'};
options.groupPhases      = {[0 1], [2 3 4]};
options.lambdaPathConfig = struct('nLambda', 60, 'lambdaMinRatio', 1e-4, ...
    'tolSelect', 0.05, 'selectionMetric', 'post_refit_train_censored_rmse');
options.sampleWeight     = S.ds.sampleWeight(:);
options.preprocessing    = S.ds.meta;
options.verbose          = false;
end

%% ========================================================================
function fails = compare_fits(fitA, lassoA, fitB, lassoB)
tol = 1e-12;
fails = {};

for g = 1:numel(lassoA.groups)
    if ~isequal(lassoA.groups(g).mask, lassoB.groups(g).mask)
        fails{end + 1} = sprintf('group %d: mask differs', g); %#ok<AGROW>
    end
    d = max(abs(lassoA.groups(g).beta(:) - lassoB.groups(g).beta(:)));
    if d > tol
        fails{end + 1} = sprintf('group %d: beta maxdiff %.3g', g, d); %#ok<AGROW>
    end
    d = max(abs(lassoA.groups(g).bias(:) - lassoB.groups(g).bias(:)));
    if d > tol
        fails{end + 1} = sprintf('group %d: bias maxdiff %.3g', g, d); %#ok<AGROW>
    end
end

d = max(abs(fitA.RMSE_refit_full(:) - fitB.RMSE_refit_full(:)));
if d > tol
    fails{end + 1} = sprintf('RMSE_refit_full maxdiff %.3g', d);
end

for ph = {'stance', 'swing'}
    d = max(abs(fitA.RMSE_refit.(ph{1})(:) - fitB.RMSE_refit.(ph{1})(:)));
    if d > tol
        fails{end + 1} = sprintf('RMSE_refit.%s maxdiff %.3g', ph{1}, d); %#ok<AGROW>
    end
    for m = 1:numel(fitA.S.(ph{1}))
        if ~isequal(fitA.S.(ph{1}){m}(:), fitB.S.(ph{1}){m}(:))
            fails{end + 1} = sprintf('S.%s{%d} differs', ph{1}, m); %#ok<AGROW>
        end
        iA = fitA.LambdaPath.(ph{1}){m}.selectedIdx;
        iB = fitB.LambdaPath.(ph{1}){m}.selectedIdx;
        if ~isequal(iA, iB)
            fails{end + 1} = sprintf('selectedIdx.%s{%d}: %d vs %d', ...
                ph{1}, m, iA, iB); %#ok<AGROW>
        end
    end
end
end
