function [betaOrig, intercept, betaStd, stats] = sgl_fit(A, y, lambda1, lambda2, opts)
%SGL_FIT Censored sparse-group-lasso regression via FISTA + backtracking.
%   [betaOrig, intercept, betaStd, stats] = sgl_fit(A, y, lambda1, lambda2, opts)
%   solves (Option 2 surrogate):
%     min_{beta,b}
%       0.5*|| y_U - (A_U*beta + b) ||_2^2
%     + (rho/2)*|| (A_C*beta + b - c)_+ ||_2^2
%     + lambda1*sum_g ||beta_g||_2 + lambda2*||beta||_1
%   where U = {i: y_i > c + tolCensor}, C = {i: y_i <= c + tolCensor},
%   and (t)_+ = max(t,0) elementwise.
%
%   IMPORTANT NOTES:
%   - y is NOT centered; threshold c is on ORIGINAL y scale and NEVER transformed.
%   - A is used as provided. If you want centering/scaling, do it before calling.
%   - c is a censoring floor on original scale: y_i = max(z_i, c) where z_i = A_i*beta_orig + b_orig
%   - Warm starts opts.beta0, opts.b0 are in the same feature scale as A.
%   - Outputs betaOrig, intercept are in the same feature scale as A.
%   - betaStd equals betaOrig when no external standardization is applied.
% 
%   Group structure is inferred from p = size(A,2), where p = 2*#mus + 8:
%     group1: 1:#mus
%     group2: #mus+1:2*#mus
%     group3: 2*#mus+1:2*#mus+8
%
%   Inputs:
%     A       - n x p predictor matrix (already masked feature matrix if needed upstream)
%     y       - n x 1 response vector
%     lambda1 - group-l2 penalty weight
%     lambda2 - l1 penalty weight
%     opts    - optional struct:
%               c             (default 0.01)
%               rho           (default 100 if empty; larger = harder censor constraint)
%               maxIter       (default 2000)
%               tol           (default 1e-6)
%               verbose       (default true)
%               backtrackBeta (default 0.5)
%               L0            (default 1.0)
%               beta0         (default []) warm start in ORIGINAL space
%               b0            (default 0) warm start intercept in ORIGINAL space
%               tolCensor     (default 1e-12)
%
%   Outputs:
%     betaOrig  - p x 1 coefficients on original A scale
%     intercept - optimized intercept bOrig on original y scale
%     betaStd   - p x 1 coefficients in internal optimization space
%     stats     - optimization traces and metadata

if nargin < 5
    opts = struct();
end
opts = apply_default_opts(opts);

y = y(:);

n = size(A, 1);
p = size(A, 2);
if numel(y) ~= n
    error('sgl_fit:BadResponseSize', 'y must have %d rows to match A.', n);
end

groups = group_indices(p);

if opts.backtrackBeta <= 0 || opts.backtrackBeta >= 1
    error('sgl_fit:BadBacktrackBeta', 'opts.backtrackBeta must be in (0, 1).');
end
if opts.L0 <= 0
    error('sgl_fit:BadL0', 'opts.L0 must be positive.');
end
%% ========== FEATURE PREPROCESSING ==========
% No internal centering/scaling. Use A as provided.
AStd = A;

c = opts.c;
rho = opts.rho;
tolCensor = opts.tolCensor;

%% ========== CENSOR SET PARTITION ==========
% Split data into uncensored (U) above threshold and censored (C) at/below threshold.
% This partition determines which objective term applies to each observation.
Uidx = find(y > c + tolCensor);        % uncensored indices: y_i strictly above c
Cidx = find(y <= c + tolCensor);       % censored indices: y_i at or below c (treated as y=c)

AU = AStd(Uidx, :);                    % |U| x p, rows of AStd for uncensored
yU = y(Uidx);                          % |U| x 1, responses for uncensored (unchanged scale)
AC = AStd(Cidx, :);                    % |C| x p, rows of AStd for censored
% Note: yC is not extracted; we use the fact that censored obs are modeled as y=c in the objective.

%% ========== INITIALIZE OPTIMIZATION VARIABLES FROM WARM STARTS ==========
% Warm starts are provided in the same feature scale as A.
betaOrig0 = zeros(p, 1);              % default: zero initialization
if ~isempty(opts.beta0)
    if ~isvector(opts.beta0) || numel(opts.beta0) ~= p
        error('sgl_fit:BadWarmStart', 'opts.beta0 must be a %d-by-1 vector.', p);
    end
    betaOrig0 = opts.beta0(:);
end
if ~(isscalar(opts.b0) && isnumeric(opts.b0) && isfinite(opts.b0))
    error('sgl_fit:BadWarmStartIntercept', 'opts.b0 must be a finite numeric scalar.');
end

bOrig0 = opts.b0;                              % intercept in original space

% Internal variables match original feature scale.
beta = betaOrig0;                               % p x 1
b = bOrig0;                                     % scalar

%% ========== FISTA INITIALIZATION ==========
% FISTA maintains current point (beta, b) and extrapolated point (yBeta, yb) for acceleration.
yBeta = beta;
yb = b;
tFista = 1;                                    % FISTA acceleration parameter (starts at 1)
step = 1 / max(opts.L0, eps);                  % initial step size (inverse of Lipschitz estimate)

% Record objective at generation 0 (initial warm start point)
smooth0 = compute_f_smooth(AU, yU, AC, c, rho, beta, b);
obj0 = smooth0 + lambda1 * group_l2_sum(beta, groups) + lambda2 * norm(beta, 1);

objHist = nan(opts.maxIter + 1, 1);
relChangeHist = nan(opts.maxIter, 1);
relObjHist = nan(opts.maxIter, 1);
stepHist = nan(opts.maxIter, 1);
smoothHist = nan(opts.maxIter + 1, 1);
objHist(1) = obj0;
smoothHist(1) = smooth0;
converged = false;
exitType = 'maxIter';
iters = 0;

%% ========== MAIN FISTA LOOP ==========
% FISTA (Fast Iterative Soft-Thresholding Algorithm) with backtracking line search.
% At each iteration:
%   1. Compute gradient at extrapolated point (yBeta, yb)
%   2. Backtrack to find a step size t satisfying descent condition
%   3. Apply proximal operator for SGL regularization at scaled step
%   4. Update extrapolation for momentum acceleration
if opts.verbose
    fprintf(['Starting censored SGL fit (FISTA): n=%d, p=%d, |U|=%d, |C|=%d, ' ...
        'c=%.6g, rho=%.6g, lambda1=%.6g, lambda2=%.6g\n'], ...
        n, p, numel(Uidx), numel(Cidx), c, rho, lambda1, lambda2);
end

for k = 1:opts.maxIter
    betaOld = beta;
    bOld = b;

    % Gradient of smooth part at extrapolated point
    [gBeta, gb, fY] = compute_grad_smooth(AU, yU, AC, c, rho, yBeta, yb);

    % ===== Backtracking line search =====
    t = step;
    accepted = false;
    for bt = 1:200
        % Proximal gradient step: prox_{SGL} applied to gradient descent direction
        betaCand = prox_sgl(yBeta - t * gBeta, t, lambda1, lambda2, groups);
        bCand = yb - t * gb;

        dBeta = betaCand - yBeta;
        db = bCand - yb;

        fCand = compute_f_smooth(AU, yU, AC, c, rho, betaCand, bCand);
        % Quadratic upper bound condition from Nesterov proximal gradient
        quadUpper = fY + gBeta.' * dBeta + gb * db + ...
            (1 / (2 * t)) * (norm(dBeta, 2)^2 + db^2);

        if fCand <= quadUpper + 1e-12 * max(1, abs(fY))
            accepted = true;
            break;
        end

        t = t * opts.backtrackBeta;
        if t <= eps
            break;
        end
    end

    if ~accepted
        exitType = 'lineSearchFailed';
        break;
    end

    % Accept the candidate solution
    beta = betaCand;
    b = bCand;
    step = t;

    smoothVal = fCand;
    % Compute full objective (smooth + nonsmooth penalties)
    obj = compute_objective_total(AU, yU, AC, c, rho, beta, b, lambda1, lambda2, groups);

    iters = k;
    objHist(k + 1) = obj;
    smoothHist(k + 1) = smoothVal;
    stepHist(k) = t;

    relChange = norm([beta - betaOld; b - bOld], 2) / max(1, norm([betaOld; bOld], 2));
    relChangeHist(k) = relChange;

    objPrev = objHist(k);
    relObj = abs(obj - objPrev) / max(1, abs(objPrev));
    relObjHist(k) = relObj;

    if opts.verbose
        fprintf('Iter %4d | obj = %.10e | relChange = %.3e | relObj = %.3e | step = %.3e | nnz = %d\n', ...
            k, obj, relChange, relObj, t, nnz(beta));
    end

    if relChange <= opts.tol || relObj <= opts.tol
        converged = true;
        if relChange <= opts.tol && relObj <= opts.tol
            exitType = 'tol_relChange_and_relObj';
        elseif relChange <= opts.tol
            exitType = 'tol_relChange';
        else
            exitType = 'tol_relObj';
        end
        break;
    end

    tFistaNew = (1 + sqrt(1 + 4 * tFista^2)) / 2;

    % ===== FISTA momentum update =====
    momentum = (tFista - 1) / tFistaNew;
    yBeta = beta + momentum * (beta - betaOld); % extrapolate beta with momentum
    yb = b + momentum * (b - bOld); % extrapolate b with momentum
    tFista = tFistaNew;
end

if iters == 0
    objHist = objHist(1);
    smoothHist = smoothHist(1);
    relChangeHist = nan(1, 1);
    relObjHist = nan(1, 1);
    stepHist = step;
else
    objHist = objHist(1:(iters + 1));
    smoothHist = smoothHist(1:(iters + 1));
    relChangeHist = relChangeHist(1:iters);
    relObjHist = relObjHist(1:iters);
    stepHist = stepHist(1:iters);
end

betaStd = beta;

%% ========== OUTPUT BACK-TRANSFORMATION ==========
% No internal centering/scaling; outputs already match input feature scale.

betaOrig = betaStd;
bStd = b;
bOrig = bStd;
intercept = bOrig;

nnzPerGroup = zeros(numel(groups), 1);
groupAllZero = false(numel(groups), 1);
nonzeroMaskPerGroup = cell(numel(groups), 1);
for g = 1:numel(groups)
    idx = groups{g};
    nzMask = abs(betaStd(idx)) > 0;
    nonzeroMaskPerGroup{g} = nzMask;
    nnzPerGroup(g) = nnz(nzMask);
    groupAllZero(g) = (nnzPerGroup(g) == 0);
end

%% ========== POST-OPTIMIZATION ANALYSIS ==========
% Compute sparsity patterns, nonzeros per group, etc.

stats = struct();
% Objective and convergence traces
stats.objHist = objHist;
stats.smoothHist = smoothHist;
stats.relChangeHist = relChangeHist;
stats.relObjHist = relObjHist;
stats.stepHist = stepHist;

% Sparsity and solution info
stats.nnz = nnz(betaStd);
stats.nnzPerGroup = nnzPerGroup;
stats.nonzeroMaskPerGroup = nonzeroMaskPerGroup;
stats.groupAllZero = groupAllZero;

% Coefficients and intercepts (both spaces)
stats.betaStd = betaStd;
stats.betaOrig = betaOrig;
stats.intercept = intercept;
stats.bStd = bStd;
stats.bOrig = bOrig;

% Transformation parameters (for reproducibility and manual verification)
stats.groups = groups;

% Censoring and penalty parameters
stats.c = c;
stats.rho = rho;
stats.Uidx = Uidx;
stats.Cidx = Cidx;
stats.Ucount = numel(Uidx);
stats.Ccount = numel(Cidx);

% Convergence info
stats.converged = converged;
stats.iters = iters;
stats.lambda1 = lambda1;
stats.lambda2 = lambda2;
stats.exitType = exitType;

fprintf('\n===== CENSORED SGL FIT SUMMARY =====\n');
fprintf('Iterations: %d\n', iters);
fprintf('Converged: %d\n', converged);
fprintf('Stop reason: %s\n', exitType);
fprintf('Uncensored |U| = %d, censored |C| = %d\n', numel(Uidx), numel(Cidx));
fprintf('c = %.12g, rho = %.12g\n', c, rho);
fprintf('Final objective: %.12g\n', objHist(end));
fprintf('Overall nonzeros: %d / %d\n', nnz(betaStd), p);
for g = 1:numel(groups)
    idx = groups{g};
    fprintf('Group %d (cols %d:%d) number of non zeros = %d / %d, all-zero = %d\n', ...
        g, idx(1), idx(end), nnzPerGroup(g), numel(idx), groupAllZero(g));
end

if opts.verbose
    figure('Name', 'Censored SGL Objective History', 'Color', 'w');
    generations = 0:(numel(objHist) - 1);
    plot(generations, objHist, '-o', 'LineWidth', 1.25, 'MarkerSize', 4);
    xlabel('Generation');
    ylabel('Objective value');
    title(['Censored Sparse Group Lasso Objective History of ', opts.muscleName], 'Interpreter','none');
    grid on;
end

% Quick sanity checks (manual):
% 1) If C is empty (or rho=0), this reduces to squared-loss + intercept + SGL.
% 2) stats.objHist should generally decrease with backtracking-enabled updates.

% % test
% z1 = A * betaOrig + intercept;
% z2 = AStd * betaStd + stats.bStd;
% fprintf('max abs diff in linear predictor = %.3g\n', max(abs(z1-z2)));
% yhat1 = max(z1, c);
% yhat2 = max(z2, c);
% fprintf('max abs diff in censored prediction = %.3g\n', max(abs(yhat1-yhat2)));

end

%% functions

function opts = apply_default_opts(opts)
% Apply default option values and validate user inputs.
%
% BACKWARD COMPATIBILITY HANDLING:
if isfield(opts, 'maxIterOuter') && (~isfield(opts, 'maxIter') || isempty(opts.maxIter))
    opts.maxIter = opts.maxIterOuter;
end
if isfield(opts, 'tolOuter') && (~isfield(opts, 'tol') || isempty(opts.tol))
    opts.tol = opts.tolOuter;
end

defaults = struct();
% Censoring and regularization parameters
defaults.c = 0.01; % censor floor
defaults.rho = [];
% Optimization parameters
defaults.maxIter = 20000;
defaults.tol = 1e-8;
defaults.verbose = true;
% FISTA + backtracking parameters
defaults.backtrackBeta = 0.5;
defaults.L0 = 1.0;
% Warm start (in original space)
defaults.beta0 = []; % first guess of beta (in original space) for warm start
defaults.b0 = 0; % first guess of intercept (in original space) for warm start
% Numerical tolerance for censoring partition
defaults.tolCensor = 1e-12;

% Fill in missing options with defaults
fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}) || isempty(opts.(fn{i}))
        opts.(fn{i}) = defaults.(fn{i});
    end
end

% Set rho default based on problem hardness if not provided
if isempty(opts.rho)
    % rho controls hardness of A_C*beta + b <= c; larger means harder enforcement.
    % 40-70% of data often censored, so rho=100 provides reasonable constraint enforcement
    opts.rho = 100;
end

end

function groups = group_indices(p)
% Define group structure for sparse group lasso regularization.
% Groups capture domain structure: muscle activations, post-activations, and other features.
%
% Input: p = size(A,2) = total number of features
% Output: groups = cell array of 3 groups, each containing column indices
%
% Assumes: p = 2*nMus + 8 (e.g., 26 = 2*9 + 8 for 9 muscles)
%   Group 1: muscle activations (columns 1 to nMus)
%   Group 2: post-activation effects (columns nMus+1 to 2*nMus)
%   Group 3: other features (columns 2*nMus+1 to 2*nMus+8)

nMus = (p - 8) / 2;
if nMus < 1 || abs(nMus - round(nMus)) > eps(max(1, nMus))
    error('sgl_fit:InvalidPredictorCount', ...
        'size(A,2) must satisfy p = 2*#mus + 8. Got p = %d.', p);
end
nMus = round(nMus);

groups = {
    1:nMus                          % group 1: fiber length feedbacks
    (nMus + 1):(2 * nMus)          % group 2: mtu force feedbacks
    (2 * nMus + 1):(2 * nMus + 8)  % group 3: kinematics / proprioception feedbacks
    };
end

function [gradBeta, gradb, fVal] = compute_grad_smooth(AU, yU, AC, c, rho, beta, b)
% Compute gradient of smooth (non-regularized) part of the objective.
%
% Smooth objective:
%   f(beta, b) = 0.5*||yU - (AU*beta + b)||_2^2 + 0.5*rho*||(AC*beta+b-c)_+||_2^2
%
% Outputs:
%   gradBeta: p x 1 gradient w.r.t. beta
%   gradb: scalar gradient w.r.t. b
%   fVal: current value of f(beta, b)

% Uncensored residual and gradient contributions
rU = AU * beta + b - yU;
% Censored term: only the positive part is active
tC = AC * beta + b - c;
pC = max(tC, 0);

% Gradient contribution from uncensored and censored terms
gradBeta = AU.' * rU + rho * (AC.' * pC);
gradb = sum(rU) + rho * sum(pC);
% Function value
fVal = 0.5 * (rU.' * rU) + 0.5 * rho * (pC.' * pC);
end

function f = compute_f_smooth(AU, yU, AC, c, rho, beta, b)
% Compute value of smooth objective (wrapper for single value query).
rU = AU * beta + b - yU;
tC = AC * beta + b - c;
pC = max(tC, 0);
f = 0.5 * (rU.' * rU) + 0.5 * rho * (pC.' * pC);
end

function obj = compute_objective_total(AU, yU, AC, c, rho, beta, b, lambda1, lambda2, groups)
% Compute total objective: smooth part + SGL penalties (L1 + group-L2).
%
% Total objective:
%   f(beta,b) + lambda1*sum_g ||beta_g||_2 + lambda2*||beta||_1
% where f is the smooth part and penalties are the non-smooth regularization.

f = compute_f_smooth(AU, yU, AC, c, rho, beta, b);
obj = f + lambda1 * group_l2_sum(beta, groups) + lambda2 * norm(beta, 1);
end

function s = group_l2_sum(beta, groups)
% Compute sum of L2 norms of beta within each group: sum_g ||beta_g||_2.
s = 0;
for g = 1:numel(groups)
    s = s + norm(beta(groups{g}), 2);
end
end

function betaProx = prox_sgl(v, t, lambda1, lambda2, groups)
% Proximal operator for sparse group LASSO (SGL) regularization.
%
% Solves: argmin_beta  0.5*||beta - v||^2 + t*lambda1*sum_g||beta_g||_2 + t*lambda2*||beta||_1
%
% Algorithm (per-group soft-thresholding followed by group shrinkage):
% For each group g with indices idx:
%   1. Soft-threshold v(idx) by t*lambda2: u = soft_threshold(v(idx), t*lambda2)
%   2. Compute group L2 norm: nu = ||u||_2
%   3. Shrink by group penalty: shrink = max(1 - t*lambda1/nu, 0)
%   4. Update: betaProx(idx) = shrink * u
%
% This implements the composition of L1 and group-L2 proximal operators.

betaProx = zeros(size(v));
for g = 1:numel(groups)
    idx = groups{g};
    u = soft_threshold(v(idx), t * lambda2);  % L1 shrinkage per group
    nu = norm(u, 2);
    if nu == 0
        betaProx(idx) = 0;                   % entire group stays zero
    else
        shrink = max(1 - (t * lambda1) / nu, 0);  % group-L2 shrinkage factor
        betaProx(idx) = shrink * u;          % only shrink if norm survives
    end
end
end

function y = soft_threshold(v, lam)
% Soft-thresholding (proximal operator for L1 regularization).
%
% Solves: argmin_y  0.5*||y - v||^2 + lam*||y||_1
% Solution: y = sign(v) . * max(|v| - lam, 0)  (element-wise)
%
% Application: shrinks small coefficients to zero when |v_i| < lam.

y = sign(v) .* max(abs(v) - lam, 0);
end