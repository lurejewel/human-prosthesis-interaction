function [lambda2Max, b0, info] = compute_lambda2_max_censored(X, y, c, rho, tolCensor, sampleWeight)
%COMPUTE_LAMBDA2_MAX_CENSORED Compute lambda2 max for censored LASSO.
%   [lambda2Max, b0, info] = compute_lambda2_max_censored(X, y, c, rho, tolCensor)
%   [lambda2Max, b0, info] = compute_lambda2_max_censored(X, y, c, rho, tolCensor, sampleWeight)
%   returns the smallest lambda2 such that beta=0 is optimal for the
%   censored least-squares objective used in sgl_fit.m (without 1/n scaling).
%
%   Inputs:
%     X, y, c, rho, tolCensor  - as in sgl_fit.m.
%     sampleWeight             - optional nonnegative n x 1 per-sample
%                                weights (default [] = all ones).  Pass the
%                                SAME vector that is given to sgl_fit.m as
%                                opts.sampleWeight; the weighted gradient
%                                then matches sgl_fit's weighted objective
%                                (uncensored squared loss + censored hinge).
%
%   KNOWN INCONSISTENCY (documented, not silently changed): the production
%   caller fit_adaptive_lasso_core.m calls this function WITHOUT
%   sampleWeight, while it passes opts.sampleWeight = wPhase to sgl_fit.m.
%   With non-uniform or zero (seam-excluded) weights the unweighted
%   lambda2Max is therefore only an approximation of the weighted
%   objective's null gradient.  The optional parameter is provided so the
%   caller can pass the same wPhase to remove the mismatch.  Existing
%   callers are unaffected: the default (all ones) is bit-identical to the
%   previous unweighted computation.

if nargin < 5 || isempty(tolCensor)
    tolCensor = 1e-12;
end

if ~isnumeric(X) || ~ismatrix(X)
    error('compute_lambda2_max_censored:BadDesign', 'X must be a numeric 2D matrix.');
end
if ~isnumeric(y)
    error('compute_lambda2_max_censored:BadResponse', 'y must be numeric.');
end
if ~isscalar(c) || ~isfinite(c)
    error('compute_lambda2_max_censored:BadC', 'c must be a finite scalar.');
end
if ~isscalar(rho) || ~isfinite(rho)
    error('compute_lambda2_max_censored:BadRho', 'rho must be a finite scalar.');
end
if ~isscalar(tolCensor) || ~isfinite(tolCensor) || tolCensor < 0
    error('compute_lambda2_max_censored:BadTolCensor', 'tolCensor must be a finite nonnegative scalar.');
end

y = y(:);
[n, p] = size(X);
if numel(y) ~= n
    error('compute_lambda2_max_censored:BadResponseSize', 'y must have %d rows to match X.', n);
end

% ---- optional sample weights (default: all ones = previous behavior) ----
if nargin < 6 || isempty(sampleWeight)
    sw = ones(n, 1);
else
    if ~isnumeric(sampleWeight) || ~isvector(sampleWeight)
        error('compute_lambda2_max_censored:BadSampleWeight', ...
            'sampleWeight must be a numeric n x 1 vector.');
    end
    sw = double(sampleWeight(:));
    if numel(sw) ~= n
        error('compute_lambda2_max_censored:BadSampleWeight', ...
            'sampleWeight must have %d entries, got %d.', n, numel(sw));
    end
    if any(~isfinite(sw)) || any(sw < 0)
        error('compute_lambda2_max_censored:BadSampleWeight', ...
            'sampleWeight must be finite and nonnegative.');
    end
end

Uidx = find(y > c + tolCensor);
Cidx = find(y <= c + tolCensor);

if isempty(Uidx)
    b0 = c;
    lambda2Max = 0;
    gradBetaAtZero = zeros(p, 1);
    info = struct();
    info.Uidx = Uidx;
    info.Cidx = Cidx;
    info.Ucount = numel(Uidx);
    info.Ccount = numel(Cidx);
    info.b0 = b0;
    info.gradBetaAtZero = gradBetaAtZero;
    info.lambda2Max = lambda2Max;
    info.smoothAtZero = 0;
    info.sampleWeighted = any(sw ~= 1);
    info.sumWeightU = 0;
    info.sumWeightC = sum(sw(Cidx));
    return;
end

nU = numel(Uidx);
nC = numel(Cidx);
swU = sw(Uidx);
swC = sw(Cidx);
sumWU = sum(swU);
sumWC = sum(swC);
if sumWU + rho * sumWC <= 0
    error('compute_lambda2_max_censored:NoPositiveWeight', ...
        'sampleWeight has no positive weight on the uncensored/censored samples.');
end

b0 = (sum(swU .* y(Uidx)) + rho * sumWC * c) / (sumWU + rho * sumWC);

rU = b0 - y(Uidx);
pC = max(b0 - c, 0);

gradBetaAtZero = X(Uidx, :).' * (swU .* rU);
if nC > 0 && pC > 0
    gradBetaAtZero = gradBetaAtZero + ...
        rho * (X(Cidx, :).' * (swC .* (pC * ones(nC, 1))));
end

lambda2Max = max(abs(gradBetaAtZero));
if ~isfinite(lambda2Max)
    error('compute_lambda2_max_censored:BadLambda2Max', 'lambda2Max is not finite.');
end

smoothAtZero = 0.5 * (rU.' * (swU .* rU));
if nC > 0 && pC > 0
    smoothAtZero = smoothAtZero + 0.5 * rho * (pC ^ 2 * sumWC);
end

info = struct();
info.Uidx = Uidx;
info.Cidx = Cidx;
info.Ucount = nU;
info.Ccount = nC;
info.b0 = b0;
info.gradBetaAtZero = gradBetaAtZero;
info.lambda2Max = lambda2Max;
info.smoothAtZero = smoothAtZero;
info.sampleWeighted = any(sw ~= 1);
info.sumWeightU = sumWU;
info.sumWeightC = sumWC;
end