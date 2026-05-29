function [lambda2Max, b0, info] = compute_lambda2_max_censored(A, y, c, rho, tolCensor)
%COMPUTE_LAMBDA2_MAX_CENSORED Compute lambda2 max for censored LASSO.
%   [lambda2Max, b0, info] = compute_lambda2_max_censored(A, y, c, rho, tolCensor)
%   returns the smallest lambda2 such that beta=0 is optimal for the
%   censored least-squares objective used in sgl_fit.m (without 1/n scaling).

if nargin < 5 || isempty(tolCensor)
    tolCensor = 1e-12;
end

if ~isnumeric(A) || ndims(A) ~= 2
    error('compute_lambda2_max_censored:BadDesign', 'A must be a numeric 2D matrix.');
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
[n, p] = size(A);
if numel(y) ~= n
    error('compute_lambda2_max_censored:BadResponseSize', 'y must have %d rows to match A.', n);
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
    return;
end

nU = numel(Uidx);
nC = numel(Cidx);

b0 = (sum(y(Uidx)) + rho * nC * c) / (nU + rho * nC);

rU = b0 - y(Uidx);
pC = max(b0 - c, 0);

gradBetaAtZero = A(Uidx, :).'* rU;
if nC > 0 && pC > 0
    gradBetaAtZero = gradBetaAtZero + rho * (A(Cidx, :).'* (pC * ones(nC, 1)));
end

lambda2Max = max(abs(gradBetaAtZero));
if ~isfinite(lambda2Max)
    error('compute_lambda2_max_censored:BadLambda2Max', 'lambda2Max is not finite.');
end

smoothAtZero = 0.5 * (rU.' * rU);
if nC > 0 && pC > 0
    smoothAtZero = smoothAtZero + 0.5 * rho * (pC ^ 2 * nC);
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
end
