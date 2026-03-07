% DEMO_SGL_FIT Example usage for Sparse Group Lasso solver.
% Assumes X (10000x18 double) and y (10000x1 double) are in workspace.

if ~exist('X', 'var') || ~exist('y', 'var')
    error('demo_sgl_fit:MissingData', 'X and y must exist in the workspace.');
end

if ~isa(X, 'double') || ~isa(y, 'double')
    error('demo_sgl_fit:TypeMismatch', 'X and y must be double.');
end

if ~isequal(size(X), [10000, 18])
    warning('demo_sgl_fit:UnexpectedSizeX', 'Expected X size [10000, 18], got [%d, %d].', size(X, 1), size(X, 2));
end

if ~isequal(size(y), [10000, 1])
    warning('demo_sgl_fit:UnexpectedSizeY', 'Expected y size [10000, 1], got [%d, %d].', size(y, 1), size(y, 2));
end

lambda1 = 0.5;
lambda2 = 0.1;
opts = struct('verbose', true);

[beta, intercept, betaStd, stats] = sgl_fit(X, y, lambda1, lambda2, opts);

fprintf('\nDemo complete. Final objective: %.12g\n', stats.objHist(end));
fprintf('Converged: %d in %d iterations\n', stats.converged, stats.iters);
fprintf('Intercept: %.12g\n', intercept);
disp('beta (original scale):');
disp(beta);
disp('betaStd (standardized scale):');
disp(betaStd);
