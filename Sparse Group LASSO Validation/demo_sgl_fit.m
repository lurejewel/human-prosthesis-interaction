% DEMO_SGL_FIT Example usage for Sparse Group Lasso solver.
% Assumes X and y are already in workspace.

if ~exist('X', 'var') || ~exist('y', 'var')
    error('demo_sgl_fit:MissingData', 'X and y must exist in the workspace.');
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
