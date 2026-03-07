% Sparse Group LASSO runner script.
% Requires X (10000x18 double) and y (10000x1 double) in workspace.

if ~exist('X', 'var') || ~exist('y', 'var')
	error('SparseGroupLASSO:MissingData', 'X and y must exist in the workspace.');
end

lambda1 = 0.5;
lambda2 = 0.1;

opts = struct();
opts.maxIterOuter = 200;
opts.maxIterInner = 200;
opts.tolOuter = 1e-6;
opts.tolInner = 1e-6;
opts.tol1D = 1e-10;
opts.maxIter1D = 50;
opts.verbose = true;

[beta, intercept, betaStd, stats] = sgl_fit(X, y, lambda1, lambda2, opts);

fprintf('\nScript run complete. Final objective = %.12g\n', stats.objHist(end));
fprintf('Converged = %d, iterations = %d\n', stats.converged, stats.iters);
fprintf('Intercept = %.12g\n', intercept);

disp('beta (original scale):');
disp(beta);
disp('betaStd (standardized scale):');
disp(betaStd);
