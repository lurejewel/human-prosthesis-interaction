function [betaOrig, intercept, betaStd, stats] = sgl_fit(X, y, lambda1, lambda2, opts)
%SGL_FIT Sparse Group Lasso for linear regression via block coordinate descent.
%   [betaOrig, intercept, betaStd, stats] = sgl_fit(X, y, lambda1, lambda2, opts)
%   solves in standardized space:
%     min_beta 0.5*||yStd - XStd*beta||_2^2 + lambda1*sum_g ||beta_g||_2 + lambda2*||beta||_1
%
%   Group structure is inferred from p = size(X,2), where p = 3*#mus + 12:
%     group1: 1:2*#mus
%     group2: 2*#mus+1:3*#mus
%     group3: 3*#mus+1:3*#mus+12
%
%   Inputs:
%     X       - n x p predictor matrix
%     y       - n x 1 (or 1 x n) response vector
%     lambda1 - group-l2 penalty weight
%     lambda2 - l1 penalty weight
%     opts    - optional struct with fields:
%               maxIterOuter (default 2000)
%               maxIterInner (default 2000)
%               tolOuter     (default 1e-6)
%               tolInner     (default 1e-6)
%               tol1D        (default 1e-10)
%               maxIter1D    (default 50)
%               verbose      (default true)
%               beta0Std     (default []) warm start in standardized space
%
%   Outputs:
%     betaOrig  - p x 1 coefficients on original scale
%     intercept - scalar intercept on original scale
%     betaStd   - p x 1 coefficients on standardized scale
%     stats     - struct with optimization traces and metadata

    if nargin < 5
        opts = struct();
    end
    opts = apply_default_opts(opts);

    y = y(:);

    n = size(X, 1);
    p = size(X, 2);
    groups = group_indices(p);

    % ----- Preprocessing: center and scale X, center y -----
    muX = mean(X, 1);
    Xc = X - muX;
    sX = std(Xc, 1, 1); % population std, divide by n
    constMask = (sX == 0);
    sX(constMask) = 1;
    XStd = Xc ./ sX;

    muy = mean(y);
    yStd = y - muy;

    p = size(XStd, 2);
    xnorm2 = sum(XStd.^2, 1).';

    % ----- Initialization -----
    betaStd = zeros(p, 1);
    if ~isempty(opts.beta0Std)
        if ~isvector(opts.beta0Std) || numel(opts.beta0Std) ~= p
            error('sgl_fit:BadWarmStart', 'opts.beta0Std must be a %d-by-1 vector.', p);
        end
        betaStd = opts.beta0Std(:);
    end

    r = yStd - XStd * betaStd;

    objHist = zeros(opts.maxIterOuter, 1);
    relChangeHist = nan(opts.maxIterOuter, 1);
    residualNormHist = nan(opts.maxIterOuter, 1);
    converged = false;
    tiny = 1e-14;

    if opts.verbose
        fprintf('Starting SGL fit: n=%d, p=%d, lambda1=%.6g, lambda2=%.6g\n', n, p, lambda1, lambda2);
    end

    for iterOuter = 1:opts.maxIterOuter
        betaOld = betaStd;

        for g = 1:numel(groups)
            idx = groups{g};

            % Partial residual excluding current group.
            rg = r + XStd(:, idx) * betaStd(idx);

            % Group-level KKT screening.
            a = XStd(:, idx).' * rg;
            u = soft_threshold(a, lambda2);
            if norm(u, 2) <= lambda1
                betaStd(idx) = 0;
                r = rg;
                continue;
            end

            % Within-group coordinate descent.
            for iterInner = 1:opts.maxIterInner
                maxDelta = 0;
                for jLocal = 1:numel(idx)
                    j = idx(jLocal);
                    xj = XStd(:, j);
                    bjOld = betaStd(j);

                    % Coordinate-wise partial residual.
                    w = r + xj * bjOld;
                    c = xj.' * w;

                    if abs(c) <= lambda2 || xnorm2(j) <= tiny
                        bjNew = 0;
                    else
                        bg = betaStd(idx);
                        s2 = sum(bg.^2) - bjOld^2;
                        s2 = max(s2, 0);
                        bjNew = solve_1d_sgl(c, xnorm2(j), lambda1, lambda2, s2, opts);
                    end

                    if bjNew ~= bjOld
                        betaStd(j) = bjNew;
                        r = r + xj * (bjOld - bjNew);
                        maxDelta = max(maxDelta, abs(bjNew - bjOld));
                    end
                end

                if maxDelta <= opts.tolInner * max(1, norm(betaStd(idx), inf))
                    break;
                end
            end
        end

        obj = compute_objective(r, betaStd, lambda1, lambda2, groups);
        objHist(iterOuter) = obj;
        residualNormHist(iterOuter) = norm(r);

        relChange = norm(betaStd - betaOld, 2) / max(1, norm(betaOld, 2));
        relChangeHist(iterOuter) = relChange;

        if iterOuter == 1
            objRel = inf;
        else
            objPrev = objHist(iterOuter - 1);
            objRel = abs(obj - objPrev) / max(1, abs(objPrev));
        end

        if opts.verbose
            fprintf('Iter %3d | obj = %.10e | relBeta = %.3e | relObj = %.3e | nnz = %d\n', ...
                iterOuter, obj, relChange, objRel, nnz(betaStd));
        end

        if relChange <= opts.tolOuter || objRel <= opts.tolOuter
            converged = true;
            break;
        end
    end

    if converged
        iters = iterOuter;
    else
        iters = opts.maxIterOuter;
    end

    objHist = objHist(1:iters);
    relChangeHist = relChangeHist(1:iters);
    residualNormHist = residualNormHist(1:iters);

    % ----- Back-transform to original scale -----
    betaOrig = betaStd ./ sX(:);
    intercept = muy - muX * betaOrig;

    % ----- Group summaries -----
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

    stats = struct();
    stats.objHist = objHist;
    stats.relChangeHist = relChangeHist;
    stats.residualNormHist = residualNormHist;
    stats.nnzPerGroup = nnzPerGroup;
    stats.nonzeroMaskPerGroup = nonzeroMaskPerGroup;
    stats.groupAllZero = groupAllZero;
    stats.betaStd = betaStd;
    stats.betaOrig = betaOrig;
    stats.intercept = intercept;
    stats.muX = muX;
    stats.sX = sX;
    stats.muy = muy;
    stats.groups = groups;
    stats.constColMask = constMask;
    stats.converged = converged;
    stats.iters = iters;
    stats.lambda1 = lambda1;
    stats.lambda2 = lambda2;

    % ----- Required result display -----
    fprintf('\n===== SGL FIT SUMMARY =====\n');
    fprintf('Outer iterations: %d\n', iters);
    fprintf('Converged: %d\n', converged);
    fprintf('Final objective: %.12g\n', objHist(end));
    fprintf('Overall nonzeros: %d / %d\n', nnz(betaStd), p);
    for g = 1:numel(groups)
        idx = groups{g};
        fprintf('Group %d (cols %d:%d) number of non zeros = %d / %d, all-zero = %d\n', ...
            g, idx(1), idx(end), nnzPerGroup(g), numel(idx), groupAllZero(g));
    end

    if opts.verbose
        figure('Name', 'SGL Objective History', 'Color', 'w');
        plot(1:iters, objHist, '-o', 'LineWidth', 1.25, 'MarkerSize', 4);
        xlabel('Outer iteration');
        ylabel('Objective value');
        title('Sparse Group Lasso Objective History');
        grid on;
    end
end

function opts = apply_default_opts(opts)
    defaults = struct();
    defaults.maxIterOuter = 200;
    defaults.maxIterInner = 200;
    defaults.tolOuter = 1e-6;
    defaults.tolInner = 1e-6;
    defaults.tol1D = 1e-10;
    defaults.maxIter1D = 50;
    defaults.verbose = true;
    defaults.beta0Std = [];

    fn = fieldnames(defaults);
    for i = 1:numel(fn)
        if ~isfield(opts, fn{i}) || isempty(opts.(fn{i}))
            opts.(fn{i}) = defaults.(fn{i});
        end
    end
end

function groups = group_indices(p)
    nMus = (p - 12) / 3;
    if nMus < 1 || abs(nMus - round(nMus)) > eps(max(1, nMus))
        error('sgl_fit:InvalidPredictorCount', ...
            'size(X,2) must satisfy p = 3*#mus + 12. Got p = %d.', p);
    end
    nMus = round(nMus);

    groups = {
        1:(2 * nMus)
        (2 * nMus + 1):(3 * nMus)
        (3 * nMus + 1):(3 * nMus + 12)
    };
end

function obj = compute_objective(r, beta, lambda1, lambda2, groups)
    groupPenalty = 0;
    for g = 1:numel(groups)
        groupPenalty = groupPenalty + norm(beta(groups{g}), 2);
    end
    obj = 0.5 * (r.' * r) + lambda1 * groupPenalty + lambda2 * norm(beta, 1);
end

function y = soft_threshold(v, lam)
    y = sign(v) .* max(abs(v) - lam, 0);
end

function t = solve_1d_sgl(c, xnorm2, lambda1, lambda2, s2, opts)
% Solve:
% min_t 0.5*xnorm2*t^2 - c*t + lambda2*|t| + lambda1*sqrt(t^2 + s2)
% with monotone root solve on fixed sign branch.

    tiny = 1e-14;
    if xnorm2 <= tiny
        t = 0;
        return;
    end

    if s2 < tiny
        % Degenerate case: group-L2 part reduces to lambda1*|t|.
        t = soft_threshold(c, lambda1 + lambda2) / xnorm2;
        return;
    end

    if abs(c) <= lambda2
        t = 0;
        return;
    end

    sgn = sign(c);
    cabs = abs(c);

    f = @(q) xnorm2 * q - cabs + lambda2 + lambda1 * q ./ sqrt(q.^2 + s2);
    df = @(q) xnorm2 + lambda1 * s2 ./ (q.^2 + s2).^(3/2);

    qLo = 0;
    qHi = max(abs(soft_threshold(c, lambda2) / xnorm2), 1);
    while f(qHi) < 0
        qHi = 2 * qHi;
        if qHi > 1e12
            break;
        end
    end

    q = min(max(abs(soft_threshold(c, lambda2) / xnorm2), qLo), qHi);

    for k = 1:opts.maxIter1D
        fq = f(q);
        if abs(fq) <= opts.tol1D
            break;
        end

        qNew = q - fq / df(q);
        if ~isfinite(qNew) || qNew <= qLo || qNew >= qHi
            qNew = 0.5 * (qLo + qHi);
        end

        fNew = f(qNew);
        if fNew < 0
            qLo = qNew;
        else
            qHi = qNew;
        end

        if abs(qNew - q) <= opts.tol1D * max(1, q)
            q = qNew;
            break;
        end

        q = qNew;
    end

    t = sgn * q;
end
