classdef CMAES_optimization < Base_optimization
    % Class name: CMAES_optimization (inherits from Base_optimization)
    % Description: Optimize parameters with the (mu/mu_w, lambda)-CMA-ES
    %     ( Hansen, 2016, "The CMA Evolution Strategy: A Tutorial", Algorithm 4 ).
    %     Active negative-weighted rank-mu update is enabled by default
    %     (a.k.a. aCMA-ES), and IPOP-CMA-ES restart is supported via the
    %     restart() method.
    %
    % PROPERTY TREE
    % - optParaNum                : int (= N)
    % - core
    %   |- xmean        : Nx1 ; current distribution mean
    %   |- sigma        : 1x1 ; global step size
    %   |- B            : NxN ; eigenvectors of C (orthonormal columns)
    %   |- diagD        : Nx1 ; STANDARD DEVIATIONS = sqrt(eigenvalues(C))
    %   |- C            : NxN ; covariance matrix
    %   |- invsqrtC     : NxN ; C^{-1/2} ; cached for active-CMA rescale
    %   |- pc , ps      : Nx1 ; evolution paths for C and sigma
    %   |- cc, cs, c1, cmu, damps, mueff, mueff_minus, chiN
    %   |- weights      : lambda x 1 ; positive (top mu) AND negative (bottom)
    %   |- mu, lambda   : ints
    % - recordOfBig3            : 3x1 struct  ( fit / arx / arz )   <- record only
    % - recordForBestParticle   : 1x1 struct  ( fit / arx )          <- record only
    % - arx , arz               : NxLambda    ; samples of current generation
    % - gMax                    : int         ; max #generations
    % - nParticles              : Nworkersx1  ; particles distributed per worker
    %
    % - lambda0, sigma0, xmean0     : initial values cached for restart
    % - popsizeFactor, restartCount : IPOP bookkeeping
    % - eigenEval, counteval        : lazy eigen-decomposition bookkeeping
    %
    % METHODS
    % - CMAES_optimization(initPara, sigma, Nworkers)
    % - [arx, arz] = generate_parameters(g)
    % - update_elite_fit(fits, arx, arz)         % updates RECORDS only
    % - update_core(fits, g)                     % standard CMA-ES update
    % - restart()                                % IPOP-CMA-ES restart

    properties
        recordOfBig3              % 3x1 struct (fit/arx/arz) — RECORD ONLY
        recordForBestParticle     % 1x1 struct (fit/arx)     — RECORD ONLY
        arx                       % current-generation samples in x-space
        arz                       % current-generation samples in z-space
        gMax                      % maximum #generation
        nParticles                % how many particles each worker evaluates
        nWorkers                  % cached parallel-pool size

        % --- IPOP / restart bookkeeping ---
        lambda0                   % initial population size
        sigma0                    % initial step size
        xmean0                    % initial mean (refreshed to "best so far" on restart)
        popsizeFactor             % multiplied by 2 on every IPOP restart
        restartCount              % number of restarts that have occurred

        % --- lazy eigen-decomposition bookkeeping ---
        eigenEval                 % counteval at the last B/diagD refresh
        counteval                 % cumulative number of fitness evaluations
    end

    methods
        function obj = CMAES_optimization(initPara, sigma, Nworkers)
            % Construct optimizer and initialize all strategy parameters.
            initPara = initPara(:);
            obj.lambda0       = 4 + floor(3*log(length(initPara)));
            obj.sigma0        = sigma;
            obj.xmean0        = initPara;
            obj.popsizeFactor = 1;
            obj.restartCount  = 0;
            obj.nWorkers      = Nworkers;
            obj.gMax          = 1000;

            obj.initStrategy();
            obj.initRecords();
        end

        function initStrategy(obj)
            % (Re)build all lambda-dependent strategy parameters and reset
            % all dynamic state. Called from the constructor and from restart().
            nPara  = length(obj.xmean0);
            lambda = obj.lambda0 * obj.popsizeFactor;
            mu     = floor(lambda/2);

            % ---- raw weights for ALL lambda points (Hansen 2016 Eq. 49) ----
            weights_raw = log(lambda/2 + 0.5) - log(1:lambda)';     % decreasing
            mueff       = sum(weights_raw(1:mu))^2     / sum(weights_raw(1:mu).^2);
            mueff_minus = sum(weights_raw(mu+1:end))^2 / sum(weights_raw(mu+1:end).^2);

            % ---- learning rates ----
            cc        = (4 + mueff/nPara) / (nPara + 4 + 2*mueff/nPara);
            cs        = (mueff + 2)       / (nPara + mueff + 5);
            c1        = 2                 / ((nPara + 1.3)^2 + mueff);
            alpha_cov = 2;
            cmu       = min(1 - c1, alpha_cov * (mueff - 2 + 1/mueff) / ...
                                    ((nPara + 2)^2 + alpha_cov*mueff/2));
            damps     = 1 + 2*max(0, sqrt((mueff-1)/(nPara+1)) - 1) + cs;

            % ---- normalize weights (Hansen 2016 Eq. 53; active-CMA negatives) ----
            sum_pos = sum(weights_raw(weights_raw > 0));
            sum_neg = -sum(weights_raw(weights_raw < 0));
            alpha_mu_minus     = 1 + c1/cmu;
            alpha_mueff_minus  = 1 + 2*mueff_minus/(mueff + 2);
            alpha_posdef_minus = (1 - c1 - cmu) / (nPara * cmu);
            weights = weights_raw;
            weights(weights > 0) = weights(weights > 0) / sum_pos;
            weights(weights < 0) = weights(weights < 0) * ...
                min([alpha_mu_minus, alpha_mueff_minus, alpha_posdef_minus]) / sum_neg;

            % expected length of N(0,I)
            chiN = sqrt(nPara) * (1 - 1/(4*nPara) + 1/(21*nPara^2));

            % ---- core state (clean slate on every restart) ----
            obj.core             = struct();
            obj.core.xmean       = obj.xmean0;
            obj.core.sigma       = obj.sigma0;
            obj.core.cc          = cc;
            obj.core.cs          = cs;
            obj.core.c1          = c1;
            obj.core.cmu         = cmu;
            obj.core.mu          = mu;
            obj.core.mueff       = mueff;
            obj.core.mueff_minus = mueff_minus;
            obj.core.lambda      = lambda;
            obj.core.damps       = damps;
            obj.core.pc          = zeros(nPara, 1);
            obj.core.ps          = zeros(nPara, 1);
            obj.core.B           = eye(nPara);
            obj.core.diagD       = ones(nPara, 1);            % sqrt(eigenvalues(C))
            obj.core.C           = eye(nPara);
            obj.core.invsqrtC    = eye(nPara);
            obj.core.chiN        = chiN;
            obj.core.weights     = weights;

            obj.optParaNum = nPara;
            obj.eigenEval  = 0;
            obj.counteval  = 0;

            % distribute lambda across parallel workers
            obj.nParticles = floor(lambda / obj.nWorkers) * ones(obj.nWorkers, 1);
            extra = mod(lambda, obj.nWorkers);
            if extra > 0
                obj.nParticles(1:extra) = obj.nParticles(1:extra) + 1;
            end
        end

        function initRecords(obj)
            % Initialize Big-3 / best records. NOT touched by restart() so that
            % the global best survives every IPOP restart.
            nPara = obj.optParaNum;
            obj.recordOfBig3 = struct( ...
                'fit', {Inf, Inf, Inf}, ...
                'arx', {nan(nPara,1), nan(nPara,1), nan(nPara,1)}, ...
                'arz', {nan(nPara,1), nan(nPara,1), nan(nPara,1)});
            obj.recordForBestParticle = struct( ...
                'fit', Inf, 'arx', nan(nPara,1));
        end

        function [arx, arz] = generate_parameters(obj, g)
            % Sample lambda candidate solutions:
            %     z_k ~ N(0, I)
            %     x_k = xmean + sigma * B * diag(diagD) * z_k
            BD  = obj.core.B .* obj.core.diagD.';            % == B * diag(diagD)
            arz = randn(obj.optParaNum, obj.core.lambda);
            arx = obj.core.xmean + obj.core.sigma * BD * arz;

            % First generation: keep one elitist copy of the mean for safety.
            % NOTE: arz(:,1) stays as a random sample, so the algebraic
            % invariant arx = xmean + sigma*BD*arz is intentionally NOT
            % enforced for k=1 in g==0 (same convention as Hansen's
            % purecmaes.m). That column will be excluded from the recombination
            % unless it ends up in the top-mu by fitness.
            if g == 0
                arx(:, 1) = obj.core.xmean;
            end

            % ---------------------------------------------------------------
            % Box-constraint clipping is INTENTIONALLY DISABLED here so that
            % the invariant  arx = xmean + sigma*B*diag(diagD)*arz  holds
            % exactly. If clipping is re-enabled in the future, you MUST
            % recompute the corresponding  arz_clipped =
            %     diag(1./diagD)*B'*(arx_clipped - xmean)/sigma
            % before passing arz into update_core(); otherwise the rank-one
            % and rank-mu updates will be evaluated in an inconsistent
            % coordinate frame.
            % ---------------------------------------------------------------
            % for i = 1 : obj.core.lambda
            %     arx(:, i) = obj.fix_parameters(arx(:, i));
            % end

            obj.arx = arx;
            obj.arz = arz;
        end

        function update_elite_fit(obj, fits, arx, arz)
            % RECORD ONLY — Big-3 / best are tracked for monitoring and IPOP
            % warm-restart, but they are NEVER injected into the selection
            % pool used by update_core (that would break the (mu/mu_w, lambda)
            % non-elitist contract of CMA-ES).
            for k = 1 : length(fits)
                fit = fits(k);
                if fit < obj.recordOfBig3(end).fit
                    obj.recordOfBig3(end).fit = fit;
                    obj.recordOfBig3(end).arx = arx(:, k);
                    obj.recordOfBig3(end).arz = arz(:, k);
                    [~, idx] = sort([obj.recordOfBig3(1).fit, ...
                                     obj.recordOfBig3(2).fit, ...
                                     obj.recordOfBig3(3).fit]);
                    obj.recordOfBig3 = obj.recordOfBig3(idx);

                    if fit < obj.recordForBestParticle.fit
                        obj.recordForBestParticle.fit = fit;
                        obj.recordForBestParticle.arx = arx(:, k);
                    end
                    disp(['fitness of elite particles: ', num2str([ ...
                        obj.recordOfBig3(1).fit, ...
                        obj.recordOfBig3(2).fit, ...
                        obj.recordOfBig3(3).fit])]);
                end
            end
        end

        function update_core(obj, fits, g)
            % Standard CMA-ES update over the CURRENT generation only.
            % References: Hansen 2016 Tutorial, Algorithm 4 (active CMA).
            nPara  = obj.optParaNum;
            lambda = obj.core.lambda;
            mu     = obj.core.mu;
            obj.counteval = obj.counteval + lambda;

            % ---- 1) selection : sort the lambda samples by fitness ----
            [~, idx]   = sort(fits);
            arx_sorted = obj.arx(:, idx);
            arz_sorted = obj.arz(:, idx);

            % ---- 2) recombination : weighted mean of TOP mu only ----
            wpos          = obj.core.weights(1:mu);
            xmean_new     = arx_sorted(:, 1:mu) * wpos;
            zmean         = arz_sorted(:, 1:mu) * wpos;
            obj.core.xmean = xmean_new;

            % ---- 3) cumulation : evolution paths ps and pc ----
            cs = obj.core.cs;  cc = obj.core.cc;
            obj.core.ps = (1 - cs) * obj.core.ps + ...
                sqrt(cs * (2 - cs) * obj.core.mueff) * (obj.core.B * zmean);

            hsig = norm(obj.core.ps) / ...
                   sqrt(1 - (1 - cs)^(2*(g + 1))) / obj.core.chiN ...
                   < 1.4 + 2/(nPara + 1);

            BD = obj.core.B .* obj.core.diagD.';
            obj.core.pc = (1 - cc) * obj.core.pc + ...
                hsig * sqrt(cc * (2 - cc) * obj.core.mueff) * (BD * zmean);

            % ---- 4) covariance matrix : rank-one + active rank-mu ----
            c1   = obj.core.c1;
            cmu  = obj.core.cmu;
            w    = obj.core.weights;                      % includes negatives

            % active-CMA rescaling for the NEGATIVE weights:
            %     wo(i) = w(i)                                , w(i) >= 0
            %           = w(i) * N / || C^{-1/2} y_i ||^2     , w(i) <  0
            % Note: || C^{-1/2} y_i ||^2 = || z_i ||^2 (B is orthonormal).
            wo     = w;
            negIdx = find(w < 0);
            for j = 1 : numel(negIdx)
                k     = negIdx(j);
                z     = arz_sorted(:, k);
                wo(k) = w(k) * nPara / max(sum(z.^2), eps);
            end

            delta_hsig = (1 - hsig) * cc * (2 - cc);
            ymat       = BD * arz_sorted;                  % y_i = B*D*z_i
            obj.core.C = (1 + c1*delta_hsig - c1 - cmu*sum(w)) * obj.core.C ...
                       + c1  * (obj.core.pc * obj.core.pc') ...
                       + cmu * (ymat * diag(wo) * ymat');

            % ---- 5) global step size sigma ----
            obj.core.sigma = obj.core.sigma * ...
                exp((cs / obj.core.damps) * ...
                    (norm(obj.core.ps) / obj.core.chiN - 1));

            % ---- 6) lazy eigen-decomposition of C (Hansen 2016, eq. 50) ----
            if obj.counteval - obj.eigenEval > lambda / (c1 + cmu) / nPara / 10
                obj.eigenEval = obj.counteval;
                obj.core.C    = triu(obj.core.C) + triu(obj.core.C, 1)';   % enforce symmetry
                [B, Dmat]     = eig(obj.core.C);
                d             = diag(Dmat);
                if any(d <= 0) || ~all(isfinite(d))
                    d = max(d, eps);                       % numerical guard
                end
                obj.core.B        = B;
                obj.core.diagD    = sqrt(d);
                obj.core.invsqrtC = B * diag(1 ./ sqrt(d)) * B';
            end
        end

        function restart(obj)
            % IPOP-CMA-ES restart: double lambda, reset all dynamic state,
            % keep recordOfBig3 / recordForBestParticle. The new mean is
            % seeded with the best-so-far parameter vector (warm start).
            obj.popsizeFactor = obj.popsizeFactor * 2;
            obj.restartCount  = obj.restartCount + 1;

            if isfinite(obj.recordForBestParticle.fit)
                obj.xmean0 = obj.recordForBestParticle.arx;
            end

            obj.initStrategy();    % rebuild lambda-dependent params + reset state

            fprintf('[CMA-ES] IPOP restart #%d : lambda=%d, sigma=%.3g\n', ...
                obj.restartCount, obj.core.lambda, obj.core.sigma);
        end

        function sigma_boost(obj, factor)
            % Multiplicative kick of the global step size, used as a soft
            % "wake-up" before triggering a full IPOP restart when fitness
            % improvement stalls. Default factor = 2.
            if nargin < 2 || isempty(factor)
                factor = 2.0;
            end
            obj.core.sigma = obj.core.sigma * factor;
            fprintf('[CMA-ES] sigma boosted x%.2f -> %.3g (stall recovery)\n', ...
                factor, obj.core.sigma);
        end
    end
end


