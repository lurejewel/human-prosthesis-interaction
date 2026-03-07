classdef CMAES_optimization < Base_optimization
    % Class name: CMAES_optimization (inherented from Base_optimization)
    % Description: optimize parameters though the Covariance Matrix Adaptation
    % Evolutionary Strategy (CMA-ES). Parameters are muscle-reflex coefficients
    % in the musculoskeletal model.
    % PROPERTY TREE:
    % - optParaNum: int
    % - optParaRecord
    % | - bestFit: 1xK (K <= #iterations)
    % | - bestPara: NxK
    % | - bestIter: 1xK
    % - core
    % | - xmean: Nx1 (N = #dimensions)
    % | - zmean: Nx1
    % | - arx: Nxlambda (lambda = #particles)
    % | - arz: Nxlambda
    % | - sigma: double (usually 0.01 or 0.02)
    % | - cc: double
    % | - cs: double
    % | - c1: double
    % | - cmu: double
    % | - mu: int
    % | - mueff: double
    % | - lambda: int
    % | - damps: double
    % | - pc: Nx1
    % | - ps: Nx1
    % | - B: NxN
    % | - C: NxN
    % | - D: NxN
    % | - chiN: double
    % | - noise: double (usually 0.05)
    % | - weights: : mux1
    % - fitnessForCurrentIteration: lambdax1
    % - fitnessWithNoiseForCurrentIteration: lambdax1
    % METHODS:
    % - CMAES_optimization(<Nx1>para, <double>sigma)
    % - [<Nxlambda>xvec, <Nxlambda>zvec] = generate_parameters(<obj, <int>iteration)
    % - xxx = update(obj)
    % - <Nx1>fixedPara = fix_parameters(obj, <Nx1>para)
    % - <double>x_ = min_max(obj, <double>x, <double>lowerTh, <double>upperTh)

    properties
        % optParaNum % inherented from Parent Class
        % optParaRecord % inherented from Parent Class 这个目前没有被用到？
        % core % inherented from Parent Class

        % fitnesses
        recordOfBig3 % 3x1 struct
        recordForBestParticle % 1x1 struct
        arx
        arz
        gMax % maximum #generation possible for optimization
        nParticles % #particle to be simulated and evaluated within each worker (batch)

    end

    methods
        function obj = CMAES_optimization(initPara, sigma, Nworkers)
            % Name: CMAES_optimization
            % Description: construction method.

            % read input parameters
            obj.core.xmean = initPara; % parameters (treated as mean value of a vector of random variables) to be optimized
            obj.core.sigma = sigma; % describes the step length of the optimizaiton

            % calculate hyper parameters
            % 1. selection of points
            nPara = length(initPara);
            lambda = 4 + floor(3*log(nPara)); % population size [=14 for 29 pars; =13 for 28 pars]
            mu = lambda / 2;
            weights = log(mu+1 / 2) - log(1:mu)'; % recombination weights, assigned from big to small
            mu = floor(mu); % number of parents/points for recombination
            weights = weights / sum(weights); % normalize recombination weights array
            mueff = sum(weights)^2 / sum(weights.^2); % variance-effective size of mu
            % 2. adaptation
            cc = (4 + mueff/nPara) / (nPara + 4 + 2*mueff/nPara); % time constant for cumulation for C
            cs = (mueff+2) / (nPara+mueff+5); % t-const for cumulation for sigma control
            c1 = 2 / ((nPara+1.3)^2 + mueff); % learning rate for rank-one update of C
            cmu = 2 * (mueff - 2 + 1/mueff) / ((nPara+2)^2 + 2*mueff/2); % and for rank-mu update
            damps = 1 + 2 * max(0, sqrt((mueff-1)/(nPara+1))-1) + cs; % damping for sigma
            % 3. Initialize dynamic (internal) strategy parameters and constants
            pc = zeros(nPara,1); ps = zeros(nPara,1); % evolution paths for C and sigma
            B = eye(nPara); % B defines the coordinate system
            D = eye(nPara); % diagonal matrix D defines the scaling
            C = B * D * (B*D)'; % covariance matrix
            chiN = nPara^0.5*(1-1/(4*nPara)+1/(21*nPara^2)); % expected length of a N(0,I)

            % assign parameters
            obj.core.cc = cc;
            obj.core.cs = cs;
            obj.core.c1 = c1;
            obj.core.cmu = cmu;
            obj.core.mu = mu;
            obj.core.mueff = mueff;
            obj.core.lambda = lambda;
            obj.core.damps = damps;
            obj.core.pc = pc;
            obj.core.ps = ps;
            obj.core.B = B;
            obj.core.C = C;
            obj.core.D = D;
            obj.core.chiN = chiN;
            obj.core.noise = 0.05;
            obj.core.weights = weights;

            obj.optParaNum = nPara;
            obj.recordOfBig3 = struct('fit', {999,999,999},...  % 'para', {}, ...
                'arx', {zeros(nPara,1), zeros(nPara,1), zeros(nPara,1)}, ...
                'arz', {zeros(nPara,1), zeros(nPara,1), zeros(nPara,1)});
            obj.recordForBestParticle = struct('fit', 999, 'arx', []);

            obj.gMax = 1000;
            obj.nParticles = floor(lambda/Nworkers) * ones(Nworkers,1);
            obj.nParticles(1 : mod(lambda,Nworkers)) = obj.nParticles(1 : mod(lambda,Nworkers)) + 1;

        end

        function [arx, arz] = generate_parameters(obj, g)
            % Description: generate muscle reflex parameters with lambda
            % loops.
            % Output:
            % - [arx] an nPara x lambda matrix. arx is where the muscle
            % reflex parameters are placed in order:
            % [para1 | para2 | ... | paraN].

            % deviation/exploration according to the covariance matrix C
            arz = randn(obj.optParaNum, obj.core.lambda); % random initialization of noise
            arx = obj.core.xmean + obj.core.sigma * obj.core.B * obj.core.D * arz; % add mutation

            % for the 1st iteration: no noise added to the 1st particle
            if g == 0
                arx(:,1) = obj.core.xmean;
            end

            % clip the parameters
            for i = 1 : obj.core.lambda
                arx(:,i) = obj.fix_parameters(arx(:,i));
            end

            obj.arx = arx;
            obj.arz = arz;

        end

        function update_elite_fit(obj, fits, arx, arz)
            % fitnesses to be updated:
            % - obj.recordOfBig3: 3x1 struct
            % - obj.recordForBestParticle: 1x1 struct

            for k = 1 : length(fits)
                fit = fits(k);

                if fit < obj.recordOfBig3(end).fit % greater than the elite
                    obj.recordOfBig3(end).fit = fit;
                    obj.recordOfBig3(end).arx = arx(:,k);
                    obj.recordOfBig3(end).arz = arz(:,k);
                    % sort the big 3
                    [~, idx] =  sort([obj.recordOfBig3(1).fit, obj.recordOfBig3(2).fit, obj.recordOfBig3(3).fit]);
                    obj.recordOfBig3 = obj.recordOfBig3(idx);

                    if fit < obj.recordForBestParticle.fit % greater than the best
                        obj.recordForBestParticle.fit = fit;
                        obj.recordForBestParticle.arx = arx(:,k);

                    end
                    disp(['fitness of elite particles: ', num2str([obj.recordOfBig3(1).fit, obj.recordOfBig3(2).fit, obj.recordOfBig3(3).fit])]);
                end
            end


        end

        function update_core(obj, fits, g)
            % TODO: 后续sigma和noise可能会修改为与g（准确来说，是与优化停滞的代数）相关

            % update fits for the generation
            fitsWithNoiseForCurrentIter = fits + obj.core.noise * randn(height(fits), width(fits));
            fits_all = [fitsWithNoiseForCurrentIter', obj.recordOfBig3(1).fit, obj.recordOfBig3(2).fit, obj.recordOfBig3(3).fit];
            arx_all = [obj.arx, obj.recordOfBig3(1).arx, obj.recordOfBig3(2).arx, obj.recordOfBig3(3).arx];
            arz_all = [obj.arz, obj.recordOfBig3(1).arz, obj.recordOfBig3(2).arz, obj.recordOfBig3(3).arz];
            [~, idx] = sort(fits_all);
            obj.core.xmean = arx_all(:, idx(1:obj.core.mu)) * obj.core.weights;
            zmean = arz_all(:, idx(1:obj.core.mu)) * obj.core.weights;

            % cumulation: update evolution paths
            obj.core.ps = (1-obj.core.cs) * obj.core.ps + (sqrt(obj.core.cs*(2-obj.core.cs)*obj.core.mueff)) * (obj.core.B * zmean);
            hsig = norm(obj.core.ps)/sqrt(1-(1-obj.core.cs)^(2*g/obj.core.lambda))/obj.core.chiN < 1.4+2/(obj.optParaNum+1);
            obj.core.pc = (1-obj.core.cc)*obj.core.pc + hsig * sqrt(obj.core.cc*(2-obj.core.cc)*obj.core.mueff) * (obj.core.B*obj.core.D*zmean);

            % adapt covariance matrix C
            obj.core.C = (1-obj.core.c1-obj.core.cmu) * obj.core.C ... % regard old matrix
                + obj.core.c1 * (obj.core.pc*obj.core.pc' ... % plus rand one update
                + (1-hsig)*obj.core.cc*(2-obj.core.cc)*obj.core.C) ... % minor correction
                + obj.core.cmu ... % plus rank mu update
                * (obj.core.B * obj.core.D * arz_all(:, idx(1:obj.core.mu))) ...
                * diag(obj.core.weights) * (obj.core.B*obj.core.D*arz_all(:,idx(1:obj.core.mu)))';

            % adapt step-size sigma
            obj.core.sigma = obj.core.sigma * exp((obj.core.cs/obj.core.damps)*(norm(obj.core.ps)/obj.core.chiN-1));

            % update B & D from C
            obj.core.C = triu(obj.core.C) + triu(obj.core.C, 1)'; % enforce symmetry
            [B, D] = eig(obj.core.C); % eigen decomposition, B == normalized eigenvectors
            obj.core.B = B;
            obj.core.D = diag(sqrt(diag(D))); % D contains standard deviations now

        end

    end

end