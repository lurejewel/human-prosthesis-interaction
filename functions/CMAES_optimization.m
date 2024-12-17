classdef CMAES_optimization < Base_optimization
    % Class name: CMAES_optimization (inherented from Base_optimization)
    % Description: optimize parameters though the Covariance Matrix Adaptation
    % Evolutionary Strategy (CMA-ES). Parameters are muscle-reflex coefficients
    % in the musculoskeletal model.
    % PROPERTY TREE:
    % - optPara
    % | - xmean: Nx1 (N = #dimensions)
    % | - zmean: Nx1
    % | - xvec: Nxlambda (lambda = #particles)
    % | - zvec: Nxlambda
    % - optParaNum: int
    % - optParaRecord
    % | - bestFit: 1xK (K <= #iterations)
    % | - bestPara: NxK
    % | - bestIter: 1xK
    % - hyperPara
    % | - sigma: double (usually 0.01 or 0.02)
    % | - stage: int (1 or 2)
    % - intrinsicPara
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
    % - CMAES_optimization(<Nx1>para, <double>sigma, <int>stage)
    % - [<Nxlambda>xvec, <Nxlambda>zvec] = generate_parameters(<obj, <int>iteration)
    % - xxx = update(obj)
    % - <Nx1>fixedPara = fix_parameters(obj, <Nx1>para)
    % - <double>x_ = min_max(obj, <double>x, <double>lowerTh, <double>upperTh)

    properties
        % optPara % inherented from Parent Class
        % optParaNum % inherented from Parent Class
        % optParaRecord % inherented from Parent Class
        % hyperPara % inherented from Parent Class
        % intrinsicPara % inherented from Parent Class
        fitnessForCurrentIteration
        fitnessWithNoiseForCurrentIteration
    end

    methods
        function obj = CMAES_optimization(optParaInput, sigma, stage)
            % Name: CMAES_optimization
            % Description: construction method.

            % read input parameters
            obj.optPara.xmean = optParaInput; % parameters (treated as mean value of a vector of random variables) to be optimized
            obj.hyperPara.sigma = sigma; % describes the step length of the optimizaiton
            obj.hyperPara.stage = stage; % this is a 2-stage optimization

            % calculate hyper parameters
            % 1. selection of points
            nPara = length(optParaInput);
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

            % assign hyper & intrinsic parameters
            obj.intrinsicPara.cc = cc;
            obj.intrinsicPara.cs = cs;
            obj.intrinsicPara.c1 = c1;
            obj.intrinsicPara.cmu = cmu;
            obj.intrinsicPara.mu = mu;
            obj.intrinsicPara.mueff = mueff;
            obj.intrinsicPara.lambda = lambda;
            obj.intrinsicPara.damps = damps;
            obj.intrinsicPara.pc = pc;
            obj.intrinsicPara.ps = ps;
            obj.intrinsicPara.B = B;
            obj.intrinsicPara.C = C;
            obj.intrinsicPara.D = D;
            obj.intrinsicPara.chiN = chiN;
            obj.intrinsicPara.noise = 0.05;
            obj.intrinsicPara.weights = weights;

            obj.optParaNum = nPara;
            obj.fitnessForCurrentIteration = nan(lambda, 1);
            obj.fitnessWithNoiseForCurrentIteration = nan(lambda, 1);

        end

        function [xvec, zvec] = generate_parameters(obj, iteration)
            % Name: generate_parameters
            % Description: generate muscle reflex parameters with lambda
            % loops.
            % Output:
            % - [paraVector] an nPara x lambda matrix. paraVector is
            % where the muscle reflex parameters are placed in order:
            % [paras1 | paras2 | ... | parasN].

            zvec = randn(obj.optParaNum, obj.intrinsicPara.lambda); % deviation/exploration according to the covariance matrix C
            xvec = zeros(obj.optParaNum, obj.intrinsicPara.lambda);
            for loop = 1 : opt.intrinsicPara.lambda
                xvec(:,loop) = obj.optPara.xmean + obj.intrinsicPara.sigma * (obj.intrinsicPara.B * obj.intrinsicPara.D * zvec(:,loop)); % add mutation
            end
            if iteration == 1 % for the 1st iteration: no noise added to the 1st particle
                xvec(:,1) = obj.optPara.xmean;
            end
        end

        function xxx = update(obj)
             % 继续写pseudo code
            % fit

        end

    end
end

