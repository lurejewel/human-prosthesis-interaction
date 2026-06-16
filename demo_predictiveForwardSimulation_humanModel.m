% -------------------------------------------------------------------------
% Name: demo_predictiveForwardSimulation_humanModel.m
% Author(s): Jin Wei, Peking U. wjin24@stu.pku.edu.cn
% Description: Functional validation script of the predictive,
% experimental-data-free, muscle-reflex-based musculoskeletal model.
% Core components:
% - [Model]: 16 DOF, 14 muscles (only lower limbs)
% - [Gait phase detection]: calculate gait phases of each leg based on the
% state of the model at real time
% - [Muscle reflex mechanism]: CNS-like part where muscle excitations are
% generated according to gait phases at real time
% - [Optimization]: tune the muscle-reflex parameters in order that the
% model can walk in a smooth, human-like way
% What if...:
% - What if more muscles are included? -> What if hundreds of muscles are
% included?
% - What if current model is replaced with human-prosthesis model?
% - What if we use alternative optimization methods other than CMA-ES?
% - What if other periodic activities (running, swimming, etc.) are
% considered?
% -------------------------------------------------------------------------
% 目前的问题：
% 1. 增加独立的处理模块：导入\导出数据；显示trc数据；显示mot数据，etc.
% 2. 反射参数换成Map映射的赋值形式
% 3. 目前是2D还是3D？如何改变地形（上下坡、上下台阶等）？
clear all; close all; clc
addpath(genpath('assets\'), genpath('model\'), genpath('functions\'), genpath('results\'))
% import org.opensim.modeling.*
projName = 'coupled_human-prosthesis_model';

%% =======================================================================
%  CHECKPOINT CONFIGURATION
%  Set resumeFromCheckpoint = true to restart from the last saved snapshot.
%  The checkpoint file is overwritten every 10 generations.
% =======================================================================
resumeFromCheckpoint = false;
checkpointFile = 'results\checkpoint.mat';

if resumeFromCheckpoint && isfile(checkpointFile)
    % ====================  RESUME FROM CHECKPOINT  ====================
    fprintf('[%s] Resuming from checkpoint: %s\n', char(datetime), checkpointFile);
    ckp = load(checkpointFile);
    fprintf('[%s] Checkpoint from gen %d (bestFit=%.6g, timestamp=%s).\n', ...
        char(datetime), ckp.g, ckp.optConfig.recordForBestParticle.fit, ...
        char(ckp.timestamp));

    % --- restore loop bookkeeping ---
    g                   = ckp.g;
    gAtRestartStart     = ckp.gAtRestartStart;
    lastImprovementGen  = ckp.lastImprovementGen;
    bestFitEver         = ckp.bestFitEver;
    softBoostTriggered  = ckp.softBoostTriggered;
    restartCount        = ckp.restartCount;
    stop                = ckp.stop;

    % --- restore CMA-ES object (full internal state) ---
    optConfig = ckp.optConfig;

    % --- restore immutable configuration ---
    optCfg          = ckp.optCfg;
    projName        = ckp.projName;
    simConfig       = ckp.simConfig;
    initPose        = ckp.initPose;
    reflexParamMap  = ckp.reflexParamMap;
    reflexTemplate  = ckp.reflexTemplate;
    initPara        = ckp.initPara;
    sigma           = ckp.sigma;
    lassoFile        = '';
    if ~isempty(ckp.lassoFile)
        lassoFile = ckp.lassoFile;
    end

    % --- recompute model static properties (needed by init_infra) ---
    modelStaticProp = read_muscle_static_prop(projName, simConfig, initPose);

    % --- restore & re-pad history arrays ---
    bestFits        = ckp.bestFits(:);
    computationLoad = ckp.computationLoad(:);
    histLen = (optCfg.maxRestarts + 1) * (optCfg.patience + 50) + 200;
    bestFits(end+1:histLen)        = nan;
    computationLoad(end+1:histLen) = nan;

    % --- activate parallel pool on current machine ---
    nWorkers = feature('numcores');
    p = gcp('nocreate');
    if isempty(p) || p.NumWorkers ~= nWorkers
        delete(p);
        parpool('local', nWorkers);
    end

    % --- adjust CMA-ES worker distribution if core count changed ---
    if nWorkers ~= optConfig.nWorkers
        fprintf('[%s] Worker count changed: %d -> %d. Redistributing particles.\n', ...
            char(datetime), optConfig.nWorkers, nWorkers);
        optConfig.nWorkers = nWorkers;
        optConfig.nParticles = floor(optConfig.core.lambda / nWorkers) * ones(nWorkers, 1);
        extra = mod(optConfig.core.lambda, nWorkers);
        if extra > 0
            optConfig.nParticles(1:extra) = optConfig.nParticles(1:extra) + 1;
        end
    end

    % --- reseed random generator for reproducibility ---
    rng(optCfg.rngSeed + g, 'threefry');

else
    % ====================  FRESH START  ====================
    resumeFromCheckpoint = false;

    %% reproducibility configuration
    optCfg.rngSeed         = 2026;
    optCfg.softPatience    = 30;     % stall gens before a soft sigma-boost
    optCfg.patience        = 100;    % stall gens before an IPOP restart
    optCfg.maxRestarts     = 3;      % IPOP-CMA-ES: hard stop after this many restarts
    optCfg.sigmaBoostFactor= 2.0;    % multiplicative kick on soft stall
    optCfg.minImprovement  = 1e-12;  % absolute fitness improvement threshold to reset stall
    rng(optCfg.rngSeed, 'threefry');


    %% activate parallel computing
    nWorkers = feature('numcores');
    p = gcp('nocreate');
    if isempty(p) || p.NumWorkers ~= nWorkers
        delete(p);
        parpool('local', nWorkers);
    end

    %% simulation and model configuration
    simConfig.endTime = 10;
    simConfig.stepTime = 0.005;
    simConfig.saveSTO = 1; % 1 = true, 0 = false
    simConfig.speed = 1.0; % desired walking speed
    simConfig.slope = 0;
    initPose = [-0.105763, 0, 0.900237, 0.439316, 0.198813,-0.393922, -1.03755, 0.104714, -0.348473, ... % value of generalized coordinates
        -0.0895989, 1.0757, 0.1543, -1.35971, 3.34368, 0.267883, -3.15281, 0.840122, 1.26642]; % velocity of generalized coordinates
    modelStaticProp = read_muscle_static_prop(projName, simConfig, initPose);

    %% optimization configuration
    % initial muscle reflex parameters; not optimized yet
    % muscleReflexParaArray = [1.14339427675597;0.683965164606834;-0.654937183750042;1.44526338204981;1               ;-2.54053931319423;-0.666188349877246;0.0496818438178654;0.109481378880035;0.388522196081861;0.266059683859970;0.748599386315988;0.543366924759648;0.598978415045811;-2.81541067233053;0.377854695673883;2.18759273326043;0.657016639688066;-0.298290360392486;0.300063335235475;1.81474706127338;0.00878201304411728;0.510662826532988;1.15502124130598;1.97028072512268;0.479041225785079;0.568851978040139;-0.188673608490112];
    % initPara = [1.27644955526603
    %     0.656012010875071
    %     -0.769176933585239
    %     1.53333505270711
    %     1.34378961328866
    %     -2.62406913538501
    %     -0.762577033023425
    %     0.0192735990058949
    %     -0.0428692717081615
    %     0.543656246235899
    %     0.132817428920201
    %     0.760791612116672
    %     0.762046844735927
    %     0.485059967268258
    %     -2.80667966041193
    %     0.732218582977587
    %     2.39622956431308
    %     0.608709134041632
    %     -0.245361663910839
    %     0.216300322985780
    %     1.98790010749173
    %     -0.0924298241125715
    %     0.322570530688795
    %     0.920656755428991
    %     1.82506209626037
    %     0.260820906580090
    %     0.446681199395220
    %     -0.290013940986199];
    % initPara = load('results\bestPara_noparallel.mat').bestPara;
    % warmStartPack   = load('results\opt_result_2026-06-10_22-32-53.mat').result;
    % initPara        = warmStartPack.bestPara;
    % reflexParamMap  = warmStartPack.reflexParamMap;
    % reflexTemplate  = warmStartPack.reflexTemplate;

    %% ---- LASSO linear-phase controller initialisation ----
    lassoFile = fullfile('results', 'lasso_controller_result.mat');
    [initPara, reflexParamMap, reflexTemplate] = ...
        load_lasso_reflex_controller(lassoFile);
    % lassoFile = '';  % not used when warm-starting; set for checkpoint compatibility

    sigma = 0.02; % deviation
    optConfig = CMAES_optimization(initPara, sigma, nWorkers);

    %% generate muscle excitation according to muscle reflex mechanism

    % OUTER LOOP: for each generation
    % (history arrays are pre-allocated generously; with IPOP restarts the run
    % can be up to (maxRestarts+1) * patience long, so size accordingly)
    histLen          = (optCfg.maxRestarts + 1) * (optCfg.patience + 50) + 200;
    bestFits         = nan(histLen,1);
    computationLoad  = nan(histLen,1);

    g                = 0;          % current generation index
    gAtRestartStart  = 0;           % g at the most recent restart
    lastImprovementGen = 0;
    bestFitEver      = inf;
    softBoostTriggered = false;     % did we already kick sigma during this stall?
    restartCount     = 0;

    stop = false;

end  % ====================  END FRESH-START / RESUME BLOCK  ====================

while g < optConfig.gMax && ~stop
    tic % , profile on
    % generate parameters (particles) to be optimized
    [arx, arz] = optConfig.generate_parameters(g - gAtRestartStart);

    fits_cell = cell(nWorkers,1); % initialize fit in cell
    % MIDDLE LOOP: for each PARALLEL worker (batch of particles)
    for w = 1 : nWorkers

        fits_local = nan(optConfig.nParticles(w),1);
        [model, modelInfo] = init_infra(projName, modelStaticProp);
        modelInfo.reflexParamMap = reflexParamMap;
        modelInfo.reflexTemplate = reflexTemplate;

        % INNER LOOP: for each particle
        for p = 1 : optConfig.nParticles(w) % optConfig.core.lambda
            k = sum(optConfig.nParticles(1:w)) - optConfig.nParticles(w) + p;
            % Create a local copy of modelInfo for this parfor iteration.
            % This avoids the "temporary variable" classification error.
            % mi = modelInfo;
            modelInfo.reset_record();
            [state, modelInfo] = reset_particle_state(model, modelInfo);  % per-particle deterministic reset
            modelInfo.read_muscleReflex_array(arx(:,k)); % read muscle reflex parameters

            % forward dynamic simulation and fitness evaluation
            modelInfo = forward_simulation(model, modelInfo, state);
            fits_local(p) = measure_simResults(modelInfo);

        end % END INNER LOOP: dynamic evaluation for at most (lambda/6) particles
        fits_cell{w} = fits_local;
    end % END MIDDLE LOOP: dynamic evaluation for 1 generation (lambda particles)

    % update elite particles (best & big3) and CMA-ES parameters
    fits = vertcat(fits_cell{:});
    optConfig.update_elite_fit(fits, arx, arz);
    optConfig.update_core(fits, g - gAtRestartStart);

    g = g+1;
    if ~mod(g, 10)
        disp(['[', char(datetime), '] ', num2str(g), ' generations finished.']);
    end

    bestFits(g) = optConfig.recordForBestParticle.fit;
    computationLoad(g) = toc;

    % --- periodic checkpoint (crash recovery) ---
    if ~mod(g, 10) || stop
        save_checkpoint(checkpointFile, g, gAtRestartStart, lastImprovementGen, ...
            bestFitEver, softBoostTriggered, restartCount, stop, ...
            bestFits, computationLoad, optConfig, optCfg, projName, ...
            reflexParamMap, reflexTemplate, initPara, sigma, nWorkers, ...
            simConfig, initPose, lassoFile);
    end

    % --- track improvement -------------------------------------------------
    if optConfig.recordForBestParticle.fit < bestFitEver - optCfg.minImprovement
        bestFitEver        = optConfig.recordForBestParticle.fit;
        lastImprovementGen = g;
        softBoostTriggered = false;     % a fresh improvement -> arm soft boost again
    end
    stallGens = g - lastImprovementGen;

    % --- staged stall handling : (a) soft sigma-boost, (b) IPOP restart, (c) stop ---
    if stallGens >= optCfg.softPatience && ~softBoostTriggered
        optConfig.sigma_boost(optCfg.sigmaBoostFactor);
        softBoostTriggered = true;
    end

    if stallGens >= optCfg.patience
        if restartCount < optCfg.maxRestarts
            optConfig.restart();          % IPOP-CMA-ES : double lambda, reset state
            restartCount       = restartCount + 1;
            gAtRestartStart    = g;       % CMA-ES counter restarts at 0 internally
            lastImprovementGen = g;       % give the new run a fair chance
            softBoostTriggered = false;
        else
            fprintf('[CMA-ES] Stall after %d gens with %d restarts -> stop.\n', ...
                stallGens, restartCount);
            stop = true;
        end
    end
    % profile viewer

end % END OUTER LOOP: CMA-ES optimization (g generations)


%% save optimization results

% ---- final checkpoint (overwrites periodic one) ----
save_checkpoint(checkpointFile, g, gAtRestartStart, lastImprovementGen, ...
    bestFitEver, softBoostTriggered, restartCount, stop, ...
    bestFits, computationLoad, optConfig, optCfg, projName, ...
    reflexParamMap, reflexTemplate, initPara, sigma, nWorkers, ...
    simConfig, initPose, lassoFile);

% ---- assemble result struct (for analysis, lighter than checkpoint) ----
bestFits = bestFits(1:g);              % trim to actual generations
computationLoad = computationLoad(1:g);

result.bestFit       = optConfig.recordForBestParticle.fit;
result.bestPara      = optConfig.recordForBestParticle.arx;
result.generations   = g;
result.fitHistory    = bestFits;
result.timeHistory   = computationLoad;
result.rngSeed       = optCfg.rngSeed;
result.patience      = optCfg.patience;
result.softPatience  = optCfg.softPatience;
result.maxRestarts   = optCfg.maxRestarts;
result.restartCount  = restartCount;
result.finalLambda   = optConfig.core.lambda;
result.finalSigma    = optConfig.core.sigma;
result.timestamp     = datetime;

% ---- controller metadata (essential for reconstructing bestPara) ----
result.controllerType    = 'lasso_linear_phase';
if exist('lassoFile', 'var') && ~isempty(lassoFile)
    result.lassoFile = lassoFile;
else
    result.lassoFile = '';
end
result.reflexParamMap    = reflexParamMap;
result.reflexTemplate    = reflexTemplate;
result.bestReflexParams  = unpack_lasso_reflex_params( ...
    result.bestPara, reflexParamMap, reflexTemplate);


outFile = ['results\opt_result_', datestr(now, 'yyyy-mm-dd_HH-MM-SS'), '.mat'];
save(outFile, 'result');
fprintf('[%s] Optimization finished. %d generations, best fit = %.6g. Saved to %s\n', ...
    char(datetime), g, result.bestFit, outFile);

% plot_simulation_results(model, modelInfo, 'excitation&activation')
% plot_simulation_results(model, modelInfo, 'muscleForce')
% plot_simulation_results(model, modelInfo, 'phase')

% bestPara = optConfig.recordForBestParticle.arx;
% save('bestPara.mat', "bestPara");
% save the states to .sto file
% save_state_as_sto(1, model, modelInfo);