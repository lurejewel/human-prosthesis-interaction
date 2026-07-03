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
projName = 'human0714';

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
    a_opt           = ckp.a_opt;           % static-optimisation initial activations
    reflexParamMap  = ckp.reflexParamMap;
    reflexTemplate  = ckp.reflexTemplate;
    initPara        = ckp.initPara;
    sigma           = ckp.sigma;
    lassoFile        = '';
    if ~isempty(ckp.lassoFile)
        lassoFile = ckp.lassoFile;
    end

    % --- recompute model static properties (needed by init_infra) ---
    if isfield(ckp, 'dofNames')
        dofNames = ckp.dofNames;
    else
        % fallback for legacy checkpoints that predate dofNames storage
        dofNames = {'pelvis_tilt','pelvis_tx','pelvis_ty', ...
                    'hip_flexion_r','knee_extension_r','ankle_dorsiflexion_r', ...
                    'hip_flexion_l','knee_extension_l','ankle_dorsiflexion_l'};
    end
    modelStaticProp = read_muscle_static_prop(projName, simConfig, initPose, dofNames);

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
    s = prepare_fresh_start(projName);

    % --- extract setup fields ---
    optCfg              = s.optCfg;
    simConfig           = s.simConfig;
    initPose            = s.initPose;
    a_opt               = s.a_opt;
    modelStaticProp     = s.modelStaticProp;
    reflexParamMap      = s.reflexParamMap;
    reflexTemplate      = s.reflexTemplate;
    initPara            = s.initPara;
    sigma               = s.sigma;
    optConfig           = s.optConfig;
    nWorkers            = s.nWorkers;
    histLen             = s.histLen;
    bestFits            = s.bestFits;
    computationLoad     = s.computationLoad;
    g                   = s.g;
    gAtRestartStart     = s.gAtRestartStart;
    lastImprovementGen  = s.lastImprovementGen;
    bestFitEver         = s.bestFitEver;
    softBoostTriggered  = s.softBoostTriggered;
    restartCount        = s.restartCount;
    stop                = s.stop;
    lassoFile           = s.lassoFile;
    dofNames            = s.dofNames;

end  % ====================  END FRESH-START / RESUME BLOCK  ====================

% OUTER LOOP: for each generation
while g < optConfig.gMax && ~stop
    tic
    % generate parameters (particles) to be optimized
    [arx, arz] = optConfig.generate_parameters(g - gAtRestartStart);

    fits_cell = cell(nWorkers,1); % initialize fit in cell
    % MIDDLE LOOP: for each PARALLEL worker (batch of particles)
    for w = 1 : nWorkers

        fits_local = nan(optConfig.nParticles(w),1);
        [model, modelInfo] = init_infra(projName, modelStaticProp, a_opt);
        modelInfo.reflexParamMap = reflexParamMap;
        modelInfo.reflexTemplate = reflexTemplate;

        % INNER LOOP: for each particle
        for p = 1 : optConfig.nParticles(w) % optConfig.core.lambda
            k = sum(optConfig.nParticles(1:w)) - optConfig.nParticles(w) + p;
            modelInfo.reset_record();
            [state, modelInfo] = reset_particle_state(model, modelInfo, a_opt);  % per-particle deterministic reset
            modelInfo.read_muscleReflex_array(arx(:,k)); % read muscle reflex parameters

            % forward dynamic simulation and fitness evaluation
            modelInfo = forward_simulation(model, modelInfo, state);
            fits_local(p) = measure_simResults(modelInfo);
        end % END INNER LOOP: dynamic evaluation for at most (lambda/6) particles

        fits_cell{w} = fits_local;
    end % END MIDDLE LOOP [PARALLEL]: dynamic evaluation for 1 generation (lambda particles)

    % update elite particles (best3) and CMA-ES parameters
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
            simConfig, initPose, dofNames, lassoFile, a_opt);
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

end % END OUTER LOOP: CMA-ES optimization (g generations)


%% save optimization results

% ---- final checkpoint (overwrites periodic one) ----
save_checkpoint(checkpointFile, g, gAtRestartStart, lastImprovementGen, ...
    bestFitEver, softBoostTriggered, restartCount, stop, ...
    bestFits, computationLoad, optConfig, optCfg, projName, ...
    reflexParamMap, reflexTemplate, initPara, sigma, nWorkers, ...
    simConfig, initPose, lassoFile, a_opt);

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