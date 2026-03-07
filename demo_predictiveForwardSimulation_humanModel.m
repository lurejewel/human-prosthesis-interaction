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
addpath(genpath('assets\'), genpath('model\'), genpath('functions\'))
% import org.opensim.modeling.*
projName = 'coupled_human-prosthesis_model';

%% activate parallel computing
nWorkers = 6; % 6 P-cores in Intel i7-12700H and Intel i5-12600KF
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
initPara = [1.27644955526603
    0.656012010875071
    -0.769176933585239
    1.53333505270711
    1.34378961328866
    -2.62406913538501
    -0.762577033023425
    0.0192735990058949
    -0.0428692717081615
    0.543656246235899
    0.132817428920201
    0.760791612116672
    0.762046844735927
    0.485059967268258
    -2.80667966041193
    0.732218582977587
    2.39622956431308
    0.608709134041632
    -0.245361663910839
    0.216300322985780
    1.98790010749173
    -0.0924298241125715
    0.322570530688795
    0.920656755428991
    1.82506209626037
    0.260820906580090
    0.446681199395220
    -0.290013940986199];
% initPara = load('bestPara_parallel.mat').bestPara;

sigma = 0.02; % deviation
optConfig = CMAES_optimization(initPara, sigma, nWorkers);
clear sigma

%% generate muscle excitation according to muscle reflex mechanism

% OUTER LOOP: for each generation
bestFits = nan(9999,1); 
computationLoad = nan(9999,1);
g = 0; % #generation
while g <= optConfig.gMax && ( g<=100 || (g>100 && g-find(diff(bestFits(1:g))~=0,1,'last') <= 100) ) % not reach the max / have diff in recent 100 generations
    tic
    % generate parameters (particles) to be optimized
    [arx, arz] = optConfig.generate_parameters(g);

    fits_cell = cell(nWorkers,1); % initialize fit in cell
    % MIDDLE LOOP: for each PARALLEL worker (batch of particles)
    for w = 1 : nWorkers
        
        fits_local = nan(optConfig.nParticles(w),1);
        [model, modelInfo, state] = init_infra(projName, modelStaticProp);

        % INNER LOOP: for each particle
        for p = 1 : optConfig.nParticles(w) % optConfig.core.lambda
            k = sum(optConfig.nParticles(1:w)) - optConfig.nParticles(w) + p;
            modelInfo.reset_record(); 
            modelInfo.read_muscleReflex_array(arx(:,k), arz(:,k)); % read muscle reflex parameters

            % forward dynamic simulation and fitness evaluation
            modelInfo = forward_simulation(model, modelInfo, state);
            fits_local(p) = measure_simResults(modelInfo);

        end % END INNER LOOP: dynamic evaluation for at most (lambda/6) particles
        fits_cell{w} = fits_local;
    end % END MIDDLE LOOP: dynamic evaluation for 1 generation (lambda particles)

    % update elite particles (best & big3) and CMA-ES parameters
    fits = vertcat(fits_cell{:});
    if isnan(sum(fits)) || min(fits)<6
        pause
    end
    optConfig.update_elite_fit(fits, arx, arz);
    optConfig.update_core(fits, g);

    g = g+1;
    if ~mod(g, 10)
        disp(['[', char(datetime), '] ', num2str(g), ' generations finished.']);
    end

    bestFits(g) = optConfig.recordForBestParticle.fit;
    computationLoad(g) = toc;

end % END OUTER LOOP: CMA-ES optimization (g generations)
% toc

% plot_simulation_results(model, modelInfo, 'excitation&activation')
% plot_simulation_results(model, modelInfo, 'muscleForce')
% plot_simulation_results(model, modelInfo, 'phase')

% bestPara = optConfig.recordForBestParticle.arx;
% save('bestPara.mat', "bestPara");
% save the states to .sto file
% save_state_as_sto(1, model, modelInfo);