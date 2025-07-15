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
% 1. 膝关节过伸比较严重，是否能够施加强约束？
% 2. CMAES参数更新代码尚未完成；CMAES部分未验证
% 3. 增加独立的处理模块：导入\导出数据；显示trc数据；显示mot数据，etc.
% 4. 反射参数换成Map映射的赋值形式
% 5. 使用OpenSim API可能同样存在内存泄露的问题
clear all; close all; clc
addpath('assets\', 'model\', 'functions\')
import org.opensim.modeling.*

%% load model
model = Model('model\coupled_human-prosthesis_model.osim');
model.setUseVisualizer(0); % display animation in real-time. set to 1(true) to activate, or 0(false) to deactivate.
model = add_muscle_actuator(model, PrescribedController()); % add the central controller and muscle actuators

%% specify simulation configuration
simConfig.endTime = 10;
simConfig.stepTime = 0.002;
simConfig.saveSTO = 1; % 1 = true, 0 = false
simConfig.speed = 1.0; % desired walking speed

%% specify optimization configuration
% initial muscle reflex parameters; not optimized yet
muscleReflexParaArray = [1.14339427675597;0.683965164606834;-0.654937183750042;1.44526338204981;1               ;-2.54053931319423;-0.666188349877246;0.0496818438178654;0.109481378880035;0.388522196081861;0.266059683859970;0.748599386315988;0.543366924759648;0.598978415045811;-2.81541067233053;0.377854695673883;2.18759273326043;0.657016639688066;-0.298290360392486;0.300063335235475;1.81474706127338;0.00878201304411728;0.510662826532988;1.15502124130598;1.97028072512268;0.479041225785079;0.568851978040139;-0.188673608490112];
sigma = 0.02; % deviation
optConfig = CMAES_optimization(muscleReflexParaArray, sigma);
clear sigma

%% first-frame initialization
modelInfo = ModelInfo(model, simConfig, muscleReflexParaArray);
% this is ugly, i know :( may find better way to do this later
[model, modelInfo, state] = init_first_frame(model, modelInfo); % including model.initSystem()
% model.updVisualizer().show(modelInfo.state);

%% generate muscle excitation according to muscle reflex mechanism

% OUTER LOOP: for each generation
bestFits = [];
while modelInfo.g <= modelInfo.gMax || (modelInfo.g > 3 && modelInfo.g - find(diff(bestFits)~=0,1,'last') < 500) % not reach the max / have diff in recent 500 generations

    % generate parameters (particles) to be optimized
    [arx, arz] = optConfig.generate_parameters(modelInfo.g);

    % initialize fits
    fits = nan(1, optConfig.core.lambda);
    

    % INNER LOOP: for each particle
    for k = 1 : optConfig.core.lambda

        % read muscle reflex parameters
        modelInfo.read_muscleReflex_array(arx(:,k));
        modelInfo.arz = arz(:,k);
        modelInfo.reset_record();

        % forward dynamic simulation of walking for the neuro-musculoskeletal model
        modelInfo = forward_simulation(model, modelInfo, state);

        % evaluation for the simulational results
        fits(k) = measure_simResults(optConfig, simConfig, modelInfo);

        % stage 1 -> stage 2
        if modelInfo.stage == 1 && (fits(k) < 0 || modelInfo.lastTime >= simConfig.endTime)
            disp(['stage 1 optimization completed, fit: ' num2str(fits(k))]);
            modelInfo.stage = 2;
            fits(k) = measure_simResults(t, optConfig, simConfig, modelInfo);
            disp(['stage 2 fit: ' num2str(fits(k))]);
            pause

            % 这里要多写一些，因为变为二阶段时，一阶段的相关精英fit记录就要清除了
            % 把当前的粒子设为best
            % 以及要不要continue
        end

        % recognize & update elite particles (best & big3)
        optConfig.update_elite_fit(fits(k), modelInfo)

    end % END LOOP: dynamic evaluation for 1 generation (lambda particles)

    % update CMA-ES parameters
    optConfig.update_core(fits, modelInfo.g);

    modelInfo.g = modelInfo.g + 1;
    bestFits = [bestFits, optConfig.recordForBestParticle.fit];
    
    if ~mod(modelInfo.g, 10)
        disp(['[', char(datetime), '] ', num2str(modelInfo.g), ' generations finished.']);
    end

end % END LOOP: CMA-ES optimization (g generations)
% toc

% plot_simulation_results(model, modelInfo, 'excitation&activation')
% plot_simulation_results(model, modelInfo, 'muscleForce')
% plot_simulation_results(model, modelInfo, 'phase')

% [after CMA-ES optimization (or other methods), muscleReflexParaArray is changed...]
% QUESTION：是否可以把均值到均值的步长，改为最优到最优的步长？
% 考虑到最优到最优步长很可能为零，会出现什么后果？是否有解决方法？
% modelInfo.muscleReflex = muscleReflexPara_array2struct(muscleReflexParaArray); % muscle reflex parameters assignment


%% save the states to .sto file

if simConfig.saveSTO

    state = model.initSystem();
    labels = ArrayStr();
    labels.append('time');
    for stateIndex = 0 : state.getNQ-1
        labels.append([char(model.getCoordinateSet.get(stateIndex)) '/value']);
    end
    for stateIndex = 0 : state.getNU-1
        labels.append([char(model.getCoordinateSet.get(stateIndex)) '/speed']);
    end
    for stateIndex = state.getNQ + state.getNU : state.getNQ + state.getNU + state.getNZ-1
        labels.append([char(model.getStateVariableNames.get(stateIndex))]); % NOTE that the order of state in Model class is not the same as that in State class
    end

    stoFile = Storage();
    stoFile.setName(['human0914_sim_' char(datetime("today")) '.sto']);
    stoFile.setColumnLabels(labels);

    for frameIndex = 1 : length(modelInfo.time)
        Y = mat_2_vec(modelInfo.stateHistory(:,frameIndex));
        stateVector = StateVector(modelInfo.time(frameIndex), Y);
        stoFile.append(stateVector);
    end

    stoFile.print(['human0914_sim_' char(datetime("today")) '.sto']);
    read_sto_file('model/coupled_human-prosthesis_model.osim', ['human0914_sim_' char(datetime("today")) '.sto'], 1);
end


clear Y labels stateIndex stateVector frameIndex numTags stoFile