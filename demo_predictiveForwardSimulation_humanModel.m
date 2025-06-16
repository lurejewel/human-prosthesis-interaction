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
clear all; close all; clc
addpath('assets\', 'model\', 'functions\')
import org.opensim.modeling.*

%% load model
model = Model('model\coupled_human-prosthesis_model.osim');
model.setUseVisualizer(0); % display animation in real-time. set to 1(true) to activate, or 0(false) to deactivate.
model = add_muscle_actuator(model, PrescribedController()); % add the central controller and muscle actuators

%% specify simulation configuration
simConfig.startTime = 0;
simConfig.endTime = 10;
simConfig.stepTime = 0.002;
simConfig.saveSTO = 1; % 1 = true, 0 = false
simConfig.speed = 1.0; % desired walking speed

%% specify optimization configuration
% initial muscle reflex parameters; not optimized yet
muscleReflexParaArray = [1.14339427675597;0.683965164606834;-0.654937183750042;1.44526338204981;1               ;-2.54053931319423;-0.666188349877246;0.0496818438178654;0.109481378880035;0.388522196081861;0.266059683859970;0.748599386315988;0.543366924759648;0.598978415045811;-2.81541067233053;0.377854695673883;2.18759273326043;0.657016639688066;-0.298290360392486;0.300063335235475;1.81474706127338;0.00878201304411728;0.510662826532988;1.15502124130598;1.97028072512268;0.479041225785079;0.568851978040139;-0.188673608490112];
sigma = 0.02; % deviation
stage = 1; % 1 or 2, for the 2-stage CMA-ES optimization (1 by default)
optConfig = CMAES_optimization(muscleReflexParaArray, sigma, stage);
clear sigma stage speed

%% first-frame initialization
modelInfo = ModelInfo(model, simConfig, muscleReflexParaArray);
% this is ugly, i know :( may find better way to do this later
[model, modelInfo, state] = init_first_frame(model, modelInfo); % including model.initSystem()
% model.updVisualizer().show(modelInfo.state);

%% generate muscle excitation according to muscle reflex mechanism

tic
frameIndex = 1;
for t = modelInfo.time % time series

    % calculate ground reaction forces, including normal forces and
    % friction forces, exerted by ground on each foot (calcn_r/l)
    modelInfo = cal_grf(model, modelInfo, state, frameIndex); % [QUESTION: IS THIS UGLY AND INEFFICIENT???]

    % calculate gait phase based on the state of the model
    modelInfo = cal_gait_phase(model, modelInfo, state, frameIndex);

    % calculate muscle excitations based on muscle states and gait phases
    modelInfo.muscleExcitations(:,frameIndex + modelInfo.muscleReflexDelay) = cal_muscle_excitation(modelInfo, state, frameIndex);
    if frameIndex <= modelInfo.muscleReflexDelay % the muscle excitations in the first few frames = the frames right after the delay (10 frames if delay=10ms)
        modelInfo.muscleExcitations(:,frameIndex) = modelInfo.muscleExcitations(:,modelInfo.muscleReflexDelay+1);
    end

    % assign muscle excitations (control signals) through brain (controller) to muscles (actuators)
    brain = model.updControllerSet.get(0);
    brain = PrescribedController.safeDownCast(brain);
    for muscleIndex = 0 : model.getMuscles.getSize - 1
        brain.prescribeControlForActuator(muscleIndex, Constant(modelInfo.muscleExcitations(muscleIndex+1,frameIndex)));
    end

    % forward simulation 这里的逻辑对吗？需要单独检查一下
    manager = Manager(model);
    manager.initialize(state);
    finalState = manager.integrate(t); % 这里是从哪个初始状态（是否会因为model.initialize而清零），向前推进了多少时间？
    finalState.setTime(t); 

    % update & record
    state = finalState;
    model.realizeDynamics(state);
    % modelInfo.stateHistory(:,frameIndex) = vec_2_mat(state.getY);
    frameIndex = frameIndex + 1;
    if frameIndex <= width(modelInfo.time)
        modelInfo = update_modelInfo(model, modelInfo, state, frameIndex);
    end

    if model.getCoordinateSet.get('pelvis_ty').getValue(state) < 0.6
        disp(['Model fell down. Maximum distance reached: ' num2str(model.getCoordinateSet.get('pelvis_tx').getValue(state)) ' m.']);
        break;
    end

end

fit = measure_simResults(t, optConfig, simConfig, modelInfo);

% stage 1 -> stage 2
if optConfig.hyperPara.stage == 1 && (fit < -simConfig.endTime*simConfig.speed || t == simConfig.endTime)
    disp(['stage 1 optimization completed, fit: ' num2str(fit)]);
    optConfig.hyperPara.stage = 2;
    fit = measure_simResults(t, optConfig, simConfig, modelInfo);
    disp(['stage 2 fit: ' num2str(fit)]);
end

clear t frameIndex muscleIndex brain finalState state
toc

% plot_simulation_results(model, modelInfo, 'excitation&activation')
% plot_simulation_results(model, modelInfo, 'muscleForce')
% plot_simulation_results(model, modelInfo, 'phase')

% [after CMA-ES optimization (or other methods), muscleReflexParaArray is changed...]
% 对于CMAES这个方法，最后别忘了optConfig.optPara = mean(这lambda个参数值）
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