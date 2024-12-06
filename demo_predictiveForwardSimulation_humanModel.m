% -------------------------------------------------------------------------
% Name: demo_predictiveForwardSimulation_humanModel.m
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
% - What if more muscles are included?
% - What if current model is replaced with human-prosthesis model?
% - What if we use alternative optimization methods other than CMA-ES?
% - What if hundreds of muscles should be included?
% - What if other periodic activities (running, swimming, etc.) are
% considered?
% -------------------------------------------------------------------------

clear all; close all; clc
addpath('assets\', 'model\', 'functions\')
import org.opensim.modeling.*

%% load model, motion and GRF data
model = Model('model\coupled_human-prosthesis_model.osim');
model.setUseVisualizer(0);
model = add_muscle_actuator(model, PrescribedController()); % add the central controller and muscle actuators

% initial muscle reflex parameters; not optimized yet
muscleReflexParaArray = [1.14339427675597;0.683965164606834;-0.654937183750042;1.44526338204981;1               ;-2.54053931319423;-0.666188349877246;0.0496818438178654;0.109481378880035;0.388522196081861;0.266059683859970;0.748599386315988;0.543366924759648;0.598978415045811;-2.81541067233053;0.377854695673883;2.18759273326043;0.657016639688066;-0.298290360392486;0.300063335235475;1.81474706127338;0.00878201304411728;0.510662826532988;1.15502124130598;1.97028072512268;0.479041225785079;0.568851978040139;-0.188673608490112];

% older version ---
% motData = Storage('assets\subject01_walk1_ik.mot');
% firstTime = motData.getFirstTime; % t_start = 0
% lastTime = motData.getLastTime; % t_end = 2.5
% frameNum = motData.getSize;
% rate = frameNum / (lastTime - firstTime);
% coordNum = motData.getStateVector(0).getData.size; % number of coordinates (DOFs)
% motDataGRF = Storage('assets\subject01_walk1_grf.mot'); % GRF rate is 10x than ik rate

%% specify simulation configuration
simConfig.startTime = 0;
simConfig.endTime = 1;
simConfig.stepTime = 0.001;
simConfig.saveSTO = 0; % 1 = true, 0 = false

%% init modelInfo Class
modelInfo = ModelInfo(model, simConfig, muscleReflexParaArray);

%% first-frame initialization
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
    if frameIndex <= modelInfo.muscleReflexDelay % the muscle excitations in the first 10 frames = those in the 11th frames (delay = 10 ms)
        modelInfo.muscleExcitations(:,frameIndex) = modelInfo.muscleExcitations(:,11);
    end
    brain = model.updControllerSet.get(0);
    brain = PrescribedController.safeDownCast(brain);
    for muscleIndex = 0 : model.getMuscles.getSize - 1
        brain.prescribeControlForActuator(muscleIndex, Constant(modelInfo.muscleExcitations(muscleIndex+1,frameIndex)));
    end

    % forward simulation
    manager = Manager(model);
    manager.initialize(state);
    finalState = manager.integrate(t);
    finalState.setTime(t);

    % update & record
    state = finalState;
    model.realizeDynamics(state);
    modelInfo.stateHistory(:,frameIndex) = vec_2_mat(state.getY);
    frameIndex = frameIndex + 1;
    if frameIndex <= width(modelInfo.time)
        modelInfo = update_modelInfo(model, modelInfo, state, frameIndex);
    end

end

clear t frameIndex muscleIndex brain finalState state
toc

% plot_simulation_results(model, modelInfo, 'excitation&activation')
% plot_simulation_results(model, modelInfo, 'muscleForce')
% plot_simulation_results(model, modelInfo, 'phase')

% [after CMA-ES optimization (or other methods), muscleReflexParaArray is changed...]
% modelInfo.muscleReflex = muscleReflexPara_array2struct(muscleReflexParaArray); % muscle reflex parameters assignment

%% save the states to .sto file

if simConfig.saveSTO

    labels = ArrayStr();
    labels.append('time');
    numTags = height(modelInfo.stateHistory);
    for stateIndex = 0 : numTags-1
        labels.append(model.getStateVariableNames.get(stateIndex));
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

end

clear Y labels stateIndex stateVector frameIndex numTags stoFile