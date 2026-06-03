% Name: determine_gait_phase_from_motion
% Author(s): Jin Wei, Peking U.
% Description: The function of this script is validate the availability of
% detecting the gait phase of each leg of the musculoskeletal model at real
% time based on the state of the model.
% There are five gait phases for each leg:
% - [Early Stance]: starting from the heel strike event, and ending at the
% time when the COM of pelvis surpasses the COM of foot (calcn).
% - [Late Stance]: ending at the oppsite foot clearance, OR the time when
% the COM of foot is left far behind the COM of pelvis (> 1 meter).
% - [Liftoff]: ending at foot clearance.
% - [Swing]: ending at the time when the COM of foot surpasses the COM of
% pelvis.
% - [Landing]: ending at the next heel strike.

addpath(genpath('..\assets\'), genpath('..\model\'), genpath('..\functions\'))
import org.opensim.modeling.*

%% load model
model = Model('model\coupled_human-prosthesis_model_scaledFinal.osim');
model.setUseVisualizer(true);
state = model.initSystem();

%% load motion and GRFdata
motData = Storage('assets\subject01_walk1_ik.mot');
rate = motData.getSize / (motData.getLastTime - motData.getFirstTime);
phase = nan(motData.getSize, 2);

motDataGRF = Storage('assets\subject01_walk1_grf.mot');

%% determine gait phase
for t = 0 : motData.getSize - 1 % for each frame/state
    stateVec = motData.getStateVector(t).getData; % get state vector for the model at the frame
    currentState = state.updQ(); % fetch the data array of the model at *t* state

    grfR = motDataGRF.getStateVector(t*10).getData.get(1); % t*10: freq of grf data is 10 times higher than that of ik motion data
    grfL = motDataGRF.getStateVector(t*10).getData.get(7);

    % assign the state vector to the model
    for i = 0 : stateVec.size - 1 % for each coordinate (DOF)
        coordVal = stateVec.get(i);
        if i<3 || i>5 % i=3~5 is translational dof, no need for conversion
            coordVal = deg2rad(coordVal);
        end
        currentState.set(i, coordVal);
    end
    model.realizeDynamics(state);
    % calculate gait phase based on the state of the model
    [phase(t+1,1), phase(t+1,2)] = cal_gait_phase(model, state, [grfR, grfL]); % phaseR phaseL
end

% %% debug: plot the gait phase
% timeVec = 0 : 1/rate : 2.5-1/rate;
% figure, plot(timeVec, phase(:,1))
% hold on, plot(timeVec, phase(:,2))
% xlabel('t/s'); ylabel('gait phase'); legend('right leg', 'left leg')

%% play the video
for t = 0 : motData.getSize - 1 % for each frame/state
    stateVec = motData.getStateVector(t).getData; % get state vector for the model at the frame
    currentState = state.updQ(); % fetch the data array of the model at *t* state
    % for each DOF, assign the state from motion file (stateVec) to the
    % state of the model (currentState)
    for i = 0 : stateVec.size - 1
        coordVal = stateVec.get(i); 
        if i<3 || i>5 % i=3~5 is translational dof, no need for conversion
            coordVal = deg2rad(coordVal);
        end
        currentState.set(i, coordVal);
        % model.realizeDynamics(state); % not necessary, but better do it
    end
    model.updVisualizer().show(state);
    pause(1/rate)
end
