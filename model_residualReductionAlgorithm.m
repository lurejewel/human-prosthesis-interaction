% -------------------------------------------------------------------------
% Name: model_resultsReductionAlgorithim.m
% Description: calculate model joint torques according to marker data
% (.trc), while realizing as little residual forces as possible by means of
% adjusting COM of some body of the model. Residual forces include
% pelvis_tilt/list/rotation, pelvis_tx/ty/tz.
% Info that needs to be specified:
% - in rra_setup.xml:
%   * model file name/path
%   * time range (initial time, final time)
%   * replaced forces: joint actuators, xxx_RRA_Actuators.xml
%   * output directory
%   * external loads file: xxx_grf.xml (same as model_inverseDynamics.m)
%   * coordinates file: output from the inverse kinematics tool (.mot file)
%   (same as model_inverseDynamics.m)
%   * task set file: xxx_RRA_Tasks.xml
%   * constraints file: xxx_RRA_ControlConstraints.xml
%   * lowpass cutoff frequency for coordinates: 6 Hz
%   * body COM to be adjusted: torso (maybe multiple bodies can be
%   selected?)
%   * 
%   * output force file name/path
% - in xxx_grf.xml:
%   * body to be applied: calcn_r/l
%   * data file
% Output:
% - xxx_id.sto: forces and torques for each DOF of the model
% -------------------------------------------------------------------------

addpath('assets\', 'model\')
import org.opensim.modeling.*

%% Run Inverse Kinematics Tool
model = Model('model\coupled_human-prosthesis_model_scaledFinal.osim');
state = model.initSystem();
idTool = InverseDynamicsTool('rra_setup.xml'); % configure ik tool
idTool.setModel(model); 
idTool.run();

%% Read .sto file and plot the force/torque profile
forceData = Storage('assets\subject01_walk1_id.sto');
frameNum = forceData.getSize; % number of frames
labels = forceData.getColumnLabels();
labelNum = labels.getSize; % number of dofs + 1 (1st col is time)

% x (time) and y (force/torque) for the plot
timeVec = nan(frameNum, 1);
forceVec = nan(frameNum, labelNum-1);
for t = 0 : frameNum - 1
    timeVec(t+1) = forceData.getStateVector(t).getTime;
    for dofIdx = 1 : labelNum - 1
        forceVec(t+1, dofIdx) = forceData.getStateVector(t).getData().get(dofIdx-1);
    end
end

% plot
for i = 1 : labelNum-1
    figure(i), plot(timeVec, forceVec(:,i));
    grid on, xlabel('time/s')
    ylabel(char(labels.get(i)),'Interpreter','none')
    title('inverse dynamics results: joint forces/torques')
end