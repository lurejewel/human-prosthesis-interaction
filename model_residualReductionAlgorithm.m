% -------------------------------------------------------------------------
% Name: model_resultsReductionAlgorithim.m
% Description: calculate model joint torques according to marker data
% (.trc), while realizing as little residual forces as possible by means of
% adjusting COM of some body of the model and joint kinematics. Residual
% forces include pelvis_tilt/list/rotation, pelvis_tx/ty/tz.
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
%   * output force file name/path
% - in xxx_grf.xml:
%   * body to be applied: calcn_r/l
%   * data file
% - in xxx_RRA_Tasks.xml:
%   * default linear controller, coordinate actuator, point actuator, point
%   actuator, torque actuator, Thelen2003Muscle, CMC joint
%   * CMC joint, exerted on each coordinate
% - in xxx_RRA_Actuators.xml:
%   * point actuator: FX, FY, FZ, exerted on pelvis
%   * torque actuator: MX, MY, MZ, exerted on pelvis
%   * coordinate actuator, exerted on biologicla joints
% - in xxx_RRA_ControlConstraints.xml:
%   * linear controller, exerted on actuators
% Output: 
% - xxx_rra_adjusted.osim: adjusted musculoskeletal model
% - xxx_RRA_Actuation_force.sto: forces or torques of the defined actuators
% - xxx_RRA_Actuation_speed.sto: linear or angular speeds of the defined
% acutators
% - xxx_RRA_Actuation_power.sto: power of the defined actuators
% (forces/torques multiply speeds)
% - xxx_RRA_controls.sto/xml: excitations of the defined actuators
% - xxx_RRA_Kinematics_q.sto: values of generalized coordiantes
% - xxx_RRA_Kinematics_u.sto: speeds of generalized coordiantes
% - xxx_RRA_Kinematics_dudt.sto: accelerations of generalized coordinates
% - xxx_RRA_states.sto: states of the model during RRA, including
% coordinate values and speeds
% - xxx_RRA_pErr.sto: position error of generalized coordinates
% - xxx_RRA_avgResiduals.txt: average residuals of FX/Y/Z and MX/Y/Z
% Note:
% - Recommended mass modification other than torso is displayed in the
% message window of OpenSim GUI, but not provided here.
% -------------------------------------------------------------------------

addpath('assets\', 'model\')
import org.opensim.modeling.*

%% run Residual Reduction Algorithm and display logger
loggerFile = fopen('opensim.log','rt');
fseek(loggerFile,0,'eof');

rraTool = RRATool('rra_setup.xml'); % configure rra tool
rraTool.run(); % RUN RRA

while ~feof(loggerFile)
    line = fgetl(loggerFile);
    disp(line)
end
fclose(loggerFile);

% report of average residual force
type('assets\ResultsRRA\subject01_walk1_RRA_avgResiduals.txt','r');

%% Read .sto file and plot the force/torque profile
forceData = Storage('assets\ResultsRRA\subject01_walk1_RRA_Actuation_force.sto');
frameNum = forceData.getSize; % number of frames
labels = forceData.getColumnLabels();
labelNum = labels.getSize; % number of dofs + 1 (1st column is time)

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
    figure(100+i), plot(timeVec, forceVec(:,i));
    grid on, xlabel('time/s')
    ylabel(char(labels.get(i)),'Interpreter','none')
    title('RRA results: joint forces/torques')
end