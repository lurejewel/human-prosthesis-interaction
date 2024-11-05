% -------------------------------------------------------------------------
% Name: model_staticOptimization.m
% Description: 
% Info that needs to be specified:
% - in so_setup.xml
%   * model file name/path, output model after RRA
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
% -------------------------------------------------------------------------

addpath('assets\', 'model\')
import org.opensim.modeling.*

%% run Static Optimization and display logger
loggerFile = fopen('opensim.log','rt');
fseek(loggerFile,0,'eof');

% load model
model = Model('model\subject01_rra_adjusted.osim');
state = model.initSystem();

% configure static optimization tool
% soTool = StaticOptimization(); % configure so tool\
soTool = AnalyzeTool('so_setup.xml');
% soTool.setModel(model);
% soTool.set
soTool.run(); % RUN SO

while ~feof(loggerFile)
    line = fgetl(loggerFile);
    disp(line)
end
fclose(loggerFile);

%% Read .sto file and plot the force/torque profile
% forceData = Storage('assets\ResultsRRA\subject01_walk1_RRA_Actuation_force.sto');
% frameNum = forceData.getSize; % number of frames
% labels = forceData.getColumnLabels();
% labelNum = labels.getSize; % number of dofs + 1 (1st column is time)
% 
% % x (time) and y (force/torque) for the plot
% timeVec = nan(frameNum, 1);
% forceVec = nan(frameNum, labelNum-1);
% for t = 0 : frameNum - 1
%     timeVec(t+1) = forceData.getStateVector(t).getTime;
%     for dofIdx = 1 : labelNum - 1
%         forceVec(t+1, dofIdx) = forceData.getStateVector(t).getData().get(dofIdx-1);
%     end
% end
% 
% % plot
% for i = 1 : labelNum-1
%     figure(100+i), plot(timeVec, forceVec(:,i));
%     grid on, xlabel('time/s')
%     ylabel(char(labels.get(i)),'Interpreter','none')
%     title('SO results: joint forces/torques')
% end