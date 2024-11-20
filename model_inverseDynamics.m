% -------------------------------------------------------------------------
% Name: model_inverseDynamics.m
% Description: calculate model joint torques according to marker data (.trc)
% Info that needs to be specified:
% - in id_setup.xml:
%   * model file name/path
%   * time range
%   * forces to exclude: muscle
%   * external loads file: xxx_grf.xml
%   * coordinates file: output from the inverse kinematics tool (.mot file)
%   * lowpass cutoff frequency for coordinates: 6 Hz
%   * output force file name/path
% - in xxx_grf.xml:
%   * body to be applied: calcn_r/l
%   * data file
% Output:
% - xxx_id.sto: forces and torques for each DOF of the model
% -------------------------------------------------------------------------

addpath('assets\', 'model\')
import org.opensim.modeling.*

%% run Inverse Kinematics Tool and display logger
loggerFile = fopen('opensim.log','rt');
fseek(loggerFile,0,'eof');

model = Model('model\coupled_human-prosthesis_model_scaledFinal.osim');
state = model.initSystem();
idTool = InverseDynamicsTool('id_setup.xml'); % configure ik tool
idTool.setModel(model);  
idTool.run(); % RUN INVERSE DYNAMICS

while ~feof(loggerFile)
    line = fgetl(loggerFile);
    disp(line)
end
fclose(loggerFile);

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
    figure(i), plot(timeVec, forceVec(:,i), 'LineWidth',2);
    grid on, xlabel('time/s')
    ylabel(char(labels.get(i)),'Interpreter','none')
    title('inverse dynamics results: joint forces/torques')
end

% % --- debug --- %
% forceData = Storage('C:\Users\Taryn Wong\Documents\OpenSim\4.5\Models\Gait2354_Simbody\ResultsInverseDynamics\inverse_dynamics.sto');
% frameNum = forceData.getSize; % number of frames
% labels = forceData.getColumnLabels();
% labelNum = labels.getSize; % number of dofs + 1 (1st col is time)
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
%     figure(i+100), plot(timeVec, forceVec(:,i));
%     grid on, xlabel('time/s')
%     ylabel(char(labels.get(i)),'Interpreter','none')
%     title('inverse dynamics results: joint forces/torques')
% end