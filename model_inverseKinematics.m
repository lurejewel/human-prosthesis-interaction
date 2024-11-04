% -------------------------------------------------------------------------
% Name: model_inverseKinematics.m
% Description: calculate model joint values according to marker data (.trc)
% Info that needs to be specified (in .xml file):
% - inverse kinematics tool:
%   * model file name/path
%   * constraint weight
%   * accuracy
%   * ik task set (in relation to markers)
%   * dynamic marker file name/path(.trc file) & time range
%   * output motion file name/path (.mot)
% Output:
% - xxx.mot: values of DOFs at every frame are recorded
% - xxx_ik_marker_errors.sto: model marker errors from IK
% Note:
% - the angle values in the .mot file is in degrees by default. deg2rad()
%   is needed manually.

% -------------------------------------------------------------------------

addpath('assets\', 'model\')
import org.opensim.modeling.*

%% run Inverse Kinematics Tool and display logger
loggerFile = fopen('opensim.log','rt');
fseek(loggerFile,0,'eof');

model = Model('model\coupled_human-prosthesis_model_scaledFinal.osim');
model.setUseVisualizer(true);
state = model.initSystem();
ikTool = InverseKinematicsTool('ik_setup.xml'); % configure ik tool
ikTool.setModel(model); 
ikTool.run(); % RUN INVERSE KINEMATICS

while ~feof(loggerFile)
    line = fgetl(loggerFile);
    disp(line)
end
fclose(loggerFile);

movefile('subject01_ik_marker_errors.sto','assets\'); % this .sto file is auto generated, which has to be moved manually

%% read .mot file and play it

motData = Storage('assets\subject01_walk1_ik.mot');
rate = motData.getSize / (motData.getLastTime - motData.getFirstTime);
for t = 0 : motData.getSize - 1 % for each frame/state
    stateVec = motData.getStateVector(t).getData;
    for i = 0 : stateVec.size - 1 % for each DOF
        stateVal = stateVec.get(i);
        stateName = model.getStateVariableNames().get(2*i); % 2*i: DOF value; 2*i+1: DOF speed
        if i<3 || i>5 % i=3~5 is translational dof, no need for conversion
            stateVal = deg2rad(stateVal);
        end
        model.setStateVariableValue(state, stateName, stateVal);
    end
    model.updVisualizer().show(state);
    pause(1/rate)
end