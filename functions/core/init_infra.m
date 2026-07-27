function [model, modelInfo] = init_infra(projName, modelStaticProp, act0)
% Name: init_infra
% Description: initialize model states, including coordinate position/ 
%   velocity, muscle activation, muscle fiber length, muscle force, etc.
%   If act0 (14x1) is provided, muscle activations are set to act0
%   before the initial equilibration, yielding a more accurate starting state.

if nargin < 3
    act0 = [];
end

model = org.opensim.modeling.Model(['model/' projName '.osim']);
modelInfo = ModelInfo(modelStaticProp);
model = add_muscle_actuator(model, org.opensim.modeling.PrescribedController()); % add the central controller and muscle actuators

% initial states fetched from a dynamically consistent gait state
% NOTE: initial kinematics are assigned by coordinate NAME, not by state
% vector index, because the index ordering depends on the .osim file
% layout and may differ from the initPose row ordering.

% overwrite initial kinematics (Q and U) — assign by coordinate name
state = model.initSystem();
dofNames = modelInfo.st.model.initPoseDofOrder;
initPose = modelInfo.st.model.initPose;
nDof = numel(dofNames);
coordSet = model.getCoordinateSet();
for j = 1:nDof
    coord = coordSet.get(dofNames{j});
    coord.setValue(state, initPose(2*j - 1));
    coord.setSpeedValue(state, initPose(2*j));
end

% ---- inject activations, equilibrate, and record initial fATN/lCEN ----
modelInfo = equilibrate_and_record_state(model, modelInfo, state, act0);

% capture initial state snapshot for per-particle reset
modelInfo.dy.initStateY = state.getY.getAsMat(); % vec_2_mat(state.getY());

end