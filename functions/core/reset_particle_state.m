function [state, modelInfo] = reset_particle_state(model, modelInfo)
% Name: reset_particle_state
% Description: Reset the OpenSim state and ModelInfo dynamic fields to
%   their initial values before running a new particle's forward simulation.
%   This ensures each particle starts from an identical, dynamically
%   consistent initial state while reusing the same model instance.
%
% Input:
%   model     - org.opensim.modeling.Model (already initialized)
%   modelInfo - ModelInfo object (st holds static props, dy holds the
%               initStateY snapshot and will have its record partially reset)
% Output:
%   state     - org.opensim.modeling.State, reset to t=0 with initPose

% fresh state and overwrite with initial kinematics
state = model.initSystem();
nQNU = state.getNQ + state.getNU;
initY = modelInfo.dy.initStateY;
for i = 0 : nQNU - 1
    state.updY.set(i, initY(i + 1));
end
state.setTime(0);

% re-establish dynamic consistency
model.equilibrateMuscles(state);
model.realizeDynamics(state);

% reset initial muscle fATN and lCEN (needed for first-frame reflex control)
nMus = numel(modelInfo.st.muscle.names);
for i = 1 : nMus
    modelInfo.dy.muscle.fATN(i, 1) = model.getMuscles.get(i - 1).getActiveFiberForce(state) / modelInfo.st.muscle.fopt(i);
    modelInfo.dy.muscle.lCEN(i, 1) = model.getMuscles.get(i - 1).getNormalizedFiberLength(state);
end
end