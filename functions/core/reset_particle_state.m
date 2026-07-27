function [state, modelInfo] = reset_particle_state(model, modelInfo, act0)
% Name: reset_particle_state
% Description: Reset the OpenSim state and ModelInfo dynamic fields to
%   their initial values before running a new particle's forward simulation.
%   This ensures each particle starts from an identical, dynamically
%   consistent initial state while reusing the same model instance.
%   If act0 (14x1) is provided, muscle activations are set before equilibration.
%
% Input:
%   model     - org.opensim.modeling.Model (already initialized)
%   modelInfo - ModelInfo object (st holds static props, dy holds the
%               initStateY snapshot and will have its record partially reset)
%   act0      - (optional) 14x1 double, initial muscle activations
% Output:
%   state     - org.opensim.modeling.State, reset to t=0 with initPose

if nargin < 3
    act0 = [];
end

% fresh state and overwrite with initial kinematics
state = model.initSystem();
nQNU = state.getNQ + state.getNU;
initY = modelInfo.dy.initStateY;
for i = 0 : nQNU - 1
    state.updY.set(i, initY(i + 1));
end
state.setTime(0);

% ---- inject activations, equilibrate, and record initial fATN/lCEN ----
modelInfo = equilibrate_and_record_state(model, modelInfo, state, act0);
end