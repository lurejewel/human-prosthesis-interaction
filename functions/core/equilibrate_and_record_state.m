function modelInfo = equilibrate_and_record_state(model, modelInfo, state, act0)
% equilibrate_and_record_state
% Shared helper for init_infra and reset_particle_state.
% (1) Optionally inject initial muscle activations by muscle name,
% (2) equilibrate muscles and realise dynamics,
% (3) record first-frame fATN and lCEN into modelInfo.dy.muscle.
%
% Input:
%   model     - org.opensim.modeling.Model
%   modelInfo - ModelInfo object (updated in-place)
%   state     - org.opensim.modeling.State (kinematics already set;
%               modified in-place by equilibrate/realise)
%   act0      - (optional) nMus×1 double, initial muscle activations
%               in the order given by modelInfo.st.muscle.names.
% Output:
%   modelInfo - same object, with dy.muscle.fATN(:,1), dy.muscle.lCEN(:,1)

nMus = numel(modelInfo.st.muscle.names);

% ---- optionally inject static-optimisation activations ----
% Assigned by muscle NAME (not ForceSet index) because the iteration
% order may differ from the act0 layout.
if nargin >= 4 && ~isempty(act0)
    for i = 1:nMus
        musc = org.opensim.modeling.Muscle.safeDownCast( ...
            model.getMuscles().get(modelInfo.st.muscle.names{i}));
        musc.setActivation(state, act0(i));
    end
end

% ---- execute dynamics ----
model.equilibrateMuscles(state);
model.realizeDynamics(state);

% ---- record first-frame fATN and lCEN (needed for reflex control) ----
for i = 1:nMus
    modelInfo.dy.muscle.fATN(i, 1) = ...
        model.getMuscles().get(i - 1).getActiveFiberForce(state) / modelInfo.st.muscle.fopt(i);
    modelInfo.dy.muscle.lCEN(i, 1) = ...
        model.getMuscles().get(i - 1).getNormalizedFiberLength(state);
end
end
