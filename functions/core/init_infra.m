function [model, modelInfo] = init_infra(projName, modelStaticProp)
% Name: init_infra
% Description: initialize model states, including coordinate position/ 
%   velocity, muscle activation, muscle fiber length, muscle force, etc.

model = org.opensim.modeling.Model(['model/' projName '.osim']);
modelInfo = ModelInfo(modelStaticProp);
model = add_muscle_actuator(model, org.opensim.modeling.PrescribedController()); % add the central controller and muscle actuators

% initial states fetched from a dynamically consistent gait state
% order of coordinates in model:
% pelvis_list, pelvis_rotation, pelvis_tilt, pelvis_tx, pelvis_ty,
% pelvis_tz, hip_adduction_r, hip_rotation_r, hip_flexion_r,
% hip_adduction_l, hip_rotation_l, hip_flexion_l, knee_flexion_r,
% knee_flexion_l, ankle_dorsiflexion_r, ankle_dorsiflexion_l

% overwrite initial kinematics (Q and U)
state = model.initSystem();
for i = 0 : state.getNQ+state.getNU-1
    state.updY.set(i, modelInfo.st.model.initPose(i+1));
end

% execute dynamics
model.equilibrateMuscles(state);
model.realizeDynamics(state);

% store fATN and lCEN [needed for the upcoming muscle-reflex control]
 for i = 1 : numel(modelInfo.st.muscle.names)
     modelInfo.dy.muscle.fATN(i,1) = model.getMuscles.get(i-1).getActiveFiberForce(state) / modelInfo.st.muscle.fopt(i);
     modelInfo.dy.muscle.lCEN(i,1) = model.getMuscles.get(i-1).getNormalizedFiberLength(state);
 end

% capture initial state snapshot for per-particle reset
modelInfo.dy.initStateY = state.getY.getAsMat(); % vec_2_mat(state.getY());

end