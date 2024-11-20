function [model, modelInfo] = init_first_frame(model, modelInfo)
% Name: init_first_frame
% Description: initialize model states, including coordinate position/ 
%   velocity, muscle activation, muscle fiber length, muscle force, etc.
%   using simulation data from SCONE (same as the MSD work).
%   This has to be changed in the future ---- since it is ridiculous that a
%   predictive forward simulation framework needs experimental data or
%   simulation data from other simulation platform to initialize the
%   model's state. Also this method is not practicle when more muscles or
%   more complex activities are involved.
%   
%   Several questions to be answered:
%   - how are states initialized in SCONE?
%   - can muscle states be determined roughly assuming v = 0?
%   - can muscle states be determined roughly using ID->SO ?
%   - are there huge effects of small change of initial states on the
%   simulation?
%   - is there any new muscle dynamic model?

%% initial states fetched from SCONE
% order of coordinates in model:
% pelvis_list, pelvis_rotation, pelvis_tilt, pelvis_tx, pelvis_ty,
% pelvis_tz, hip_adduction_r, hip_rotation_r, hip_flexion_r,
% hip_adduction_l, hip_rotation_l, hip_flexion_l, knee_flexion_r,
% knee_flexion_l, ankle_dorsiflexion_r, ankle_dorsiflexion_l

% coordinate position （用数组表示是不是更利于后面赋值？）
init.pose.pelvis_list = 0;
init.pose.pelvis_rotation = 0;
init.pose.pelvis_tilt = -0.1029;
init.pose.pelvis_tx = 0;
init.pose.pelvis_ty = 0.95; % 0.9044;
init.pose.pelvis_tz = 0;
init.pose.hip_adduction_r = 0;
init.pose.hip_rotation_r = 0;
init.pose.hip_flexion_r = -0.2363;
init.pose.hip_adduction_l = 0;
init.pose.hip_rotation_l = 0;
init.pose.hip_flexion_l = 0.4395;
init.pose.knee_flexion_r = -0.0950;
init.pose.knee_flexion_l = -0.1229;
init.pose.ankle_dorsiflexion_r = 0.0551;
init.pose.ankle_dorsiflexion_l = -0.0787;

% coordinate velocity
init.velocity.pelvis_list = 0;
init.velocity.pelvis_rotation = 0;
init.velocity.pelvis_tilt = 0.0724;
init.velocity.pelvis_tx = 1.1476;
init.velocity.pelvis_ty = -0.0185;
init.velocity.pelvis_tz = 0;
init.velocity.hip_adduction_r = 0;
init.velocity.hip_rotation_r = 0;
init.velocity.hip_flexion_r = -0.4484;
init.velocity.hip_adduction_l = 0;
init.velocity.hip_rotation_l = 0;
init.velocity.hip_flexion_l = -2.6960;
init.velocity.knee_flexion_r = -1.5861;
init.velocity.knee_flexion_l = 3.3155;
init.velocity.ankle_dorsiflexion_r = -1.2715;
init.velocity.ankle_dorsiflexion_l = -15.6270;
% gait phase
init.phase.phase_r = 2;
init.phase.phase_l = 0;

% order of muscles in model:
% ham_r, glut_r, ili_r, vas_r, gas_r, sol_r, tib_r, (left in same order)

% normalized fiber force of muscle along tendon
init.normalized_fiber_force_along_tendon.hamstrings_r = 0;
init.normalized_fiber_force_along_tendon.glut_max_r = 0;
init.normalized_fiber_force_along_tendon.ilipsoas_r = 0.4551;
init.normalized_fiber_force_along_tendon.vasti_r = 0.0562;
init.normalized_fiber_force_along_tendon.gastroc_r = 0.4319;
init.normalized_fiber_force_along_tendon.soleus_r = 0.2863;
init.normalized_fiber_force_along_tendon.tibia_r = 0.0061;
init.normalized_fiber_force_along_tendon.hamstrings_l = 0.0321;
init.normalized_fiber_force_along_tendon.glut_max_l = 0.1722;
init.normalized_fiber_force_along_tendon.ilipsoas_l = 0.3356;
init.normalized_fiber_force_along_tendon.vasti_l = 0;
init.normalized_fiber_force_along_tendon.gastroc_l = 0;
init.normalized_fiber_force_along_tendon.soleus_l = 0;
init.normalized_fiber_force_along_tendon.tibia_l = 0.0157;

% muscle activation
init.muscle_activation.hamstrings_r = 0.0105;
init.muscle_activation.glut_max_r = 0.1714;
init.muscle_activation.ilipsoas_r = 0.3752;
init.muscle_activation.vasti_r = 0.2832;
init.muscle_activation.gastroc_r = 0.6444;
init.muscle_activation.soleus_r = 0.4998;
init.muscle_activation.tibia_r = 0.0105;
init.muscle_activation.hamstrings_l = 0.0286;
init.muscle_activation.glut_max_l = 0.3260;
init.muscle_activation.ilipsoas_l = 0.4207;
init.muscle_activation.vasti_l = 0.5832;
init.muscle_activation.gastroc_l = 0.0198;
init.muscle_activation.soleus_l = 0.0101;
init.muscle_activation.tibia_l = 0.1408;

% normalized fiber length of muscle
init.normalized_fiber_length.hamstrings_r = 0.7281;
init.normalized_fiber_length.glut_max_r = 0.4600;
init.normalized_fiber_length.ilipsoas_r = 0.9455;
init.normalized_fiber_length.vasti_r = 0.5359;
init.normalized_fiber_length.gastroc_r = 0.9058;
init.normalized_fiber_length.soleus_r = 0.8332;
init.normalized_fiber_length.tibia_r = 0.7775;
init.normalized_fiber_length.hamstrings_l = 1.0465;
init.normalized_fiber_length.glut_max_l = 1.0853;
init.normalized_fiber_length.ilipsoas_l = 0.7032;
init.normalized_fiber_length.vasti_l = 0.7300;
init.normalized_fiber_length.gastroc_l = 1.0853;
init.normalized_fiber_length.soleus_l = 1.0298;
init.normalized_fiber_length.tibia_l = 0.8371;


%% assign init data to model
initState = model.initSystem();

% pose & velocity
initState.updY.set(0, init.pose.pelvis_list); % /jointset/ground_pelvis/pelvis_list/value	
initState.updY.set(1, init.pose.pelvis_rotation); % /jointset/ground_pelvis/pelvis_rotation/value	
initState.updY.set(2, init.pose.pelvis_tilt); % /jointset/ground_pelvis/pelvis_tilt/value	
initState.updY.set(3, init.pose.pelvis_tx); % /jointset/ground_pelvis/pelvis_tx/value	
initState.updY.set(4, init.pose.pelvis_ty); % /jointset/ground_pelvis/pelvis_ty/value	
initState.updY.set(5, init.pose.pelvis_tz); % /jointset/ground_pelvis/pelvis_tz/value	
initState.updY.set(6, init.pose.hip_adduction_r); % /jointset/hip_r/hip_adduction_r/value	
initState.updY.set(7, init.pose.hip_rotation_r); % /jointset/hip_r/hip_rotation_r/value	
initState.updY.set(8, init.pose.hip_flexion_r); % /jointset/hip_r/hip_flexion_r/value	
initState.updY.set(9, init.pose.hip_adduction_l); % /jointset/hip_l/hip_adduction_l/value	
initState.updY.set(10, init.pose.hip_rotation_l); % /jointset/hip_l/hip_rotation_l/value	
initState.updY.set(11, init.pose.hip_flexion_l); % /jointset/hip_l/hip_flexion_l/value	
initState.updY.set(12, init.pose.knee_flexion_r); % /jointset/knee_r/knee_flexion_r/value	
initState.updY.set(13, init.pose.knee_flexion_l); % /jointset/knee_l/knee_flexion_l/value	
initState.updY.set(14, init.pose.ankle_dorsiflexion_r); % /jointset/ankle_r/ankle_dorsiflexion_r/value	
initState.updY.set(15, init.pose.ankle_dorsiflexion_l); % /jointset/ankle_l/ankle_dorsiflexion_l/value	

initState.updY.set(16, init.velocity.pelvis_list); % /jointset/ground_pelvis/pelvis_list/speed	
initState.updY.set(17, init.velocity.pelvis_rotation); % /jointset/ground_pelvis/pelvis_rotation/speed	
initState.updY.set(18, init.velocity.pelvis_tilt); % /jointset/ground_pelvis/pelvis_tilt/speed	
initState.updY.set(19, init.velocity.pelvis_tx); % /jointset/ground_pelvis/pelvis_tx/speed	
initState.updY.set(20, init.velocity.pelvis_ty); % /jointset/ground_pelvis/pelvis_ty/speed	
initState.updY.set(21, init.velocity.pelvis_tz); % /jointset/ground_pelvis/pelvis_tz/speed	
initState.updY.set(22, init.velocity.hip_adduction_r); % /jointset/hip_r/hip_adduction_r/speed	
initState.updY.set(23, init.velocity.hip_rotation_r); % /jointset/hip_r/hip_rotation_r/speed	
initState.updY.set(24, init.velocity.hip_flexion_r); % /jointset/hip_r/hip_flexion_r/speed	
initState.updY.set(25, init.velocity.hip_adduction_l); % /jointset/hip_l/hip_adduction_l/speed	
initState.updY.set(26, init.velocity.hip_rotation_l); % /jointset/hip_l/hip_rotation_l/speed	
initState.updY.set(27, init.velocity.hip_flexion_l); % /jointset/hip_l/hip_flexion_l/speed	
initState.updY.set(28, init.velocity.knee_flexion_r); % /jointset/knee_r/knee_flexion_r/speed	
initState.updY.set(29, init.velocity.knee_flexion_l); % /jointset/knee_l/knee_flexion_l/speed	
initState.updY.set(30, init.velocity.ankle_dorsiflexion_r); % /jointset/ankle_r/ankle_dorsiflexion_r/speed	
initState.updY.set(31, init.velocity.ankle_dorsiflexion_l); % /jointset/ankle_l/ankle_dorsiflexion_l/speed	

% muscle activation & fiber length
optimalFiberLength.ham = model.getMuscles.get(0).getOptimalFiberLength;
optimalFiberLength.glu = model.getMuscles.get(1).getOptimalFiberLength;
optimalFiberLength.ili = model.getMuscles.get(2).getOptimalFiberLength;
optimalFiberLength.vas = model.getMuscles.get(3).getOptimalFiberLength;
optimalFiberLength.gas = model.getMuscles.get(4).getOptimalFiberLength;
optimalFiberLength.sol = model.getMuscles.get(5).getOptimalFiberLength;
optimalFiberLength.tib = model.getMuscles.get(6).getOptimalFiberLength;

initState.updY.set(32, init.muscle_activation.hamstrings_r); % /forceset/hamstrings_r/activation	
initState.updY.set(33, init.normalized_fiber_length.hamstrings_r * optimalFiberLength.ham); % /forceset/hamstrings_r/fiber_length	
initState.updY.set(34, init.muscle_activation.glut_max_r); % /forceset/glut_max_r/activation	
initState.updY.set(35, init.normalized_fiber_length.glut_max_r * optimalFiberLength.glu); % /forceset/glut_max_r/fiber_length	
initState.updY.set(36, init.muscle_activation.ilipsoas_r); % /forceset/ilipsoas_r/activation
initState.updY.set(37, init.normalized_fiber_length.ilipsoas_r * optimalFiberLength.ili); % /forceset/ilipsoas_r/fiber_length	
initState.updY.set(38, init.muscle_activation.vasti_r); % /forceset/vasti_r/activation	
initState.updY.set(39, init.normalized_fiber_length.vasti_r * optimalFiberLength.vas); % /forceset/vasti_r/fiber_length	
initState.updY.set(40, init.muscle_activation.gastroc_r); % /forceset/gastroc_r/activation	
initState.updY.set(41, init.normalized_fiber_length.gastroc_r * optimalFiberLength.gas); % /forceset/gastroc_r/fiber_length	
initState.updY.set(42, init.muscle_activation.soleus_r); % /forceset/soleus_r/activation	
initState.updY.set(43, init.normalized_fiber_length.soleus_r * optimalFiberLength.sol); % /forceset/soleus_r/fiber_length	
initState.updY.set(44, init.muscle_activation.tibia_r); % /forceset/tibia_r/activation	
initState.updY.set(45, init.normalized_fiber_length.tibia_r * optimalFiberLength.tib); % /forceset/tibia_r/fiber_length	

initState.updY.set(46, init.muscle_activation.hamstrings_l); % /forceset/hamstrings_l/activation	
initState.updY.set(47, init.normalized_fiber_length.hamstrings_l * optimalFiberLength.ham); % /forceset/hamstrings_l/fiber_length	
initState.updY.set(48, init.muscle_activation.glut_max_l); % /forceset/glut_max_l/activation	
initState.updY.set(49, init.normalized_fiber_length.glut_max_l * optimalFiberLength.glu); % /forceset/glut_max_l/fiber_length	
initState.updY.set(50, init.muscle_activation.ilipsoas_l); % /forceset/ilipsoas_l/activation	
initState.updY.set(51, init.normalized_fiber_length.ilipsoas_l * optimalFiberLength.ili); % /forceset/ilipsoas_l/fiber_length	
initState.updY.set(52, init.muscle_activation.vasti_l); % /forceset/vasti_l/activation	
initState.updY.set(53, init.normalized_fiber_length.vasti_l * optimalFiberLength.vas); % /forceset/vasti_l/fiber_length	
initState.updY.set(54, init.muscle_activation.gastroc_l); % /forceset/gastroc_l/activation	
initState.updY.set(55, init.normalized_fiber_length.gastroc_l * optimalFiberLength.gas); % /forceset/gastroc_l/fiber_length	
initState.updY.set(56, init.muscle_activation.soleus_l); % /forceset/soleus_l/activation	
initState.updY.set(57, init.normalized_fiber_length.soleus_l * optimalFiberLength.sol); % /forceset/soleus_l/fiber_length	
initState.updY.set(58, init.muscle_activation.tibia_l); % /forceset/tibia_l/activation	
initState.updY.set(59, init.normalized_fiber_length.tibia_l * optimalFiberLength.tib); % /forceset/tibia_l/fiber_length

%% model final connection
model.realizeDynamics(initState); % call model.realizeDynamics(state) first, then call model.equilibrateMuscles(state), according to gpt.
model.equilibrateMuscles(initState); % after model.equilibrateMuscles(state), the realization stage of the model goes back to Velcotiy
model.realizeDynamics(initState); % call model.realizeDynamics(state) again for Dynamics-stage realization

%% assign init data to modelInfo
modelInfo.state = initState;
% normalized fiber force along tendon
modelInfo.muscleFiberForcesATN(1,1) = init.normalized_fiber_force_along_tendon.hamstrings_r;
modelInfo.muscleFiberForcesATN(2,1) = init.normalized_fiber_force_along_tendon.glut_max_r;
modelInfo.muscleFiberForcesATN(3,1) = init.normalized_fiber_force_along_tendon.ilipsoas_r;
modelInfo.muscleFiberForcesATN(4,1) = init.normalized_fiber_force_along_tendon.vasti_r;
modelInfo.muscleFiberForcesATN(5,1) = init.normalized_fiber_force_along_tendon.gastroc_r;
modelInfo.muscleFiberForcesATN(6,1) = init.normalized_fiber_force_along_tendon.soleus_r;
modelInfo.muscleFiberForcesATN(7,1) = init.normalized_fiber_force_along_tendon.tibia_r;
modelInfo.muscleFiberForcesATN(8,1) = init.normalized_fiber_force_along_tendon.hamstrings_l;
modelInfo.muscleFiberForcesATN(9,1) = init.normalized_fiber_force_along_tendon.glut_max_l;
modelInfo.muscleFiberForcesATN(10,1) = init.normalized_fiber_force_along_tendon.ilipsoas_l;
modelInfo.muscleFiberForcesATN(11,1) = init.normalized_fiber_force_along_tendon.vasti_l;
modelInfo.muscleFiberForcesATN(12,1) = init.normalized_fiber_force_along_tendon.gastroc_l;
modelInfo.muscleFiberForcesATN(13,1) = init.normalized_fiber_force_along_tendon.soleus_l;
modelInfo.muscleFiberForcesATN(14,1) = init.normalized_fiber_force_along_tendon.tibia_l;

% normalized fiber length
modelInfo.muscleFiberLengthN(1,1) = init.normalized_fiber_length.hamstrings_r;
modelInfo.muscleFiberLengthN(2,1) = init.normalized_fiber_length.glut_max_r;
modelInfo.muscleFiberLengthN(3,1) = init.normalized_fiber_length.ilipsoas_r;
modelInfo.muscleFiberLengthN(4,1) = init.normalized_fiber_length.vasti_r;
modelInfo.muscleFiberLengthN(5,1) = init.normalized_fiber_length.gastroc_r;
modelInfo.muscleFiberLengthN(6,1) = init.normalized_fiber_length.soleus_r;
modelInfo.muscleFiberLengthN(7,1) = init.normalized_fiber_length.tibia_r;
modelInfo.muscleFiberLengthN(8,1) = init.normalized_fiber_length.hamstrings_l;
modelInfo.muscleFiberLengthN(9,1) = init.normalized_fiber_length.glut_max_l;
modelInfo.muscleFiberLengthN(10,1) = init.normalized_fiber_length.ilipsoas_l;
modelInfo.muscleFiberLengthN(11,1) = init.normalized_fiber_length.vasti_l;
modelInfo.muscleFiberLengthN(12,1) = init.normalized_fiber_length.gastroc_l;
modelInfo.muscleFiberLengthN(13,1) = init.normalized_fiber_length.soleus_l;
modelInfo.muscleFiberLengthN(14,1) = init.normalized_fiber_length.tibia_l;


end