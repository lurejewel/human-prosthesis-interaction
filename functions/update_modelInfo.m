function modelInfo = update_modelInfo(model, modelInfo, state, frameIndex)
% Name: update_MuscleInfo
% Description: update the muscle-related information in the Class modelInfo
%   for the muscle excitation in the next loop. The recorded information
%   includes:
%   - modelInfo.state: state of the model
%   - modelInfo.muscleFiberForceATN: normalized fiber force along tendon
%   - modelInfo.muscleFiberLengthN: normalized fiber length
%   Other interested information:
%   - modelInfo.muscleActivations
% Author(s): Jin Wei, Peking U. wjin24@stu.pku.edu.cn

%% state
% modelInfo.state = state;
modelInfo.stateHistory(:,frameIndex-1) = vec_2_mat(state.getY);

%% muscleFiberForceATN 后续改成数组形式
muscles = model.getMuscles;
optForces = modelInfo.OptimalFiberForces;
modelInfo.muscleFiberForcesATN(1, frameIndex) = muscles.get('hamstrings_r').getFiberForceAlongTendon(state) / optForces.hamstrings;
modelInfo.muscleFiberForcesATN(2, frameIndex) = muscles.get('glut_max_r').getFiberForceAlongTendon(state) / optForces.glut_max;
modelInfo.muscleFiberForcesATN(3, frameIndex) = muscles.get('ilipsoas_r').getFiberForceAlongTendon(state) / optForces.ilipsoas;
modelInfo.muscleFiberForcesATN(4, frameIndex) = muscles.get('vasti_r').getFiberForceAlongTendon(state) / optForces.vasti;
modelInfo.muscleFiberForcesATN(5, frameIndex) = muscles.get('gastroc_r').getFiberForceAlongTendon(state) / optForces.gastrocnemius;
modelInfo.muscleFiberForcesATN(6, frameIndex) = muscles.get('soleus_r').getFiberForceAlongTendon(state) / optForces.soleus;
modelInfo.muscleFiberForcesATN(7, frameIndex) = muscles.get('tibia_r').getFiberForceAlongTendon(state) / optForces.tibia;
modelInfo.muscleFiberForcesATN(8, frameIndex) = muscles.get('hamstrings_l').getFiberForceAlongTendon(state) / optForces.hamstrings;
modelInfo.muscleFiberForcesATN(9, frameIndex) = muscles.get('glut_max_l').getFiberForceAlongTendon(state) / optForces.glut_max;
modelInfo.muscleFiberForcesATN(10, frameIndex) = muscles.get('ilipsoas_l').getFiberForceAlongTendon(state) / optForces.ilipsoas;
modelInfo.muscleFiberForcesATN(11, frameIndex) = muscles.get('vasti_l').getFiberForceAlongTendon(state) / optForces.vasti;
modelInfo.muscleFiberForcesATN(12, frameIndex) = muscles.get('gastroc_l').getFiberForceAlongTendon(state) / optForces.gastrocnemius;
modelInfo.muscleFiberForcesATN(13, frameIndex) = muscles.get('soleus_l').getFiberForceAlongTendon(state) / optForces.soleus;
modelInfo.muscleFiberForcesATN(14, frameIndex) = muscles.get('tibia_l').getFiberForceAlongTendon(state) / optForces.tibia;

%% muscleFiberLengthN
modelInfo.muscleFiberLengthN(1, frameIndex) = muscles.get('hamstrings_r').getNormalizedFiberLength(state);
modelInfo.muscleFiberLengthN(2, frameIndex) = muscles.get('glut_max_r').getNormalizedFiberLength(state);
modelInfo.muscleFiberLengthN(3, frameIndex) = muscles.get('ilipsoas_r').getNormalizedFiberLength(state);
modelInfo.muscleFiberLengthN(4, frameIndex) = muscles.get('vasti_r').getNormalizedFiberLength(state);
modelInfo.muscleFiberLengthN(5, frameIndex) = muscles.get('gastroc_r').getNormalizedFiberLength(state);
modelInfo.muscleFiberLengthN(6, frameIndex) = muscles.get('soleus_r').getNormalizedFiberLength(state);
modelInfo.muscleFiberLengthN(7, frameIndex) = muscles.get('tibia_r').getNormalizedFiberLength(state);
modelInfo.muscleFiberLengthN(8, frameIndex) = muscles.get('hamstrings_l').getNormalizedFiberLength(state);
modelInfo.muscleFiberLengthN(9, frameIndex) = muscles.get('glut_max_l').getNormalizedFiberLength(state);
modelInfo.muscleFiberLengthN(10, frameIndex) = muscles.get('ilipsoas_l').getNormalizedFiberLength(state);
modelInfo.muscleFiberLengthN(11, frameIndex) = muscles.get('vasti_l').getNormalizedFiberLength(state);
modelInfo.muscleFiberLengthN(12, frameIndex) = muscles.get('gastroc_l').getNormalizedFiberLength(state);
modelInfo.muscleFiberLengthN(13, frameIndex) = muscles.get('soleus_l').getNormalizedFiberLength(state);
modelInfo.muscleFiberLengthN(14, frameIndex) = muscles.get('tibia_l').getNormalizedFiberLength(state);

%% muscleActivations
modelInfo.muscleActivations(1, frameIndex) = muscles.get('hamstrings_r').getActivation(state);
modelInfo.muscleActivations(2, frameIndex) = muscles.get('glut_max_r').getActivation(state);
modelInfo.muscleActivations(3, frameIndex) = muscles.get('ilipsoas_r').getActivation(state);
modelInfo.muscleActivations(4, frameIndex) = muscles.get('vasti_r').getActivation(state);
modelInfo.muscleActivations(5, frameIndex) = muscles.get('gastroc_r').getActivation(state);
modelInfo.muscleActivations(6, frameIndex) = muscles.get('soleus_r').getActivation(state);
modelInfo.muscleActivations(7, frameIndex) = muscles.get('tibia_r').getActivation(state);
modelInfo.muscleActivations(8, frameIndex) = muscles.get('hamstrings_l').getActivation(state);
modelInfo.muscleActivations(9, frameIndex) = muscles.get('glut_max_l').getActivation(state);
modelInfo.muscleActivations(10, frameIndex) = muscles.get('ilipsoas_l').getActivation(state);
modelInfo.muscleActivations(11, frameIndex) = muscles.get('vasti_l').getActivation(state);
modelInfo.muscleActivations(12, frameIndex) = muscles.get('gastroc_l').getActivation(state);
modelInfo.muscleActivations(13, frameIndex) = muscles.get('soleus_l').getActivation(state);
modelInfo.muscleActivations(14, frameIndex) = muscles.get('tibia_l').getActivation(state);

%% muscleForces
modelInfo.muscleForces(1, frameIndex) = muscles.get('hamstrings_r').getFiberForce(state);
modelInfo.muscleForces(2, frameIndex) = muscles.get('glut_max_r').getFiberForce(state);
modelInfo.muscleForces(3, frameIndex) = muscles.get('ilipsoas_r').getFiberForce(state);
modelInfo.muscleForces(4, frameIndex) = muscles.get('vasti_r').getFiberForce(state);
modelInfo.muscleForces(5, frameIndex) = muscles.get('gastroc_r').getFiberForce(state);
modelInfo.muscleForces(6, frameIndex) = muscles.get('soleus_r').getFiberForce(state);
modelInfo.muscleForces(7, frameIndex) = muscles.get('tibia_r').getFiberForce(state);
modelInfo.muscleForces(8, frameIndex) = muscles.get('hamstrings_l').getFiberForce(state);
modelInfo.muscleForces(9, frameIndex) = muscles.get('glut_max_l').getFiberForce(state);
modelInfo.muscleForces(10, frameIndex) = muscles.get('ilipsoas_l').getFiberForce(state);
modelInfo.muscleForces(11, frameIndex) = muscles.get('vasti_l').getFiberForce(state);
modelInfo.muscleForces(12, frameIndex) = muscles.get('gastroc_l').getFiberForce(state);
modelInfo.muscleForces(13, frameIndex) = muscles.get('soleus_l').getFiberForce(state);
modelInfo.muscleForces(14, frameIndex) = muscles.get('tibia_l').getFiberForce(state);

end