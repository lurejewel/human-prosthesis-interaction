function modelInfo = update_MuscleInfo(model, modelInfo, state, frameIndex)
% Name: update_MuscleInfo
% Description: update the muscle-related information in the Class modelInfo
%   for the muscle excitation in the next loop. The recorded information
%   includes:
%   - modelInfo.state: state of the model
%   - modelInfo.muscleFiberForceATN: normalized fiber force along tendon
%   - modelInfo.muscleFiberLengthN: normalized fiber length
%   Other interested information:
%   - modelInfo.muscleActivations
% Author(s): Jin Wei, Peking U.

%% state
modelInfo.state = state;

%% muscleFiberForceATN后续可以改成数组形式
muscles = model.getMuscles;
optForces = modelInfo.OptimalFiberForces;
modelInfo.muscleFiberForcesATN(1, frameIndex) = muscles.get('hamstrings_r').getFiberForce(state) / optForces.hamstrings;
modelInfo.muscleFiberForcesATN(2, frameIndex) = muscles.get('glut_max_r').getFiberForce(state) / optForces.glut_max;
modelInfo.muscleFiberForcesATN(3, frameIndex) = muscles.get('ilipsoas_r').getFiberForce(state) / optForces.ilipsoas;
modelInfo.muscleFiberForcesATN(4, frameIndex) = muscles.get('vasti_r').getFiberForce(state) / optForces.vasti;
modelInfo.muscleFiberForcesATN(5, frameIndex) = muscles.get('gastroc_r').getFiberForce(state) / optForces.gastrocnemius;
modelInfo.muscleFiberForcesATN(6, frameIndex) = muscles.get('soleus_r').getFiberForce(state) / optForces.soleus;
modelInfo.muscleFiberForcesATN(7, frameIndex) = muscles.get('tibia_r').getFiberForce(state) / optForces.tibia;
modelInfo.muscleFiberForcesATN(8, frameIndex) = muscles.get('hamstrings_l').getFiberForce(state) / optForces.hamstrings;
modelInfo.muscleFiberForcesATN(9, frameIndex) = muscles.get('glut_max_l').getFiberForce(state) / optForces.glut_max;
modelInfo.muscleFiberForcesATN(10, frameIndex) = muscles.get('ilipsoas_l').getFiberForce(state) / optForces.ilipsoas;
modelInfo.muscleFiberForcesATN(11, frameIndex) = muscles.get('vasti_l').getFiberForce(state) / optForces.vasti;
modelInfo.muscleFiberForcesATN(12, frameIndex) = muscles.get('gastroc_l').getFiberForce(state) / optForces.gastrocnemius;
modelInfo.muscleFiberForcesATN(13, frameIndex) = muscles.get('soleus_l').getFiberForce(state) / optForces.soleus;
modelInfo.muscleFiberForcesATN(14, frameIndex) = muscles.get('tibia_l').getFiberForce(state) / optForces.tibia;

%% muscleFiberLengthN
optLens = modelInfo.OptimalFiberLengths;
modelInfo.muscleFiberLengthN(1, frameIndex) = muscles.get('hamstrings_r').getFiberLength(state) / optLens.hamstrings;
modelInfo.muscleFiberLengthN(2, frameIndex) = muscles.get('glut_max_r').getFiberLength(state) / optLens.glut_max;
modelInfo.muscleFiberLengthN(3, frameIndex) = muscles.get('ilipsoas_r').getFiberLength(state) / optLens.ilipsoas;
modelInfo.muscleFiberLengthN(4, frameIndex) = muscles.get('vasti_r').getFiberLength(state) / optLens.vasti;
modelInfo.muscleFiberLengthN(5, frameIndex) = muscles.get('gastroc_r').getFiberLength(state) / optLens.gastrocnemius;
modelInfo.muscleFiberLengthN(6, frameIndex) = muscles.get('soleus_r').getFiberLength(state) / optLens.soleus;
modelInfo.muscleFiberLengthN(7, frameIndex) = muscles.get('tibia_r').getFiberLength(state) / optLens.tibia;
modelInfo.muscleFiberLengthN(8, frameIndex) = muscles.get('hamstrings_l').getFiberLength(state) / optLens.hamstrings;
modelInfo.muscleFiberLengthN(9, frameIndex) = muscles.get('glut_max_l').getFiberLength(state) / optLens.glut_max;
modelInfo.muscleFiberLengthN(10, frameIndex) = muscles.get('ilipsoas_l').getFiberLength(state) / optLens.ilipsoas;
modelInfo.muscleFiberLengthN(11, frameIndex) = muscles.get('vasti_l').getFiberLength(state) / optLens.vasti;
modelInfo.muscleFiberLengthN(12, frameIndex) = muscles.get('gastroc_l').getFiberLength(state) / optLens.gastrocnemius;
modelInfo.muscleFiberLengthN(13, frameIndex) = muscles.get('soleus_l').getFiberLength(state) / optLens.soleus;
modelInfo.muscleFiberLengthN(14, frameIndex) = muscles.get('tibia_l').getFiberLength(state) / optLens.tibia;

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

%% muscleForces(along tendon)
modelInfo.muscleForces(1, frameIndex) = muscles.get('hamstrings_r').getFiberForceAlongTendon(state);
modelInfo.muscleForces(2, frameIndex) = muscles.get('glut_max_r').getFiberForceAlongTendon(state);
modelInfo.muscleForces(3, frameIndex) = muscles.get('ilipsoas_r').getFiberForceAlongTendon(state);
modelInfo.muscleForces(4, frameIndex) = muscles.get('vasti_r').getFiberForceAlongTendon(state);
modelInfo.muscleForces(5, frameIndex) = muscles.get('gastroc_r').getFiberForceAlongTendon(state);
modelInfo.muscleForces(6, frameIndex) = muscles.get('soleus_r').getFiberForceAlongTendon(state);
modelInfo.muscleForces(7, frameIndex) = muscles.get('tibia_r').getFiberForceAlongTendon(state);
modelInfo.muscleForces(8, frameIndex) = muscles.get('hamstrings_l').getFiberForceAlongTendon(state);
modelInfo.muscleForces(9, frameIndex) = muscles.get('glut_max_l').getFiberForceAlongTendon(state);
modelInfo.muscleForces(10, frameIndex) = muscles.get('ilipsoas_l').getFiberForceAlongTendon(state);
modelInfo.muscleForces(11, frameIndex) = muscles.get('vasti_l').getFiberForceAlongTendon(state);
modelInfo.muscleForces(12, frameIndex) = muscles.get('gastroc_l').getFiberForceAlongTendon(state);
modelInfo.muscleForces(13, frameIndex) = muscles.get('soleus_l').getFiberForceAlongTendon(state);
modelInfo.muscleForces(14, frameIndex) = muscles.get('tibia_l').getFiberForceAlongTendon(state);

end