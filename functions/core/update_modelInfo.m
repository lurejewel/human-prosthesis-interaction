function modelInfo = update_modelInfo(modelInfo, state, allMuscles, frameIndex)
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

modelInfo.dy.stateHistory(:,frameIndex) = vec_2_mat(state.getY);

fopt = modelInfo.st.muscle.fopt;
for i = 1 : numel(fopt)
    muscle = allMuscles{i}; % model.getMuscles().get(i-1);
    modelInfo.dy.muscle.fATN(i,frameIndex+1) = muscle.getFiberForceAlongTendon(state) / fopt(i);
    modelInfo.dy.muscle.lCEN(i,frameIndex+1) = muscle.getNormalizedFiberLength(state);
    modelInfo.dy.muscle.vCE(i,frameIndex+1) = muscle.getFiberVelocity(state);
    modelInfo.dy.muscle.act(i,frameIndex+1) = muscle.getActivation(state);
    modelInfo.dy.muscle.fMTU(i,frameIndex+1) = muscle.getFiberForce(state);
    modelInfo.dy.muscle.fCE(i,frameIndex+1) = muscle.getActiveFiberForce(state);

end

% modelInfo.muscleFiberForcesATN(1, frameIndex) = muscles.get('hamstrings_r').getFiberForceAlongTendon(state) / optForces.hamstrings;
% modelInfo.muscleFiberForcesATN(2, frameIndex) = muscles.get('glut_max_r').getFiberForceAlongTendon(state) / optForces.glut_max;
% modelInfo.muscleFiberForcesATN(3, frameIndex) = muscles.get('iliopsoas_r').getFiberForceAlongTendon(state) / optForces.iliopsoas;
% modelInfo.muscleFiberForcesATN(4, frameIndex) = muscles.get('vasti_r').getFiberForceAlongTendon(state) / optForces.vasti;
% modelInfo.muscleFiberForcesATN(5, frameIndex) = muscles.get('gastroc_r').getFiberForceAlongTendon(state) / optForces.gastrocnemius;
% modelInfo.muscleFiberForcesATN(6, frameIndex) = muscles.get('soleus_r').getFiberForceAlongTendon(state) / optForces.soleus;
% modelInfo.muscleFiberForcesATN(7, frameIndex) = muscles.get('tibia_r').getFiberForceAlongTendon(state) / optForces.tibia;
% modelInfo.muscleFiberForcesATN(8, frameIndex) = muscles.get('hamstrings_l').getFiberForceAlongTendon(state) / optForces.hamstrings;
% modelInfo.muscleFiberForcesATN(9, frameIndex) = muscles.get('glut_max_l').getFiberForceAlongTendon(state) / optForces.glut_max;
% modelInfo.muscleFiberForcesATN(10, frameIndex) = muscles.get('iliopsoas_l').getFiberForceAlongTendon(state) / optForces.iliopsoas;
% modelInfo.muscleFiberForcesATN(11, frameIndex) = muscles.get('vasti_l').getFiberForceAlongTendon(state) / optForces.vasti;
% modelInfo.muscleFiberForcesATN(12, frameIndex) = muscles.get('gastroc_l').getFiberForceAlongTendon(state) / optForces.gastrocnemius;
% modelInfo.muscleFiberForcesATN(13, frameIndex) = muscles.get('soleus_l').getFiberForceAlongTendon(state) / optForces.soleus;
% modelInfo.muscleFiberForcesATN(14, frameIndex) = muscles.get('tibia_l').getFiberForceAlongTendon(state) / optForces.tibia;
% 
% %% muscleFiberLengthN
% modelInfo.muscleFiberLengthN(1, frameIndex) = muscles.get('hamstrings_r').getNormalizedFiberLength(state);
% modelInfo.muscleFiberLengthN(2, frameIndex) = muscles.get('glut_max_r').getNormalizedFiberLength(state);
% modelInfo.muscleFiberLengthN(3, frameIndex) = muscles.get('iliopsoas_r').getNormalizedFiberLength(state);
% modelInfo.muscleFiberLengthN(4, frameIndex) = muscles.get('vasti_r').getNormalizedFiberLength(state);
% modelInfo.muscleFiberLengthN(5, frameIndex) = muscles.get('gastroc_r').getNormalizedFiberLength(state);
% modelInfo.muscleFiberLengthN(6, frameIndex) = muscles.get('soleus_r').getNormalizedFiberLength(state);
% modelInfo.muscleFiberLengthN(7, frameIndex) = muscles.get('tibia_r').getNormalizedFiberLength(state);
% modelInfo.muscleFiberLengthN(8, frameIndex) = muscles.get('hamstrings_l').getNormalizedFiberLength(state);
% modelInfo.muscleFiberLengthN(9, frameIndex) = muscles.get('glut_max_l').getNormalizedFiberLength(state);
% modelInfo.muscleFiberLengthN(10, frameIndex) = muscles.get('iliopsoas_l').getNormalizedFiberLength(state);
% modelInfo.muscleFiberLengthN(11, frameIndex) = muscles.get('vasti_l').getNormalizedFiberLength(state);
% modelInfo.muscleFiberLengthN(12, frameIndex) = muscles.get('gastroc_l').getNormalizedFiberLength(state);
% modelInfo.muscleFiberLengthN(13, frameIndex) = muscles.get('soleus_l').getNormalizedFiberLength(state);
% modelInfo.muscleFiberLengthN(14, frameIndex) = muscles.get('tibia_l').getNormalizedFiberLength(state);
% 
% %% muscleFiberVelocity
% modelInfo.muscleFiberVelocity(1, frameIndex) = muscles.get('hamstrings_r').getFiberVelocity(state);
% modelInfo.muscleFiberVelocity(2, frameIndex) = muscles.get('glut_max_r').getFiberVelocity(state);
% modelInfo.muscleFiberVelocity(3, frameIndex) = muscles.get('iliopsoas_r').getFiberVelocity(state);
% modelInfo.muscleFiberVelocity(4, frameIndex) = muscles.get('vasti_r').getFiberVelocity(state);
% modelInfo.muscleFiberVelocity(5, frameIndex) = muscles.get('gastroc_r').getFiberVelocity(state);
% modelInfo.muscleFiberVelocity(6, frameIndex) = muscles.get('soleus_r').getFiberVelocity(state);
% modelInfo.muscleFiberVelocity(7, frameIndex) = muscles.get('tibia_r').getFiberVelocity(state);
% modelInfo.muscleFiberVelocity(8, frameIndex) = muscles.get('hamstrings_l').getFiberVelocity(state);
% modelInfo.muscleFiberVelocity(9, frameIndex) = muscles.get('glut_max_l').getFiberVelocity(state);
% modelInfo.muscleFiberVelocity(10, frameIndex) = muscles.get('iliopsoas_l').getFiberVelocity(state);
% modelInfo.muscleFiberVelocity(11, frameIndex) = muscles.get('vasti_l').getFiberVelocity(state);
% modelInfo.muscleFiberVelocity(12, frameIndex) = muscles.get('gastroc_l').getFiberVelocity(state);
% modelInfo.muscleFiberVelocity(13, frameIndex) = muscles.get('soleus_l').getFiberVelocity(state);
% modelInfo.muscleFiberVelocity(14, frameIndex) = muscles.get('tibia_l').getFiberVelocity(state);
% 
% %% muscleActivations
% modelInfo.muscleActivations(1, frameIndex) = muscles.get('hamstrings_r').getActivation(state);
% modelInfo.muscleActivations(2, frameIndex) = muscles.get('glut_max_r').getActivation(state);
% modelInfo.muscleActivations(3, frameIndex) = muscles.get('iliopsoas_r').getActivation(state);
% modelInfo.muscleActivations(4, frameIndex) = muscles.get('vasti_r').getActivation(state);
% modelInfo.muscleActivations(5, frameIndex) = muscles.get('gastroc_r').getActivation(state);
% modelInfo.muscleActivations(6, frameIndex) = muscles.get('soleus_r').getActivation(state);
% modelInfo.muscleActivations(7, frameIndex) = muscles.get('tibia_r').getActivation(state);
% modelInfo.muscleActivations(8, frameIndex) = muscles.get('hamstrings_l').getActivation(state);
% modelInfo.muscleActivations(9, frameIndex) = muscles.get('glut_max_l').getActivation(state);
% modelInfo.muscleActivations(10, frameIndex) = muscles.get('iliopsoas_l').getActivation(state);
% modelInfo.muscleActivations(11, frameIndex) = muscles.get('vasti_l').getActivation(state);
% modelInfo.muscleActivations(12, frameIndex) = muscles.get('gastroc_l').getActivation(state);
% modelInfo.muscleActivations(13, frameIndex) = muscles.get('soleus_l').getActivation(state);
% modelInfo.muscleActivations(14, frameIndex) = muscles.get('tibia_l').getActivation(state);
% 
% %% muscleForces
% modelInfo.muscleForces(1, frameIndex) = muscles.get('hamstrings_r').getFiberForce(state);
% modelInfo.muscleForces(2, frameIndex) = muscles.get('glut_max_r').getFiberForce(state);
% modelInfo.muscleForces(3, frameIndex) = muscles.get('iliopsoas_r').getFiberForce(state);
% modelInfo.muscleForces(4, frameIndex) = muscles.get('vasti_r').getFiberForce(state);
% modelInfo.muscleForces(5, frameIndex) = muscles.get('gastroc_r').getFiberForce(state);
% modelInfo.muscleForces(6, frameIndex) = muscles.get('soleus_r').getFiberForce(state);
% modelInfo.muscleForces(7, frameIndex) = muscles.get('tibia_r').getFiberForce(state);
% modelInfo.muscleForces(8, frameIndex) = muscles.get('hamstrings_l').getFiberForce(state);
% modelInfo.muscleForces(9, frameIndex) = muscles.get('glut_max_l').getFiberForce(state);
% modelInfo.muscleForces(10, frameIndex) = muscles.get('iliopsoas_l').getFiberForce(state);
% modelInfo.muscleForces(11, frameIndex) = muscles.get('vasti_l').getFiberForce(state);
% modelInfo.muscleForces(12, frameIndex) = muscles.get('gastroc_l').getFiberForce(state);
% modelInfo.muscleForces(13, frameIndex) = muscles.get('soleus_l').getFiberForce(state);
% modelInfo.muscleForces(14, frameIndex) = muscles.get('tibia_l').getFiberForce(state);
% 
% %% muscleActiveForces
% modelInfo.muscleActiveForces(1, frameIndex) = muscles.get('hamstrings_r').getActiveFiberForce(state);
% modelInfo.muscleActiveForces(2, frameIndex) = muscles.get('glut_max_r').getActiveFiberForce(state);
% modelInfo.muscleActiveForces(3, frameIndex) = muscles.get('iliopsoas_r').getActiveFiberForce(state);
% modelInfo.muscleActiveForces(4, frameIndex) = muscles.get('vasti_r').getActiveFiberForce(state);
% modelInfo.muscleActiveForces(5, frameIndex) = muscles.get('gastroc_r').getActiveFiberForce(state);
% modelInfo.muscleActiveForces(6, frameIndex) = muscles.get('soleus_r').getActiveFiberForce(state);
% modelInfo.muscleActiveForces(7, frameIndex) = muscles.get('tibia_r').getActiveFiberForce(state);
% modelInfo.muscleActiveForces(8, frameIndex) = muscles.get('hamstrings_l').getActiveFiberForce(state);
% modelInfo.muscleActiveForces(9, frameIndex) = muscles.get('glut_max_l').getActiveFiberForce(state);
% modelInfo.muscleActiveForces(10, frameIndex) = muscles.get('iliopsoas_l').getActiveFiberForce(state);
% modelInfo.muscleActiveForces(11, frameIndex) = muscles.get('vasti_l').getActiveFiberForce(state);
% modelInfo.muscleActiveForces(12, frameIndex) = muscles.get('gastroc_l').getActiveFiberForce(state);
% modelInfo.muscleActiveForces(13, frameIndex) = muscles.get('soleus_l').getActiveFiberForce(state);
% modelInfo.muscleActiveForces(14, frameIndex) = muscles.get('tibia_l').getActiveFiberForce(state);

end