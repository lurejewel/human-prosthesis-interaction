function modelInfo = update_modelInfo(modelInfo, state, allMuscles, frameIndex, fopt)
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

modelInfo.dy.stateHistory(:,frameIndex) = state.getY.getAsMat(); % vec_2_mat(state.getY);

for i = 1 : numel(fopt)
    muscle = allMuscles{i}; % model.getMuscles().get(i-1);
    modelInfo.dy.muscle.fATN(i,frameIndex+1) = muscle.getFiberForceAlongTendon(state) / fopt(i);
    modelInfo.dy.muscle.lCEN(i,frameIndex+1) = muscle.getNormalizedFiberLength(state);
    modelInfo.dy.muscle.vCE(i,frameIndex+1) = muscle.getFiberVelocity(state);
    modelInfo.dy.muscle.act(i,frameIndex+1) = muscle.getActivation(state);
    modelInfo.dy.muscle.fMTU(i,frameIndex+1) = muscle.getFiberForce(state);
    modelInfo.dy.muscle.fCE(i,frameIndex+1) = muscle.getActiveFiberForce(state);

end

end