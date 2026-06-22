function modelInfo = update_modelInfo(modelInfo, state, allMuscles, idx, fopt)
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

map = modelInfo.st.model.map;
Y = state.getY.getAsMat();

for i = 1 : numel(fopt)
    muscle = allMuscles{i}; % model.getMuscles().get(i-1);
    modelInfo.dy.muscle.fATN(i,idx+1) = muscle.getFiberForceAlongTendon(state) / fopt(i);
    modelInfo.dy.muscle.lCEN(i,idx+1) = muscle.getNormalizedFiberLength(state);
    % modelInfo.dy.muscle.vCE(i,idx+1) = muscle.getFiberVelocity(state);
    % modelInfo.dy.muscle.act(i,idx+1) = muscle.getActivation(state);
    % modelInfo.dy.muscle.fMTU(i,idx+1) = muscle.getFiberForce(state);
    % modelInfo.dy.muscle.fCE(i,idx+1) = muscle.getActiveFiberForce(state);

    key = ['/forceset/' char(allMuscles{i}.getName) '/activation'];
    Y(map(key)) = muscle.getActivation(state);

end

modelInfo.dy.stateHistory(:,idx) = Y;

end