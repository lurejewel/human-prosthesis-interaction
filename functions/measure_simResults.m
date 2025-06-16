function fit = measure_simResults(t, optConfig, simConfig, modelInfo)
% Name: measure_simResults
% Description: evaluate the simulation results of the model from the
% aspects of walking speed, metabolic expenditure, joint hyperextension,
% etc.
% Author(s): Jin Wei, Peking U. wjin24@stu.pku.edu.cn

if optConfig.hyperPara.stage == 1
    lastFrameIndex = find(all(~isnan(modelInfo.stateHistory)), 1, 'last'); % 为什么要加all？
    fit = -modelInfo.stateHistory(2, lastFrameIndex); % -pelvis_tx in meters 这里为什么不是正的？

elseif optConfig.hyperPara.stage ==2
    error('stage 2 measurement not defined yet.')

else
    error('invalid optimization stage information.')
end

end