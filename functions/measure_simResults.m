function fit = measure_simResults(optConfig, simConfig, modelInfo)
% Name: measure_simResults
% Description: evaluate the simulation results of the model from the
% aspects of walking speed, metabolic expenditure, joint hyperextension,
% etc.
% Author(s): Jin Wei, Peking U. wjin24@stu.pku.edu.cn
if modelInfo.lastTime == -1
    error('[CUSTOMIZED ERROR] simulation did not operate normally.')
end

if modelInfo.stage == 1
    lastFrameIndex = find(all(~isnan(modelInfo.stateHistory)), 1, 'last');
    fit = simConfig.endTime*simConfig.speed - modelInfo.stateHistory(2, lastFrameIndex); % desired distance - pelvis_tx in meters

elseif modelInfo.stage ==2
    error('[CUSTOMIZED ERROR] stage 2 measurement not defined yet.')

else
    error('[CUSTOMIZED ERROR] invalid optimization stage information.')
end

end