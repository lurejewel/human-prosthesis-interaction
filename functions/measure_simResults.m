function fit = measure_simResults(t, optConfig, simConfig, modelInfo)
% Name: measure_simResults
% Description: evaluate the simulation results of the model from the
% aspects of walking speed, metabolic expenditure, joint hyperextension,
% etc.

if optConfig.hyperPara.stage == 1
    lastFrameIndex = find(all(~isnan(modelInfo.stateHistory)), 1, 'last');
    fit = -modelInfo.stateHistory(3, lastFrameIndex); % -pelvis_tx in meters

elseif optConfig.hyperPara.stage ==2
    error('stage 2 measurement not defined yet.')

else
    error('invalid optimization stage information.')
end

end

