function modelInfo = cal_grf(simCache, modelInfo, state, frameIndex)

% heelRX = simCache.frcHeelR.getRecordValues(state).get(0);
% heelLX = simCache.frcHeelL.getRecordValues(state).get(0);
% toeRX  = simCache.frcToeR.getRecordValues(state).get(0);
% toeLX  = simCache.frcToeL.getRecordValues(state).get(0);
% heelRY = simCache.frcHeelR.getRecordValues(state).get(1);
% heelLY = simCache.frcHeelL.getRecordValues(state).get(1);
% toeRY  = simCache.frcToeR.getRecordValues(state).get(1);
% toeLY  = simCache.frcToeL.getRecordValues(state).get(1);
% 
% modelInfo.dy.grf.fxr(frameIndex) = heelRX + toeRX;
% modelInfo.dy.grf.fxl(frameIndex) = heelLX + toeLX;
% modelInfo.dy.grf.fyr(frameIndex) = heelRY + toeRY;
% modelInfo.dy.grf.fyl(frameIndex) = heelLY + toeLY;

modelInfo.dy.grf.fxr(frameIndex) = -simCache.grfR.getRecordValues(state).get(0);
modelInfo.dy.grf.fxl(frameIndex) = -simCache.grfL.getRecordValues(state).get(0);
modelInfo.dy.grf.fyr(frameIndex) = -simCache.grfR.getRecordValues(state).get(1);
modelInfo.dy.grf.fyl(frameIndex) = -simCache.grfL.getRecordValues(state).get(1);


end

