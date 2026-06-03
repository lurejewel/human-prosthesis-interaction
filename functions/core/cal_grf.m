function modelInfo = cal_grf(simCache, modelInfo, state, frameIndex)

heelRX = simCache.frcHeelR.getRecordValues(state).get(0);
heelLX = simCache.frcHeelL.getRecordValues(state).get(0);
toeRX  = simCache.frcToeR.getRecordValues(state).get(0);
toeLX  = simCache.frcToeL.getRecordValues(state).get(0);
heelRY = simCache.frcHeelR.getRecordValues(state).get(1);
heelLY = simCache.frcHeelL.getRecordValues(state).get(1);
toeRY  = simCache.frcToeR.getRecordValues(state).get(1);
toeLY  = simCache.frcToeL.getRecordValues(state).get(1);

modelInfo.dy.grf.fxr(frameIndex) = heelRX + toeRX;
modelInfo.dy.grf.fxl(frameIndex) = heelLX + toeLX;
modelInfo.dy.grf.fyr(frameIndex) = heelRY + toeRY;
modelInfo.dy.grf.fyl(frameIndex) = heelLY + toeLY;
% modelInfo.grf.frictionR(frameIndex) = heelRX + toeRX;
% modelInfo.grf.frictionL(frameIndex) = heelLX + toeLX;
% modelInfo.grf.normalR(frameIndex) = heelRY + toeRY;
% modelInfo.grf.normalL(frameIndex) = heelLY + toeLY;

end

