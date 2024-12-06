function modelInfo = cal_grf(model, modelInfo, state, frameIndex)

heelRX = model.getForceSet.get('heelR_ground_contact_force').getRecordValues(state).get(0);
heelLX = model.getForceSet.get('heelL_ground_contact_force').getRecordValues(state).get(0);
toeRX  = model.getForceSet.get('toeR_ground_contact_force').getRecordValues(state).get(0);
toeLX  = model.getForceSet.get('toeL_ground_contact_force').getRecordValues(state).get(0);
heelRY = model.getForceSet.get('heelR_ground_contact_force').getRecordValues(state).get(1);
heelLY = model.getForceSet.get('heelL_ground_contact_force').getRecordValues(state).get(1);
toeRY  = model.getForceSet.get('toeR_ground_contact_force').getRecordValues(state).get(1);
toeLY  = model.getForceSet.get('toeL_ground_contact_force').getRecordValues(state).get(1);

modelInfo.grf.frictionR(frameIndex) = heelRX + toeRX;
modelInfo.grf.frictionL(frameIndex) = heelLX + toeLX;
modelInfo.grf.normalR(frameIndex) = heelRY + toeRY;
modelInfo.grf.normalL(frameIndex) = heelLY + toeLY;

end

