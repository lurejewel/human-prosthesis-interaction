function modelInfo = forward_simulation(model, modelInfo, state)

frameIndex = 1;

% ---- build simulation cache once (avoid per-frame Java lookups / allocations) ----
nMus = numel(modelInfo.st.muscle.names);
muscleSet = model.getMuscles();

% muscle handles
simCache.allMuscles = cell(1, nMus);
for i = 1 : nMus
    simCache.allMuscles{i} = muscleSet.get(i - 1);
end

% controller + per-muscle Constant objects (reuse each frame)
simCache.brain = org.opensim.modeling.PrescribedController.safeDownCast(model.updControllerSet.get(0));
simCache.controlConsts = cell(1, nMus);
for i = 0 : nMus - 1
    c = org.opensim.modeling.Constant(0);
    simCache.brain.prescribeControlForActuator(i, c);
    simCache.controlConsts{i + 1} = c;
end

% body handles for gait phase detection
simCache.bodyPelvis  = model.getBodySet().get('pelvis');
simCache.bodyCalcnR  = model.getBodySet().get('calcn_r');
simCache.bodyCalcnL  = model.getBodySet().get('calcn_l');
simCache.pelvisCOMLocal = simCache.bodyPelvis.getMassCenter();
simCache.calcnRCOMLocal = simCache.bodyCalcnR.getMassCenter();
simCache.calcnLCOMLocal = simCache.bodyCalcnL.getMassCenter();

% force handles for GRF reading
simCache.frcHeelR = model.getForceSet().get('heelR_ground_contact_force');
simCache.frcHeelL = model.getForceSet().get('heelL_ground_contact_force');
simCache.frcToeR  = model.getForceSet().get('toeR_ground_contact_force');
simCache.frcToeL  = model.getForceSet().get('toeL_ground_contact_force');

% precomputed constants
simCache.totalMass = model.getTotalMass(state);
simCache.gravity   = model.getGravity().get(1);
simCache.stanceTh  = 0.23137978;

% initialize manager
manager = org.opensim.modeling.Manager(model);
manager.initialize(state);

for t = modelInfo.st.simInfo.timeSeries % for every frame

    % calculate ground reaction forces and gait phases
    modelInfo = cal_grf(simCache, modelInfo, state, frameIndex);
    modelInfo = cal_gait_phase(simCache, modelInfo, state, frameIndex);

    % calculate muscle excitations based on muscle states and gait phases
    modelInfo.dy.muscle.exc(:, frameIndex + modelInfo.st.muscle.delay) = cal_muscle_excitation(modelInfo, frameIndex);
    if frameIndex <= modelInfo.st.muscle.delay
        modelInfo.dy.muscle.exc(:, frameIndex) = modelInfo.dy.muscle.exc(:, modelInfo.st.muscle.delay + 1);
    end

    % assign muscle excitations via cached Constant objects (in-place value update)
    for i = 0 : nMus - 1
        simCache.controlConsts{i + 1}.setValue(modelInfo.dy.muscle.exc(i + 1, frameIndex));
    end

    % forward simulation
    state = manager.integrate(t);
    state.setTime(t);

    % update & record
    model.realizeDynamics(state);
    if frameIndex < width(modelInfo.st.simInfo.timeSeries)
        modelInfo = update_modelInfo(modelInfo, state, simCache.allMuscles, frameIndex);
    end

    % fall detection
    if modelInfo.dy.stateHistory(modelInfo.st.model.map('pelvis_ty/value'), frameIndex) < 0.6
        break;
    end
    frameIndex = frameIndex + 1;

end
modelInfo.dy.lastTime = t;

end