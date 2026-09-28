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

% ---- handles for GRF / CoP reading ----
% human0918.osim splits each foot's ground contact into three independent
% HuntCrossleyForce components (heel / toe / tip); foot_contact_defs() is the
% single place where those names live.  For every contact point we cache the
% force handle plus the parent PhysicalFrame and the local offset of its
% ContactSphere, so that cal_grf can sum the forces and weight their
% application points into a CoP without any further Java lookups per frame.
contactDef = foot_contact_defs();
simCache.contactDef  = contactDef;
simCache.footContact = cell(1, numel(contactDef.sides));   % {side}(point)
for i = 1 : numel(contactDef.sides)
    pts = struct('force', cell(1, contactDef.nPoints), ...
                 'frame', cell(1, contactDef.nPoints), ...
                 'loc',   cell(1, contactDef.nPoints));
    for j = 1 : contactDef.nPoints
        pts(j).force = model.getForceSet().get(contactDef.forceNames{i, j});
        cg           = model.getContactGeometrySet().get(contactDef.geomNames{i, j});
        pts(j).frame = cg.getFrame();
        pts(j).loc   = cg.get_location();
    end
    simCache.footContact{i} = pts;
end

% force handles for knee limit force reading
simCache.frcKneeLimitR = model.getForceSet().get('knee_lim_r');
simCache.frcKneeLimitL = model.getForceSet().get('knee_lim_l');

% precomputed constants
simCache.bw = model.getTotalMass(state) * abs(model.getGravity().get(1)); % body weight (in newtons) = total mass (kg) * g (m/s^2)
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
        modelInfo.update(state, simCache.allMuscles, ...
            simCache.frcKneeLimitR, simCache.frcKneeLimitL, frameIndex);
    end

    % fall detection
    if modelInfo.dy.labelHistory(modelInfo.st.model.map('pelvis_ty/value'), frameIndex) < 0.8
        break;
    end
    frameIndex = frameIndex + 1;

end

modelInfo.dy.lastTime = t;
modelInfo.dy.hasRun  = true;
modelInfo.dy.totalMetabolic = model.getProbeSet.get(0).getProbeOutputs(state).get(0);
end