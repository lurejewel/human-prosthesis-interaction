function modelInfo = forward_simulation(model, modelInfo, state)

frameIndex = 1;
% initialize muscle objects
allMuscles = cell(1, numel(modelInfo.st.muscle.names));
muscleSet = model.getMuscles();
for i = 1 : numel(allMuscles)
    allMuscles{i} = muscleSet.get(i-1);
end
% initialize brain
brain = org.opensim.modeling.PrescribedController.safeDownCast(model.updControllerSet.get(0));
for muscleIndex = 0 : numel(allMuscles)-1
    brain.prescribeControlForActuator(muscleIndex, org.opensim.modeling.Constant(0));
end
% initialize manager
manager = org.opensim.modeling.Manager(model);
manager.initialize(state);

for t = modelInfo.st.simInfo.timeSeries % for every frame

    % calculate ground reaction forces and gait phases
    modelInfo = cal_grf(model, modelInfo, state, frameIndex);
    modelInfo = cal_gait_phase(model, modelInfo, state, frameIndex);

    % calculate muscle excitations based on muscle states and gait phases
    modelInfo.dy.muscle.exc(:,frameIndex + modelInfo.st.muscle.delay) = cal_muscle_excitation(modelInfo, frameIndex);
    if frameIndex <= modelInfo.st.muscle.delay % the muscle excitations in the first few frames = the frames right after the delay
        modelInfo.dy.muscle.exc(:,frameIndex) = modelInfo.dy.muscle.exc(:,modelInfo.st.muscle.delay+1);
    end

    % assign muscle excitations (control signals) through brain (controller) to muscles (actuators)
    for muscleIndex = 0 : numel(allMuscles)-1
        brain.prescribeControlForActuator(muscleIndex, org.opensim.modeling.Constant(modelInfo.dy.muscle.exc(muscleIndex+1,frameIndex)));
    end

    % forward simulation
    state = manager.integrate(t);
    state.setTime(t);

    % update & record
    model.realizeDynamics(state);
    if frameIndex < width(modelInfo.st.simInfo.timeSeries) % 每一帧都update modelInfo,计算开销太大
        modelInfo = update_modelInfo(modelInfo, state, allMuscles, frameIndex);
    end
 
    if modelInfo.dy.stateHistory(modelInfo.st.model.map('pelvis_ty/value'), frameIndex) < 0.6
        % disp(['Model fell down. Maximum distance reached: ' num2str(model.getCoordinateSet.get('pelvis_tx').getValue(state)) ' m.']);
        break;
    end
    frameIndex = frameIndex+1;

end
modelInfo.dy.lastTime = t;

end