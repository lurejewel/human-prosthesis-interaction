function modelInfo = forward_simulation(model, modelInfo, state)

import org.opensim.modeling.*

frameIndex = 1;
for t = modelInfo.timeSeries % for every frame

    % calculate ground reaction forces, including normal forces and
    % friction forces, exerted by ground on each foot (calcn_r/l)
    modelInfo = cal_grf(model, modelInfo, state, frameIndex); % [QUESTION: IS THIS UGLY AND INEFFICIENT???->May be integrated in a big Class later.]

    % calculate gait phase based on the state of the model
    modelInfo = cal_gait_phase(model, modelInfo, state, frameIndex);

    % calculate muscle excitations based on muscle states and gait phases
    modelInfo.muscleExcitations(:,frameIndex + modelInfo.muscleReflexDelay) = cal_muscle_excitation(modelInfo, state, frameIndex);
    if frameIndex <= modelInfo.muscleReflexDelay % the muscle excitations in the first few frames = the frames right after the delay (10 frames if delay=10ms)
        modelInfo.muscleExcitations(:,frameIndex) = modelInfo.muscleExcitations(:,modelInfo.muscleReflexDelay+1);
    end

    % assign muscle excitations (control signals) through brain (controller) to muscles (actuators)
    brain = model.updControllerSet.get(0);
    brain = PrescribedController.safeDownCast(brain);
    for muscleIndex = 0 : model.getMuscles.getSize - 1
        brain.prescribeControlForActuator(muscleIndex, Constant(modelInfo.muscleExcitations(muscleIndex+1,frameIndex)));
    end

    % forward simulation 这里的逻辑对吗？需要单独检查一下
    manager = Manager(model);
    manager.initialize(state);
    finalState = manager.integrate(t); % 这里是从哪个初始状态（是否会因为model.initialize而清零），向前推进了多少时间？
    finalState.setTime(t);

    % update & record
    state = finalState;
    model.realizeDynamics(state);
    % modelInfo.stateHistory(:,frameIndex) = vec_2_mat(state.getY);
    frameIndex = frameIndex + 1;
    if frameIndex <= width(modelInfo.timeSeries)
        modelInfo = update_modelInfo(model, modelInfo, state, frameIndex);
    end

    if model.getCoordinateSet.get('pelvis_ty').getValue(state) < 0.6
        % disp(['Model fell down. Maximum distance reached: ' num2str(model.getCoordinateSet.get('pelvis_tx').getValue(state)) ' m.']);
        break;
    end

end
modelInfo.lastTime = t;

end