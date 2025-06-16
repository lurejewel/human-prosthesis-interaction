addpath('assets\', 'model\', 'functions\');
import org.opensim.modeling.*

%% load model
model = Model('model\coupled_human-prosthesis_model_scaledFinal.osim');
% model = Model('model\NOT YET TO USE\Human-prosthesis body.osim');
% model.setUseVisualizer(true);


%% simulate preparation
% simulation configuration
simConfig.startTime = 0;
simConfig.endTime = 5;
simConfig.stepTime = 0.001;
timeSeries = simConfig.startTime+simConfig.stepTime : simConfig.stepTime : simConfig.endTime;

% data to be recorded
heelRFY = nan(1, ceil((simConfig.endTime-simConfig.startTime)/simConfig.stepTime));
groundFY = nan(1, length(heelRFY));
% pelvis_ty = nan(1, length(heelRFY));

% model configuration
model.setUseVisualizer(true);
% initState = read_state_from_file('InitStateGait.xlsx', model); % model.initSystem() included ;)
initState = model.initSystem();
model.updCoordinateSet().get('pelvis_ty').setValue(initState, 1.2); % a little bit higher
% visualizer = model.updVisualizer;
% model.updVisualizer.updSimbodyVisualizer.setShutdownWhenDestructed(true);

%% simulation: step by step
i = 1;
for time = timeSeries
    % finalState = opensimSimulation.simulate(model, initState, time);
    manager = Manager(model); % the Manager Class needs to be called repeatedly in the loop, according to the OpenSim API Guide
    manager.initialize(initState); % this is necessary before simulation (manager.integrate)
    finalState = manager.integrate(time); % note that the para "time" is not duration time but final time

    finalState.setTime(time);
    model.realizeDynamics(finalState);
    
    heelRFY(i) = model.getForceSet.get('heelR_ground_contact_force').getRecordValues(finalState).get(1);
    % groundFY(i) = model.getForceSet.get('heelR_ground_contact_force').getRecordValues(finalState).get(7);
    % pelvis_ty(i) = initState.getQ.get(4);

    initState = finalState;
    i = i + 1;
end
figure, plot(timeSeries, heelRFY), title('GRF of right heel -- normal force')
% figure, plot(timeSeries, pelvis_ty), title('vertical altitude of COM of pelvis')