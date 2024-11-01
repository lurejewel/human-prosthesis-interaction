addpath('assets\', 'model\');
import org.opensim.modeling.*

%% load model
% model = Model('model\coupled_human-prosthesis_model_scaledFinal.osim');
model = Model('model\NOT YET TO USE\Human-prosthesis body.osim');

%% simulate
model.setUseVisualizer(true);
% initState = read_state_from_file('InitStateGait.xlsx', model); % model.initSystem() included ;)
initState = model.initSystem();
finalState = opensimSimulation.simulate(model, initState, 5);