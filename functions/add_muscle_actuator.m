function model = add_muscle_actuator(model, brain)
% Name: add_muscle_actuator.m
% Description: assign a *PrescribedController* to the model, and an 
%   *Acutator* for each muscle. The controller can generate control signals
%   (i.e., muscle exctations) in the form of *StepFunction*) responded by 
%   muscles. 
% Reference Code: <installation path>\OpenSim\4.5\Code\Matlab\build_and_
%   simulate_simple_arm.m
% Author(s): Jin Wei, Peking U.

muscleNum = model.getMuscles.getSize;

% brain = PrescribedController(); % central neural system
brain.setName('central_neural_system');
for muscleIndex = 0 : muscleNum-1
    brain.addActuator(model.getMuscles.get(muscleIndex));
end
model.addController(brain);

end