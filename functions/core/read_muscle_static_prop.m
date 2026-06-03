function st = read_muscle_static_prop(projName, simConfig, initPose)
% st: static properties of the model, corresponding to obj.st in the
% ModelInfo class.
% st
% ├── muscle
% │   ├── names
% │   ├── lambda                    muscle type I fiber composition
% │   ├── mass
% │   ├── lopt                      optimal fiber length
% │   ├── fopt                      optimal fiber force
% │   └── map                       abbr -> index 目前只支持lambda取值，其他数组需要改写
% │   └── delay                     time delay of electrical signals from the central nervous system to the muscle fibers
% ├── model
% │   ├── totalMass
% │   ├── initPose
% │   └── map
% └── simInfo
%     ├── stepTime
%     ├── endTime
%     ├── speed
%     └── timeSeries


%% muscular-level static property
model = org.opensim.modeling.Model(['model/' projName '.osim']);
nMus = model.getMuscles.getSize;
names = cell(nMus,1);
lopt = nan(nMus,1);
fopt = nan(nMus,1);
for i = 0 : nMus-1
    names{i+1} = char(model.getMuscles.get(i).getName);
    lopt(i+1) = model.getMuscles.get(names{i+1}).getOptimalFiberLength;
    fopt(i+1) = model.getMuscles.get(names{i+1}).getMaxIsometricForce;
end
[lambda, mass, map] = read_muscle_lambda_and_mass(['assets/' projName '_muscleProp.csv']);
delay = round(10/1000/simConfig.stepTime); % 10 ms delay for muscle reflex mechanism to take effect

st.muscle.names = names;
st.muscle.lambda = lambda;
st.muscle.mass = mass;
st.muscle.lopt = lopt;
st.muscle.fopt = fopt;
st.muscle.map = map;
st.muscle.delay = delay;

%% model-level static property
state = model.initSystem(); % this is a must before calling model.getTotalMass() and model.getNumStateVariables()
totalMass = model.getTotalMass(state);
keys = strings(1,state.getNQ +state.getNU+state.getNZ);
for i = 0 : state.getNQ-1 % value & speed of generalized coordinates
    keys(i+1) = string([char(model.getCoordinateSet.get(i)),'/value']);
    keys(i+state.getNQ+1) = string([char(model.getCoordinateSet.get(i)),'/speed']);
end
for i = state.getNQ+state.getNU : state.getNQ+state.getNU+state.getNZ-1 % value of variables in forceset
    keys(i+1) = string(model.getStateVariableNames.get(i));
end
map = containers.Map(keys, 1:numel(keys));

st.model.totalMass = totalMass;
st.model.initPose = initPose;
st.model.map = map;

%% simulation-level static property

st.simInfo.stepTime = simConfig.stepTime;
st.simInfo.endTime = simConfig.endTime;
st.simInfo.speed = simConfig.speed;
st.simInfo.timeSeries = 0 : simConfig.stepTime : simConfig.endTime;

end