function st = read_muscle_static_prop(projName, simConfig, initPose, dofNames)
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
% │   ├── initPoseDofOrder          cell array, DOF names in initPose row order
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
% Build muscle name -> index map directly from model order (1-based)
map = containers.Map('KeyType', 'char', 'ValueType', 'int32');
for i = 1:nMus
    map(names{i}) = i;
end
delay = round(10/1000/simConfig.stepTime); % 10 ms delay for muscle reflex mechanism to take effect

st.muscle.names = names;
% st.muscle.lambda = lambda;
% st.muscle.mass = mass;
st.muscle.lopt = lopt;
st.muscle.fopt = fopt;
st.muscle.map = map;
st.muscle.delay = delay;

% Build muscle-name → within-leg position map (1-based, used by cal_muscle_excitation)
nMusPerLeg = nMus / 2;
keys = map.keys;
vals = zeros(1, numel(keys));
for j = 1:numel(keys)
    vals(j) = map(keys{j});
end
[~, sortIdx] = sort(vals);
rightKeys = keys(sortIdx(1:nMusPerLeg));
legIdx = containers.Map('KeyType', 'char', 'ValueType', 'int32');
for j = 1:nMusPerLeg
    baseName = regexprep(rightKeys{j}, '_[rl]$', '');
    legIdx(baseName) = j;
end
st.muscle.legIdx = legIdx;

%% model-level static property
state = model.initSystem(); % this is a must before calling model.getTotalMass() and model.getNumStateVariables()
totalMass = model.getTotalMass(state);
keys = strings(1,state.getNQ +state.getNU+state.getNZ);
for i = 0 : state.getNQ-1 % value & speed of generalized coordinates
    keys(i+1) = string([char(model.getCoordinateSet.get(i)),'/value']);
    keys(i+state.getNQ+1) = string([char(model.getCoordinateSet.get(i)),'/speed']);
end
for i = state.getNQ+state.getNU : state.getNQ+state.getNU+state.getNZ-2 % value of variables in forceset
    keys(i+1) = string(model.getStateVariableNames.get(i));
end
keys(state.getNQ+state.getNU+state.getNZ) = "totalMetabolic"; % metabolic expenditure (integrated) of the total body
map = containers.Map(keys, 1:numel(keys));

st.model.totalMass = totalMass;
st.model.initPose = initPose;
st.model.initPoseDofOrder = dofNames(:);
st.model.map = map;

%% simulation-level static property

st.simInfo.stepTime = simConfig.stepTime;
st.simInfo.endTime = simConfig.endTime;
st.simInfo.speed = simConfig.speed;
st.simInfo.timeSeries = 0 : simConfig.stepTime : simConfig.endTime;

end