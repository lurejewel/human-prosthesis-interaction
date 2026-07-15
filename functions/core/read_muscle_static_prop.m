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
% │   ├── gravity
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
st.muscle.act0Order = names;  % muscle-name order matching act0 (from iterative_static_optimization)
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
gravity = abs(model.getGravity.get(1));
keys = strings(1,state.getNQ +state.getNU+state.getNZ);
for i = 0 : state.getNQ-1 % value & speed of generalized coordinates (interleaved)
    keys(2*i + 1) = string([char(model.getCoordinateSet.get(i)),'/value']);
    keys(2*i + 2) = string([char(model.getCoordinateSet.get(i)),'/speed']);
end
for i = state.getNQ+state.getNU : state.getNQ+state.getNU+state.getNZ-2 % value of variables in forceset
    keys(i+1) = string(model.getStateVariableNames.get(i));
end
keys(state.getNQ+state.getNU+state.getNZ) = "totalMetabolic"; % metabolic expenditure (integrated) of the total body
map = containers.Map(keys, 1:numel(keys));

st.model.totalMass = totalMass;
st.model.gravity = gravity;
st.model.initPose = initPose;
st.model.initPoseDofOrder = dofNames(:);
st.model.map = map;

% ---- build initPoseMap: coordinate name -> initPose index (dofNames order) ----
% Unlike st.model.map (which indexes into the state vector using .osim order),
% this map indexes into initPose using the user-specified dofNames order.
% It is used by extract_lasso_reflex_features to safely read the initial pose.
nDof = numel(dofNames);
initPoseMap = containers.Map('KeyType', 'char', 'ValueType', 'int32');
for j = 1:nDof
    initPoseMap([dofNames{j} '/value']) = 2*j - 1;
    initPoseMap([dofNames{j} '/speed']) = 2*j;
end
st.model.initPoseMap = initPoseMap;

% ---- build permutation between internal state Y order and label order ----
% state.getY() stores variables in SimTK internal order; the map (above)
% uses label order (coordinate-set + getStateVariableNames + totalMetabolic).
% These two orders DIFFER for at least some variables.
% We determine the mapping by perturbing each state variable by name and
% observing which position in state.getY() changes.  No realizeDynamics()
% is called, so exactly 1 element should change per perturbation.
statePerm = model.initSystem();
Y0 = statePerm.getY().getAsMat();
nY = numel(Y0);  % NQ + NU + NZ

permInternalToLabel = zeros(1, nY);  % internalIdx -> labelIdx (1-based)
permLabelToInternal = zeros(1, nY);  % labelIdx   -> internalIdx (1-based)

% Q, U and Z state variables (all except the last entry which is totalMetabolic)
for labelIdx = 1:nY-1
    % full path name in label order (0-based input to getStateVariableNames)
    fullPath = char(model.getStateVariableNames().get(labelIdx - 1));
    testVal  = labelIdx + 0.357;  % unique test value per label

    model.setStateVariableValue(statePerm, fullPath, testVal);
    Y1 = statePerm.getY().getAsMat();

    diffMask = abs(Y1 - Y0) > 1e-12;
    nChanged = nnz(diffMask);
    assert(nChanged == 1, ...
        'Permutation: expected 1 changed entry for "%s", got %d.', fullPath, nChanged);

    internalIdx = find(diffMask, 1);  % 1-based
    permInternalToLabel(internalIdx) = labelIdx;
    permLabelToInternal(labelIdx)    = internalIdx;

    % ---- reset the perturbed variable back to its original value ----
    model.setStateVariableValue(statePerm, fullPath, Y0(internalIdx));
    Y0 = statePerm.getY().getAsMat();  % should be identical to before
end

% totalMetabolic is the last element (identity mapping)
permInternalToLabel(nY) = nY;
permLabelToInternal(nY) = nY;

st.model.permInternalToLabel = permInternalToLabel;
st.model.permLabelToInternal = permLabelToInternal;

%% simulation-level static property

st.simInfo.stepTime = simConfig.stepTime;
st.simInfo.endTime = simConfig.endTime;
st.simInfo.speed = simConfig.speed;
st.simInfo.timeSeries = 0 : simConfig.stepTime : simConfig.endTime;

end