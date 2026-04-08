function [excitation, featureMatrix, featureMap, maskMatrix, phaseData] = extractStoMuscleFeatures(stoFilePath, muscleName)
% extractStoMuscleFeatures
% Read a SCONE/OpenSim .sto file with OpenSim API and extract:
% 1) T x 1 excitation of a target muscle
% 2) T x N feature matrix for one body side
% 3) 1 x N containers.Map from feature name to column index
% 4) stance/swing split data based on gait state label
%
% Inputs:
%   stoFilePath : path to .sto file
%   muscleName  : e.g. "soleus_r" or "soleus_l"
%
% Outputs:
%   excitation   : T x 1 vector, from tag "<muscleName>.excitation"
%   featureMatrix: T x N matrix containing
%                  - all "_<side>.fiber_length_norm"
%                  - all "_<side>.mtu_force_norm"
%                  - side sagittal joint angles (hip/knee/ankle)
%                  - side sagittal joint angular velocities with _u suffix (hip/knee/ankle)
%                  - pelvis_tilt and pelvis_tilt_u
%   featureMap   : containers.Map, featureMap(featureName) = column index in featureMatrix
%   phaseData    : struct containing split results
%                  - stance.excitation / stance.featureMatrix / stance.state
%                  - swing.excitation / swing.featureMatrix / swing.state
%                  - stanceIdx / swingIdx / state / stateLabel

if ~(ischar(stoFilePath) || isstring(stoFilePath))
    error('stoFilePath must be a char or string path.');
end
if ~(ischar(muscleName) || isstring(muscleName))
    error('muscleName must be a char or string, e.g. "soleus_r".');
end

stoFilePath = char(stoFilePath);
muscleName = char(muscleName);

if ~isfile(stoFilePath)
    error('STO file does not exist: %s', stoFilePath);
end

table = org.opensim.modeling.TimeSeriesTable(stoFilePath);
numRows = table.getNumRows();

labels = stdVectorStringToCell(table.getColumnLabels());

excitationLabel = [muscleName, '.excitation'];
excitation = getColumnOrError(table, excitationLabel, numRows);

% Infer side from muscle name suffix, e.g. soleus_r -> side = 'r'.
sideToken = regexp(muscleName, '_(r|l)$', 'tokens', 'once');
if isempty(sideToken)
    error('muscleName must end with _r or _l, e.g. soleus_r.');
end
side = sideToken{1};
sideSuffix = ['_', side];

stateLabel = pickStateLabelFromSide(side);
state = getColumnOrError(table, stateLabel, numRows);
if any(state < 0 | state > 4)
    error('Gait state column %s contains values outside [0, 4].', stateLabel);
end

isStance = state >= 0 & state <= 1; % Early Stance (0) & Late Stance (1)
isSwing = state >= 2 & state <= 4; % Liftoff (2), Swing (3), & Landing (4)

fiberLabels = labels(endsWith(labels, [sideSuffix, '.fiber_length_norm']));
forceLabels = labels(endsWith(labels, [sideSuffix, '.mtu_force_norm']));

if isempty(fiberLabels)
    error('No fiber_length_norm features found for side %s.', side);
end
if isempty(forceLabels)
    error('No mtu_force_norm features found for side %s.', side);
end

hipLabel = pickExistingLabel(labels, {['hip_flexion_', side], 'hip_flexion', 'hip_flexion_x'}, 'hip angle');
kneeLabel = pickExistingLabel(labels, {['knee_angle_', side], 'knee_angle', 'knee_angle_x'}, 'knee angle');
ankleLabel = pickExistingLabel(labels, {['ankle_angle_', side], 'ankle_angle', 'ankle_angle_x'}, 'ankle angle');
hipULabel = pickExistingLabel(labels, {['hip_flexion_', side, '_u']}, 'hip angle velocity');
kneeULabel = pickExistingLabel(labels, {['knee_angle_', side, '_u']}, 'knee angle velocity');
ankleULabel = pickExistingLabel(labels, {['ankle_angle_', side, '_u']}, 'ankle angle velocity');
pelvisLabel = pickExistingLabel(labels, {'pelvis_tilt'}, 'pelvis tilt');
pelvisULabel = pickExistingLabel(labels, {'pelvis_tilt_u'}, 'pelvis tilt velocity');

featureLabels = [fiberLabels, forceLabels, {pelvisLabel, pelvisULabel, hipLabel, hipULabel, kneeLabel, kneeULabel, ankleLabel, ankleULabel}];
numFeatures = numel(featureLabels);

featureMatrix = zeros(numRows, numFeatures);
for i = 1:numFeatures
    featureMatrix(:, i) = getColumnOrError(table, featureLabels{i}, numRows);
end

featureMap = containers.Map(featureLabels, num2cell(1:numFeatures));

% Build mask matrix M
maskMatrix = buildMask(numel(forceLabels));

% Split all time-series outputs by gait phase for downstream processing.
phaseData = struct();
phaseData.stateLabel = stateLabel;
phaseData.state = state;
phaseData.stanceIdx = find(isStance);
phaseData.swingIdx = find(isSwing);

phaseData.stance = struct();
phaseData.stance.excitation = excitation(isStance, :);
phaseData.stance.featureMatrix = featureMatrix(isStance, :);
phaseData.stance.state = state(isStance, :);

phaseData.swing = struct();
phaseData.swing.excitation = excitation(isSwing, :);
phaseData.swing.featureMatrix = featureMatrix(isSwing, :);
phaseData.swing.state = state(isSwing, :);

end

%% functions

function stateLabel = pickStateLabelFromSide(side)
if strcmp(side, 'l')
    stateLabel = 'leg0_l.state';
elseif strcmp(side, 'r')
    stateLabel = 'leg1_r.state';
else
    error('muscleName must end with _r or _l to infer gait state label.');
end
end

function label = pickExistingLabel(allLabels, candidates, labelDesc)
label = ''; %#ok<NASGU>
for i = 1:numel(candidates)
    if any(strcmp(allLabels, candidates{i}))
        label = candidates{i};
        return;
    end
end
error('Cannot find %s column. Tried: %s', labelDesc, strjoin(candidates, ', '));
end

function data = getColumnOrError(table, label, numRows)
try
    colVec = table.getDependentColumn(label);
catch
    error('Column not found in .sto file: %s', label);
end

data = simtkVectorToArray(colVec, numRows);
end

function out = simtkVectorToArray(simtkVec, expectedRows)
out = zeros(expectedRows, 1);
for i = 1:expectedRows
    out(i) = simtkVec.get(i - 1);
end
end

function labels = stdVectorStringToCell(stdVec)
labels = cell(1, stdVec.size());
for i = 1:stdVec.size()
    labels{i} = char(stdVec.get(i - 1));
end
end
