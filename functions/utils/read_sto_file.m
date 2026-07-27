function varargout = read_sto_file(model, stoFilePath, showVideo)
% Name: read_sto_file
% Description: visualize the .sto file (state history of model).
%   Sets state variables by name (via coord.setValue / coord.setSpeedValue /
%   model.setStateVariableValue), NOT by raw Y-index, so the STO column
%   order does not need to match the internal state-vector order.
stoFile = org.opensim.modeling.Storage(stoFilePath, false);
model.setUseVisualizer(true);

labels = stoFile.getColumnLabels;
numLabels = labels.getSize; % including the label 'time'
npts = stoFile.getSize;
time = nan(npts, 1);
data = nan(npts, numLabels-1);
state = model.initSystem();

% ---- Build column handlers: map each STO column to the right setter ----
% colHandler type: 'coordValue' | 'coordSpeed' | 'stateVar' | 'skip'
colHandlers = cell(1, numLabels - 1);
for iLabel = 1 : numLabels - 1
    colName = char(labels.get(iLabel));  % labels are 0‑based; iLabel → column iLabel

    if strcmp(colName, 'totalMetabolic')
        colHandlers{iLabel}.type = 'skip';
    elseif startsWith(colName, '/forceset/')
        % muscle state variable – already a full path
        colHandlers{iLabel}.type = 'stateVar';
        colHandlers{iLabel}.fullPath = colName;
    elseif endsWith(colName, '/value')
        coordName = extractBefore(colName, '/value');
        colHandlers{iLabel}.type = 'coordValue';
        colHandlers{iLabel}.coord = model.getCoordinateSet().get(coordName);
    elseif endsWith(colName, '/speed')
        coordName = extractBefore(colName, '/speed');
        colHandlers{iLabel}.type = 'coordSpeed';
        colHandlers{iLabel}.coord = model.getCoordinateSet().get(coordName);
    else
        warning('read_sto_file: unknown column type "%s" – skipped.', colName);
        colHandlers{iLabel}.type = 'skip';
    end
end

% store data, which can be plot in MATLAB if needed (showVideo = 0/1).
for frameIndex = 1 : npts
    stateVector = stoFile.getStateVector(frameIndex - 1);
    time(frameIndex) = stateVector.getTime;
    for labelIndex = 1 : numLabels - 1
        data(frameIndex, labelIndex) = stateVector.getData.get(labelIndex - 1);
    end
end

if showVideo
    for frameIndex = 1 : npts
        state.setTime(time(frameIndex));

        % ---- set state variables by name (bypass Y-index) ----
        for labelIndex = 1 : numLabels - 1
            val = data(frameIndex, labelIndex);
            handler = colHandlers{labelIndex};

            switch handler.type
                case 'coordValue'
                    handler.coord.setValue(state, val);
                case 'coordSpeed'
                    handler.coord.setSpeedValue(state, val);
                case 'stateVar'
                    model.setStateVariableValue(state, handler.fullPath, val);
                case 'skip'
                    % totalMetabolic – ignore
                otherwise
                    % unknown – ignore
            end
        end

        model.realizeDynamics(state);
        model.updVisualizer.show(state);

        if frameIndex < npts
            pause(time(frameIndex + 1) - time(frameIndex));
        end
    end
end

switch nargout
    case 0

    case 1
        varargout{1} = data;
    case 2
        varargout{1} = time;
        varargout{2} = data;
    case 3
        varargout{1} = time;
        varargout{2} = data;
        labelNames = cell(numLabels, 1);
        for i = 0 : numLabels - 1
            labelNames{i + 1} = char(labels.get(i));
        end
        varargout{3} = labelNames;

end