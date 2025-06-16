function varargout = read_sto_file(modelFilePath, stoFilePath, showVideo)
% Name: read_sto_file
% Description: visualize the .sto file (state history of model).
import org.opensim.modeling.*
stoFile = Storage(stoFilePath, false);
model = Model(modelFilePath);
model.setUseVisualizer(true);

labels = stoFile.getColumnLabels;
numLabels = labels.getSize; % including the label 'time'
npts = stoFile.getSize;
time = nan(npts, 1);
data = nan(npts, numLabels-1);
state = model.initSystem();

% store data, which can be plot in MATLAB if needed (showVideo = 0/1).
for frameIndex = 1 : npts
    stateVector = stoFile.getStateVector(frameIndex-1);
    time(frameIndex) = stateVector.getTime;
    for labelIndex = 1 : numLabels-1
        data(frameIndex, labelIndex) = stateVector.getData.get(labelIndex-1);
    end

end

if showVideo
    for frameIndex = 1 : npts
        state.setTime(time(frameIndex));
        Y = state.updY;
        for labelIndex = 1 : numLabels-1
            Y.set(labelIndex-1, data(frameIndex, labelIndex));
        end
        model.realizeDynamics(state);
        model.updVisualizer.show(state);

        if frameIndex < npts
            pause(time(frameIndex+1) - time(frameIndex));
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
        labelNames = nan(numLabels, 1);
        for i = 0 : numLabels
            labelNames(i+1) = labels.get(i);
        end
        varargout{3} = labelNames;

end