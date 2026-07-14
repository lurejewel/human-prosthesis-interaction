function save_state_as_sto(model, modelInfo, showVideoFlag)

    state = model.initSystem();
    labels = org.opensim.modeling.ArrayStr();
    labels.append('time');
    for stateIndex = 0 : state.getNQ-1 % value & speed of jointset (interleaved)
        labels.append([char(model.getCoordinateSet.get(stateIndex)) '/value']);
        labels.append([char(model.getCoordinateSet.get(stateIndex)) '/speed']);
    end
    for stateIndex = state.getNQ + state.getNU : state.getNQ + state.getNU + state.getNZ-1 - 1 % value of forceset
        labels.append([char(model.getStateVariableNames.get(stateIndex))]);
    end
    labels.append('totalMetabolic');  % integrated metabolic probe output

    stoFile = org.opensim.modeling.Storage();
    filename = ['results\sim_result_' char(datetime("now","Format","yyyy-MM-dd_HH-mm-ss")) '.sto'];
    stoFile.setName(filename);
    stoFile.setColumnLabels(labels);

    for frameIndex = 1 : find(all(~isnan(modelInfo.dy.labelHistory)), 1, 'last')
        % labelHistory rows are in label order, matching Storage column labels.
        nDataCols = labels.getSize() - 1;   % columns excluding 'time'
        Y = mat_2_vec(modelInfo.dy.labelHistory(1:nDataCols, frameIndex));
        stateVector = org.opensim.modeling.StateVector(modelInfo.st.simInfo.timeSeries(frameIndex), Y);
        stoFile.append(stateVector);
    end

    stoFile.print(filename);
    disp(['Successfully save the data to: ' filename]);
    
    if showVideoFlag
        % play the video
        read_sto_file(model, filename, 1);
    end
    
end