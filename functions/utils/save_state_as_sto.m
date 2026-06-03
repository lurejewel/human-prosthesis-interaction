function save_state_as_sto(flag, model, modelInfo)

if flag

    state = model.initSystem();
    labels = org.opensim.modeling.ArrayStr();
    labels.append('time');
    for stateIndex = 0 : state.getNQ-1 % value of jointset
        labels.append([char(model.getCoordinateSet.get(stateIndex)) '/value']);
    end
    for stateIndex = 0 : state.getNU-1 % speed of jointset
        labels.append([char(model.getCoordinateSet.get(stateIndex)) '/speed']);
    end
    for stateIndex = state.getNQ + state.getNU : state.getNQ + state.getNU + state.getNZ-1 % value of forceset
        labels.append([char(model.getStateVariableNames.get(stateIndex))]); % NOTE that the order of state in Model class (model.getState...) is not the same as that in State class (state from model.initSystem())
    end

    stoFile = org.opensim.modeling.Storage();
    filename = ['results\sim_result_' char(datetime("now","Format","yyyy-MM-dd_HH-mm-ss")) '.sto'];
    stoFile.setName(filename);
    stoFile.setColumnLabels(labels);

    for frameIndex = 1 : find(all(~isnan(modelInfo.dy.stateHistory)), 1, 'last')
        Y = mat_2_vec(modelInfo.dy.stateHistory(:,frameIndex));
        stateVector = org.opensim.modeling.StateVector(modelInfo.st.simInfo.timeSeries(frameIndex), Y);
        stoFile.append(stateVector);
    end

    stoFile.print(filename);
    disp(['Successfully save the data to: ' filename]);
    
    % play the video
    % read_sto_file('model/coupled_human-prosthesis_model.osim', filename, 1);

end

end