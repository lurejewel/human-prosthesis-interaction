function initState = read_state_from_file(stateFileName, model)

initState = model.initSystem();
stateCell = readcell(stateFileName);

for i = 0 : model.getCoordinateSet().getSize()-1
    [row, col] = find(strcmp(stateCell, model.getCoordinateSet().get(i).getName));
    model.updCoordinateSet().get(i).setValue(initState, stateCell{row, col+1});
    model.updCoordinateSet().get(i).setSpeedValue(initState, stateCell{row, col+3});
end

model.updCoordinateSet().get('pelvis_tx').setValue(initState, 0); % set pelvis_tx = 0
model.updCoordinateSet().get('pelvis_ty').setValue(initState, 0.93); % a little bit higher

end