% Compare OLD vs NEW fitness on existing opt_result_*.mat files.
% Only the effort term changed:
%   old: effort = C / pelvisX            (signed)
%   new: effort = C / max(pelvisX, 0.5)  (always positive)
% so  old_total = new_total - effort_new + effort_old.
clear; close all; clc
here = fileparts(mfilename('fullpath'));
repoRoot = fileparts(here);
cd(repoRoot);
addpath(genpath('assets'), genpath('model'), genpath('functions'))

projName = 'human0918';
simConfig.endTime = 10; simConfig.stepTime = 0.005;
simConfig.speed = 1.0;  simConfig.slope = 0;

files = { ...
    'results\opt_result_2026-07-15_14-34-18.mat', ...
    'results\opt_result_2026-07-29_17-51-41.mat', ...
    'results\opt_result_2026-08-03_07-16-07.mat', ...
    'results\opt_result_2026-09-29_17-54-54.mat'};

%% init pose from UN.sto
stoPath = 'Sparse Group LASSO Validation\UN.sto';
fid = fopen(stoPath, 'r');
for k = 1:6, fgetl(fid); end
headerLine = fgetl(fid);
dataLine = str2double(strsplit(fgetl(fid), '\t'));
fclose(fid);
colNames = strsplit(headerLine, '\t');
colMap = containers.Map('KeyType', 'char', 'ValueType', 'int32');
for i = 1:numel(colNames), colMap(colNames{i}) = i; end
dofNames = {'pelvis_tilt','pelvis_tx','pelvis_ty', ...
            'hip_flexion_r','knee_extension_r','ankle_dorsiflexion_r', ...
            'hip_flexion_l','knee_extension_l','ankle_dorsiflexion_l'};
nDof = numel(dofNames);
initPose = zeros(2*nDof, 1);
for j = 1:nDof
    initPose(2*j-1) = dataLine(colMap(dofNames{j}));
    initPose(2*j)   = dataLine(colMap([dofNames{j} '_u']));
end
modelStaticProp = read_muscle_static_prop(projName, simConfig, initPose, dofNames);
soCoordNames = {'hip_flexion_r','hip_flexion_l','knee_extension_r', ...
                'knee_extension_l','ankle_dorsiflexion_r','ankle_dorsiflexion_l'};
tauTarget = zeros(numel(soCoordNames),1);
for j = 1:numel(soCoordNames)
    tauTarget(j) = dataLine(colMap([soCoordNames{j} '.moment']));
end
tauLimit = zeros(numel(soCoordNames),1);
tauLimit(strcmp(soCoordNames,'knee_extension_r')) = dataLine(colMap('knee_r.torque'));
tauLimit(strcmp(soCoordNames,'knee_extension_l')) = dataLine(colMap('knee_l.torque'));
a_opt = iterative_static_optimization(projName, initPose, dofNames, tauTarget, tauLimit, soCoordNames);

fprintf('\n%-26s %8s %8s %9s %10s %10s %10s\n', ...
    'file', 'oldFit', 'x[m]', 'C[J/kg]', 'effOld', 'effNew', 'newFit');
for i = 1:numel(files)
    S = load(files{i}); R = S.result;
    [model, modelInfo] = init_infra(projName, modelStaticProp, a_opt);
    modelInfo.reset_record();
    [state, modelInfo] = reset_particle_state(model, modelInfo, a_opt);
    modelInfo.reflexParamMap = R.reflexParamMap;
    modelInfo.reflexTemplate = R.reflexTemplate;
    modelInfo.read_muscleReflex_array(R.bestPara);
    modelInfo = forward_simulation(model, modelInfo, state);

    [fit, parts] = measure_simResults(modelInfo);
    M = Measurements(modelInfo);
    x = M.distance;
    C = modelInfo.dy.totalMetabolic / M.mass;
    if abs(x) < 1e-9
        effOld = sign(C) * inf;
    else
        effOld = C / x;
    end
    if x >= 0.5
        effNew = C / x;        % unchanged (floor inactive)
    else
        effNew = C / 0.5;
    end
    oldTotal = fit - parts.effort + effOld;
    [~, nm] = fileparts(files{i});
    fprintf('%-26s %8.3f %8.3f %9.2f %10.3f %10.3f %10.3f\n', ...
        nm, oldTotal, x, C, effOld, effNew, fit);
end
