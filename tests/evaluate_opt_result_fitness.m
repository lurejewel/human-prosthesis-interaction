% -------------------------------------------------------------------------
% Name: evaluate_opt_result_fitness.m
% Description: Re-run ONE forward simulation for the best parameters stored in
%   an opt_result_*.mat and print the full fitness breakdown (completeness,
%   velocity, knee/ankle/GRF limits, effort).  Lightweight counterpart of
%   test_single_parameter.m (no plots, no 3D playback).
%
% Usage: set paraSourceFile below, then run this script.
% -------------------------------------------------------------------------
clear; close all; clc
here = fileparts(mfilename('fullpath'));
repoRoot = fileparts(here);
cd(repoRoot);
addpath(genpath('assets'), genpath('model'), genpath('functions'))

projName = 'human0918';
paraSourceFile = 'results\opt_result_2026-11-20_10-37-17.mat';
simConfig.endTime = 10; simConfig.stepTime = 0.005;
simConfig.speed = 1.0;  simConfig.slope = 0;

loaded = load(paraSourceFile);
R = loaded.result;
para = R.bestPara;
reflexParamMap = R.reflexParamMap;
reflexTemplate = R.reflexTemplate;
fprintf('loaded %s | bestFit=%.6g | lassoFile=%s\n', paraSourceFile, R.bestFit, R.lassoFile);

%% init pose from UN.sto (same as test_single_parameter)
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
fprintf('running static optimisation ...\n');
a_opt = iterative_static_optimization(projName, initPose, dofNames, tauTarget, tauLimit, soCoordNames);

[model, modelInfo] = init_infra(projName, modelStaticProp, a_opt);
modelInfo.reset_record();
[state, modelInfo] = reset_particle_state(model, modelInfo, a_opt);
modelInfo.reflexParamMap = reflexParamMap;
modelInfo.reflexTemplate = reflexTemplate;
modelInfo.read_muscleReflex_array(para);
fprintf('running forward simulation ...\n');
modelInfo = forward_simulation(model, modelInfo, state);

%% breakdown
M = Measurements(modelInfo);
nHS_R = numel(M.heelstrikeEventR); nHS_L = numel(M.heelstrikeEventL);
fprintf('\n--- state ---\n');
fprintf('pelvisX(final) = %.3f m | lastTime = %.2f s | npts = %d\n', M.distance, M.lastTime, M.npts);
fprintf('heel strikes R/L = %d/%d | totalMetabolic = %.1f J\n', nHS_R, nHS_L, modelInfo.dy.totalMetabolic);
[fit, parts] = measure_simResults(modelInfo);
fprintf('\n--- fitness terms (weighted) ---\n');
fprintf('completeness (w=1)    : %+.4f\n', parts.completeness);
fprintf('velocity     (w=100)  : %+.4f\n', parts.velocity);
fprintf('knee_limit   (w=0.01) : %+.4f\n', parts.knee_limit);
fprintf('ankle_limit  (w=0.01) : %+.4f\n', parts.ankle_limit);
fprintf('grf_limit    (w=10)   : %+.4f\n', parts.grf_limit);
fprintf('effort       (w=1)    : %+.4f\n', parts.effort);
fprintf('TOTAL                 : %+.4f\n', fit);
