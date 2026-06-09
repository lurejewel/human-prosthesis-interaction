% -------------------------------------------------------------------------
% Name: test_single_parameter.m
% Author(s): Jin Wei, Peking U. wjin24@stu.pku.edu.cn
% Description: Test a single set of muscle-reflex parameters by running one
%   forward-dynamics simulation and visualizing the resulting gait.
%   Intended for rapidly evaluating the output of a previous CMA-ES
%   optimization without re-running the whole optimization loop.
%
% Usage:
%   1. Set `paraSourceFile` to a previous opt_result_*.mat file.
%   2. Run the script. It will load the best parameters, execute one
%      forward simulation, print the fitness breakdown, and display
%      multi-panel diagnostic plots.
% -------------------------------------------------------------------------

clear all; close all; clc
addpath(genpath('assets\'), genpath('model\'), genpath('functions\'))

%% ---------- user configuration ----------
projName       = 'coupled_human-prosthesis_model';
paraSourceFile = 'results\opt_result_2026-06-02_00-10-00.mat';  % <-- change to your result file

simConfig.endTime = 10;
simConfig.stepTime = 0.005;
simConfig.speed   = 1.0;
simConfig.slope   = 0;
showVideo         = true;  % set to false to skip 3D visual playback
% -----------------------------------------

%% load the parameter set to test
loaded = load(paraSourceFile);
if isfield(loaded, 'result')
    para = loaded.result.bestPara;
    fprintf('Loaded bestPara from: %s  (bestFit = %.6g, %d generations)\n', ...
        paraSourceFile, loaded.result.bestFit, loaded.result.generations);
else
    para = loaded.bestPara;  % legacy format
    fprintf('Loaded bestPara from: %s (legacy format)\n', paraSourceFile);
end

%% initialise model infrastructure (same pipeline as the main demo)
initPose = [-0.105763, 0, 0.900237, 0.439316, 0.198813, -0.393922, -1.03755, 0.104714, -0.348473, ...
            -0.0895989, 1.0757, 0.1543, -1.35971, 3.34368, 0.267883, -3.15281, 0.840122, 1.26642];
modelStaticProp = read_muscle_static_prop(projName, simConfig, initPose);
[model, modelInfo, ~] = init_infra(projName, modelStaticProp);
model.setUseVisualizer(true);  % enable 3D visualisation

%% set parameters and reset state
modelInfo.reset_record();
[state, modelInfo] = reset_particle_state(model, modelInfo);
modelInfo.read_muscleReflex_array(para, zeros(size(para)));  % arz not needed for a single run

%% run forward simulation
modelInfo = forward_simulation(model, modelInfo, state);

%% fitness evaluation
fit = measure_simResults(modelInfo);
fprintf('\n========== Fitness Summary ==========\n');
fprintf('Total fitness              : %.6g\n', fit);
fprintf('Last simulation time       : %.3f s (of %.0f s)\n', ...
    modelInfo.dy.lastTime, simConfig.endTime);
fprintf('Distance travelled         : %.3f m\n', ...
    modelInfo.dy.stateHistory(modelInfo.st.model.map('pelvis_tx/value'), ...
    find(all(~isnan(modelInfo.dy.stateHistory)), 1, 'last')));
fprintf('======================================\n\n');

%% prepare data for visualisation
stateMap   = modelInfo.st.model.map;
stateHist  = modelInfo.dy.stateHistory;
nValid     = find(all(~isnan(stateHist)), 1, 'last');  % last fully-recorded frame
tVec       = modelInfo.st.simInfo.timeSeries(1:nValid);

% ---- joint angles (deg) ----
hipR_ang   = rad2deg(stateHist(stateMap('hip_flexion_r/value'),      1:nValid));
hipL_ang   = rad2deg(stateHist(stateMap('hip_flexion_l/value'),      1:nValid));
kneeR_ang  = rad2deg(stateHist(stateMap('knee_flexion_r/value'),     1:nValid));
kneeL_ang  = rad2deg(stateHist(stateMap('knee_flexion_l/value'),     1:nValid));
ankleR_ang = rad2deg(stateHist(stateMap('ankle_dorsiflexion_r/value'),1:nValid));
ankleL_ang = rad2deg(stateHist(stateMap('ankle_dorsiflexion_l/value'),1:nValid));

% ---- pelvis trajectory ----
pelvis_tx  = stateHist(stateMap('pelvis_tx/value'), 1:nValid);
pelvis_ty  = stateHist(stateMap('pelvis_ty/value'), 1:nValid);

% ---- gait phases ----
phaseR = modelInfo.dy.phase.r(1:nValid);
phaseL = modelInfo.dy.phase.l(1:nValid);

% ---- muscle excitations ----
delay = modelInfo.st.muscle.delay;
exc = modelInfo.dy.muscle.exc(:, 1:nValid);
muscleNames = modelInfo.st.muscle.names;

% ---- ground reaction forces ----
grfYR = modelInfo.dy.grf.fyr(1:nValid);
grfYL = modelInfo.dy.grf.fyl(1:nValid);

% ---- walking speed (m/s, from pelvis_tx gradient) ----
walkSpeed = gradient(pelvis_tx, simConfig.stepTime);

%% =====================  VISUALISATION  =====================

figure('Name', 'Gait Simulation Results', 'Position', [50, 50, 1400, 900]);

% -------- Panel A: Joint Angles --------
subplot(3, 3, 1)
plot(tVec, hipR_ang, 'LineWidth', 1.5); hold on
plot(tVec, hipL_ang, '--', 'LineWidth', 1.5);
ylabel('Hip flexion (deg)'); xlabel('Time (s)');
legend('Right', 'Left', 'Location', 'best'); grid on; title('Hip Angles');

subplot(3, 3, 2)
plot(tVec, kneeR_ang, 'LineWidth', 1.5); hold on
plot(tVec, kneeL_ang, '--', 'LineWidth', 1.5);
ylabel('Knee flexion (deg)'); xlabel('Time (s)');
legend('Right', 'Left', 'Location', 'best'); grid on; title('Knee Angles');

subplot(3, 3, 3)
plot(tVec, ankleR_ang, 'LineWidth', 1.5); hold on
plot(tVec, ankleL_ang, '--', 'LineWidth', 1.5);
ylabel('Ankle dorsi. (deg)'); xlabel('Time (s)');
legend('Right', 'Left', 'Location', 'best'); grid on; title('Ankle Angles');

% -------- Panel B: Pelvis Trajectory --------
subplot(3, 3, 4)
plot(pelvis_tx, pelvis_ty, 'LineWidth', 1.5);
xlabel('Pelvis X (m)'); ylabel('Pelvis Y (m)');
grid on; title('Pelvis Trajectory (sagittal)'); axis equal;

% -------- Panel C: Gait Phases --------
subplot(3, 3, 5)
yyaxis left
stairs(tVec, phaseR, 'LineWidth', 1.5);
ylim([-0.2, 4.2]); ylabel('Right leg phase');
yyaxis right
stairs(tVec, phaseL, '--', 'LineWidth', 1.5);
ylim([-0.2, 4.2]); ylabel('Left leg phase');
xlabel('Time (s)'); grid on; title('Gait Phases (0-ES,1-LS,2-LO,3-SW,4-LD)');
legend('Right', 'Left', 'Location', 'best');

% -------- Panel D: Muscle Excitations (4 representative muscles) --------
subplot(3, 3, 6)
musclesToPlot = {'glut_max_r', 'vasti_r', 'gastroc_r', 'tibia_r'};
colors = lines(numel(musclesToPlot));
for j = 1:numel(musclesToPlot)
    idx = find(strcmp(muscleNames, musclesToPlot{j}), 1);
    if ~isempty(idx)
        plot(tVec, exc(idx, :), 'Color', colors(j, :), 'LineWidth', 1.2); hold on
    end
end
ylabel('Excitation'); xlabel('Time (s)');
legend(musclesToPlot, 'Location', 'best', 'Interpreter', 'none');
grid on; title('Muscle Excitations (right leg)');

% -------- Panel E: Ground Reaction Forces --------
subplot(3, 3, 7)
plot(tVec, grfYR, 'LineWidth', 1.5); hold on
plot(tVec, grfYL, '--', 'LineWidth', 1.5);
ylabel('Normal GRF (N)'); xlabel('Time (s)');
legend('Right', 'Left', 'Location', 'best'); grid on; title('Vertical Ground Reaction Force');

% -------- Panel F: Walking Speed --------
subplot(3, 3, 8)
plot(tVec, walkSpeed, 'LineWidth', 1.5); hold on
yline(simConfig.speed, 'r--', 'LineWidth', 1.2);
ylabel('Speed (m/s)'); xlabel('Time (s)');
legend('Actual', 'Target', 'Location', 'best'); grid on; title('Walking Speed');

% -------- Panel G: All muscle excitations (heatmap) --------
subplot(3, 3, 9)
imagesc(tVec, 1:numel(muscleNames), exc);
xlabel('Time (s)'); ylabel('Muscle index');
yticks(1:numel(muscleNames)); yticklabels(muscleNames);
set(gca, 'YTickLabel', get(gca, 'YTickLabel'), 'FontSize', 7);
colorbar; title('Muscle Excitations (heatmap)'); colormap(jet);

sgtitle(sprintf('Gait Simulation  |  Fitness = %.4g  |  Distance = %.2f m  |  Time = %.2f s', ...
    fit, pelvis_tx(end), modelInfo.dy.lastTime), 'FontSize', 12, 'FontWeight', 'bold');

endOfSim = modelInfo.dy.lastTime;
if endOfSim < simConfig.endTime
    fprintf('Note: model fell at t = %.2f s (pelvis height < 0.6 m).\n', endOfSim);
end

%% ================  3D VISUAL PLAYBACK  ================

if showVideo
    playbackSpeed = 1.0;   % >1 = faster, <1 = slower, 1 = real-time
    frameStep     = 1;     % render every N-th frame (1 = every frame)

    stateVis = model.initSystem();
    nStates  = size(stateHist, 1);

    fprintf('Starting 3D playback (%d frames, %.1fx speed)...\n', nValid, playbackSpeed);
    for fi = 1 : frameStep : nValid
        stateVis.setTime(tVec(fi));
        Y = stateVis.updY();
        for i = 0 : nStates - 1
            Y.set(i, stateHist(i + 1, fi));
        end
        model.realizeDynamics(stateVis);
        model.updVisualizer().show(stateVis);

        if fi + frameStep <= nValid
            dt = (tVec(min(fi + frameStep, nValid)) - tVec(fi)) / playbackSpeed;
            pause(dt);
        end
    end
    fprintf('3D playback finished.\n');
end
