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
cd(fileparts(mfilename('fullpath')));
cd('..');
addpath(genpath('assets'), genpath('model'), genpath('functions'))

%% ---------- user configuration ----------
projName       = 'human0918';  % <-- change to your model name (without .osim)
paraSourceFile = 'results\opt_result_2026-07-06_20-35-01.mat';  % <-- change to your result file

simConfig.endTime = 10;
simConfig.stepTime = 0.005;
simConfig.speed   = 1.0;
simConfig.slope   = 0;
showVideo         = 1;  % set to false to skip 3D visual playback
% -----------------------------------------

%% load the parameter set to test
loaded = load(paraSourceFile);
if isfield(loaded, 'result')
    para = loaded.result.bestPara;
    % ---- LASSO controller metadata (needed for matrix-multiplication form) ----
    if isfield(loaded.result, 'reflexParamMap') && isfield(loaded.result, 'reflexTemplate')
        reflexParamMap = loaded.result.reflexParamMap;
        reflexTemplate = loaded.result.reflexTemplate;
    else
        error('Result file does not contain reflexParamMap/reflexTemplate. ' + ...
              'Re-run the optimisation with the current demo script.');
    end
    fprintf('Loaded bestPara from: %s  (bestFit = %.6g, %d generations)\n', ...
        paraSourceFile, loaded.result.bestFit, loaded.result.generations);
else
    para = loaded.bestPara;  % legacy format
    fprintf('Loaded bestPara from: %s (legacy format)\n', paraSourceFile);
    error('Legacy format does not contain reflexParamMap/reflexTemplate. ' + ...
          'Re-run the optimisation with the current demo script.');
end

%% read initial pose and joint moments from UN.sto
stoPath = 'Sparse Group LASSO Validation\UN.sto';
assert(isfile(stoPath), 'UN.sto not found: %s', stoPath);

fid = fopen(stoPath, 'r');
for k = 1:6, fgetl(fid); end
headerLine = fgetl(fid);
dataLine   = str2double(strsplit(fgetl(fid), '\t'));
fclose(fid);
colNames = strsplit(headerLine, '\t');
colMap = containers.Map('KeyType', 'char', 'ValueType', 'int32');
for i = 1:numel(colNames), colMap(colNames{i}) = i; end

% 9 DOF (Q then U → 18 elements)
dofNames = {'pelvis_tilt','pelvis_tx','pelvis_ty', ...
            'hip_flexion_r','knee_extension_r','ankle_dorsiflexion_r', ...
            'hip_flexion_l','knee_extension_l','ankle_dorsiflexion_l'};
nDof = numel(dofNames);
initPose = zeros(2 * nDof, 1);
for j = 1:nDof
    initPose(2*j - 1) = dataLine(colMap(dofNames{j}));
    initPose(2*j)     = dataLine(colMap([dofNames{j} '_u']));
end

%% initialise model infrastructure (same pipeline as the main demo)
modelStaticProp = read_muscle_static_prop(projName, simConfig, initPose, dofNames);

% ---- iterative static optimisation for first-frame muscle activations ----
fprintf('Running iterative static optimisation ...\n');
soCoordNames = {'hip_flexion_r','hip_flexion_l', ...
                'knee_extension_r','knee_extension_l', ...
                'ankle_dorsiflexion_r','ankle_dorsiflexion_l'};
nSoCoords = numel(soCoordNames);
tauTarget = zeros(nSoCoords, 1);
for j = 1:nSoCoords
    tauTarget(j) = dataLine(colMap([soCoordNames{j} '.moment']));
end
tauLimit = zeros(nSoCoords, 1);
tauLimit(strcmp(soCoordNames, 'knee_extension_r')) = dataLine(colMap('knee_r.torque'));
tauLimit(strcmp(soCoordNames, 'knee_extension_l')) = dataLine(colMap('knee_l.torque'));
a_opt = iterative_static_optimization(projName, initPose, dofNames, ...
    tauTarget, tauLimit, soCoordNames);

[model, modelInfo] = init_infra(projName, modelStaticProp, a_opt);

%% set parameters and reset state
modelInfo.reset_record();
[state, modelInfo] = reset_particle_state(model, modelInfo, a_opt);

% ---- wire LASSO controller metadata BEFORE unpacking parameters ----
modelInfo.reflexParamMap = reflexParamMap;
modelInfo.reflexTemplate = reflexTemplate;
modelInfo.read_muscleReflex_array(para);

%% run forward simulation
modelInfo = forward_simulation(model, modelInfo, state);

%% fitness evaluation
fit = measure_simResults(modelInfo);
fprintf('\n========== Fitness Summary ==========\n');
fprintf('Total fitness              : %.6g\n', fit);
fprintf('Last simulation time       : %.3f s (of %.0f s)\n', ...
    modelInfo.dy.lastTime, simConfig.endTime);
fprintf('Distance travelled         : %.3f m\n', ...
    modelInfo.dy.labelHistory(modelInfo.st.model.map('pelvis_tx/value'), ...
    find(all(~isnan(modelInfo.dy.labelHistory)), 1, 'last')));
fprintf('======================================\n\n');

%% prepare data for visualisation
stateMap   = modelInfo.st.model.map;
labelHist  = modelInfo.dy.labelHistory;   % label order (matches stateMap indices)
nValid     = find(all(~isnan(labelHist)), 1, 'last');
tVec       = modelInfo.st.simInfo.timeSeries(1:nValid);

% ---- joint angles (deg) ----
hipR_ang   = rad2deg(labelHist(stateMap('hip_flexion_r/value'),      1:nValid));
hipL_ang   = rad2deg(labelHist(stateMap('hip_flexion_l/value'),      1:nValid));
kneeR_ang  = rad2deg(labelHist(stateMap('knee_extension_r/value'),     1:nValid));
kneeL_ang  = rad2deg(labelHist(stateMap('knee_extension_l/value'),     1:nValid));
ankleR_ang = rad2deg(labelHist(stateMap('ankle_dorsiflexion_r/value'),1:nValid));
ankleL_ang = rad2deg(labelHist(stateMap('ankle_dorsiflexion_l/value'),1:nValid));

% ---- pelvis trajectory ----
pelvis_tx  = labelHist(stateMap('pelvis_tx/value'), 1:nValid);
pelvis_ty  = labelHist(stateMap('pelvis_ty/value'), 1:nValid);

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

%% compute net joint torques via muscle moment summation
% Net joint torque ≈ Σ (moment_arm_i × tendon_force_i) for all muscles.
% Uses the forward-simulation model directly — no separate ID model needed.
fprintf('Computing net joint torques via muscle moment summation (%d frames)...\n', nValid);

% ---- coordinate handles (retrieved once for efficiency) ----
coord_hipR   = model.getCoordinateSet().get('hip_flexion_r');
coord_kneeR  = model.getCoordinateSet().get('knee_extension_r');
coord_ankleR = model.getCoordinateSet().get('ankle_dorsiflexion_r');
coord_hipL   = model.getCoordinateSet().get('hip_flexion_l');
coord_kneeL  = model.getCoordinateSet().get('knee_extension_l');
coord_ankleL = model.getCoordinateSet().get('ankle_dorsiflexion_l');
coords = {coord_hipR, coord_kneeR, coord_ankleR, ...
          coord_hipL, coord_kneeL, coord_ankleL};

% ---- muscle handles ----
nMus = numel(muscleNames);
nQ = model.getNumCoordinates();  % same as state.getNQ()
muscleRefs = cell(1, nMus);
for i = 1 : nMus
    muscleRefs{i} = model.getMuscles().get(i - 1);
end

% ---- state for moment-arm queries (fresh copy, kinematics set per frame) ----
% labelHistory is in label order; permute to internal order for state.updY().set()
permLabelToInternal = modelInfo.st.model.permLabelToInternal;
stateMA = model.initSystem();
nStates = stateMA.getNQ() + stateMA.getNU();

% ---- preallocate ----
torque_hipR   = nan(1, nValid);
torque_kneeR  = nan(1, nValid);
torque_ankleR = nan(1, nValid);
torque_hipL   = nan(1, nValid);
torque_kneeL  = nan(1, nValid);
torque_ankleL = nan(1, nValid);

for fi = 1 : nValid
    % Permute from label order to internal order for state.updY()
    labelFrame = labelHist(:, fi);
    internalFrame = labelFrame(permLabelToInternal);
    Y = stateMA.updY();
    for i = 0 : nStates - 1
        Y.set(i, internalFrame(i + 1));
    end
    stateMA.setTime(tVec(fi));

    % Realize to Velocity stage (required for path-based muscle moment arms)
    model.realizeVelocity(stateMA);

    % ---- accumulate muscle moments ----
    tau = zeros(6, 1);
    for m = 1 : nMus
        fT = modelInfo.dy.muscle.fATN(m, fi) * muscleRefs{m}.getMaxIsometricForce;
        if fT == 0 || isnan(fT)
            continue;
        end
        for c = 1 : 6
            ma = muscleRefs{m}.computeMomentArm(stateMA, coords{c});
            tau(c) = tau(c) + ma * fT;
        end
    end

    % ---- add knee coordinate limit forces ----
    tau(2) = tau(2) + modelInfo.dy.limitForce.kneeR(fi);
    tau(5) = tau(5) + modelInfo.dy.limitForce.kneeL(fi);

    torque_hipR(fi)   = tau(1);
    torque_kneeR(fi)  = tau(2);
    torque_ankleR(fi) = tau(3);
    torque_hipL(fi)   = tau(4);
    torque_kneeL(fi)  = tau(5);
    torque_ankleL(fi) = tau(6);
end
fprintf('Muscle moment summation done.\n'); % 这里需要验证

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

% -------- Panel D: Muscle Excitations (6 representative muscles) --------
subplot(3, 3, 6)
musclesToPlot = {'bifemsh_r', 'glut_max_r', 'rect_fem_r', 'vasti_r', 'gastroc_r', 'tib_ant_r'};
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

% -------- Panel G: Joint Torques (right leg) --------
subplot(3, 3, 9)
plot(tVec, torque_hipR, 'LineWidth', 1.5); hold on
plot(tVec, torque_kneeR, 'LineWidth', 1.5);
plot(tVec, torque_ankleR, 'LineWidth', 1.5);
ylabel('Joint moment (Nm)'); xlabel('Time (s)');
legend('Hip', 'Knee', 'Ankle', 'Location', 'best');
grid on; title('Joint Torques (right leg)');

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

    model.setUseVisualizer(true);
    stateVis = model.initSystem();
    % labelHistory is in label order; permute to internal order for state.updY()
    permInternalToLabel = modelInfo.st.model.permInternalToLabel; %permLabelToInternal;
    nStates  = size(labelHist, 1);

    fprintf('Starting 3D playback (%d frames, %.1fx speed)...\n', nValid, playbackSpeed);
    for fi = 1 : frameStep : nValid
        stateVis.setTime(tVec(fi));
        Y = stateVis.updY();
        internalFrame = labelHist(:, fi);
        internalFrame = internalFrame(permInternalToLabel);
        for i = 0 : nStates - 1
            Y.set(i, internalFrame(i + 1));
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

%% ================  GAIT-CYCLE-NORMALISED PLOTS  ================

% ---- normalisation constants (body weight / mass) ----
bodyMass   = modelStaticProp.model.totalMass;     % kg
bodyWeight = bodyMass * 9.80665;                  % N

% ---- detect heel-strike events (phase 4 -> 0 transition) ----
hsRight = find(phaseR(1:end-1) == 4 & phaseR(2:end) == 0) + 1;
if length(hsRight) < 2
    warning('Fewer than 2 right heel strikes detected; skipping gait-cycle plots.');
else
    fprintf('Detected %d right heel-strike events. Normalising gait cycles...\n', length(hsRight));

    nCycles = length(hsRight) - 1;
    nPts    = 101;  % 0..100 %% gait cycle
    gcPct   = linspace(0, 100, nPts);

    % preallocate: nCycles × nPts
    gc_hipAngR   = nan(nCycles, nPts);
    gc_kneeAngR  = nan(nCycles, nPts);
    gc_ankleAngR = nan(nCycles, nPts);
    gc_hipTauR   = nan(nCycles, nPts);
    gc_kneeTauR  = nan(nCycles, nPts);
    gc_ankleTauR = nan(nCycles, nPts);
    % muscle forces (all 14) per cycle, raw (N)
    gc_musForce  = nan(nCycles, nPts, numel(muscleNames));

    for c = 1 : nCycles
        idxStart = hsRight(c);
        idxEnd   = hsRight(c + 1);
        nFrames  = idxEnd - idxStart + 1;
        framePct = linspace(0, 100, nFrames);

        gc_hipAngR(c, :)   = interp1(framePct, hipR_ang(idxStart:idxEnd),   gcPct, 'pchip');
        gc_kneeAngR(c, :)  = interp1(framePct, kneeR_ang(idxStart:idxEnd),  gcPct, 'pchip');
        gc_ankleAngR(c, :) = interp1(framePct, ankleR_ang(idxStart:idxEnd), gcPct, 'pchip');
        % normalise torques by body mass (Nm/kg)
        gc_hipTauR(c, :)   = interp1(framePct, torque_hipR(idxStart:idxEnd),   gcPct, 'pchip') / bodyMass;
        gc_kneeTauR(c, :)  = interp1(framePct, torque_kneeR(idxStart:idxEnd),  gcPct, 'pchip') / bodyMass;
        gc_ankleTauR(c, :) = interp1(framePct, torque_ankleR(idxStart:idxEnd), gcPct, 'pchip') / bodyMass;

        for m = 1 : numel(muscleNames)
            fATN_cycle = modelInfo.dy.muscle.fATN(m, idxStart:idxEnd);
            % normalise by body weight (dimensionless, F / BW)
            gc_musForce(c, :, m) = interp1(framePct, fATN_cycle, gcPct, 'pchip') / bodyWeight;
        end
    end

    % ---- define muscle groups (individual muscles, NOT summed) ----
    grpNames = {'Hip muscles', 'Knee muscles', 'Ankle muscles'};
    % each group: list of muscle names (right-leg only for the normalised plot)
    grpMuscles = { ...
        {'hamstrings_r','bifemsh_r','glut_max_r','iliopsoas_r'}, ...
        {'rect_fem_r','vasti_r'}, ...
        {'gastroc_r','soleus_r','tib_ant_r'} };

    % ---- plot ----
    figure('Name', 'Gait-Cycle-Normalised Kinetics & Kinematics', ...
           'Position', [100, 100, 1400, 900]);

    rowLabels = {'Joint Angles (deg)', 'Joint Torques (Nm/kg)', 'Muscle Forces (F / BW)'};
    colLabels = {'Hip', 'Knee', 'Ankle'};

    angData = {gc_hipAngR, gc_kneeAngR, gc_ankleAngR};
    tauData = {gc_hipTauR, gc_kneeTauR, gc_ankleTauR};

    for col = 1 : 3
        % Row 1: angles
        subplot(3, 3, col);
        plot_gc_shaded(gcPct, angData{col});
        if col == 1, ylabel(rowLabels{1}); end
        title(colLabels{col});

        % Row 2: torques
        subplot(3, 3, 3 + col);
        plot_gc_shaded(gcPct, tauData{col});
        if col == 1, ylabel(rowLabels{2}); end

        % Row 3: individual muscle forces (one curve per muscle)
        subplot(3, 3, 6 + col);
        musList = grpMuscles{col};
        nCurves = numel(musList);
        colors  = lines(nCurves);
        legStr  = cell(1, nCurves);
        for j = 1 : nCurves
            mIdx = find(strcmp(muscleNames, musList{j}), 1);
            if ~isempty(mIdx)
                plot_gc_shaded(gcPct, gc_musForce(:, :, mIdx), colors(j, :));
                hold on;
                % build legend: strip '_r' suffix
                nm = musList{j};
                if nm(end) == 'r' || (length(nm) > 1 && strcmp(nm(end-1:end), '_r'))
                    nm = nm(1:end-2);
                end
                legStr{j} = nm;
            end
        end
        hold off;
        if col == 1, ylabel(rowLabels{3}); end
        xlabel('Gait cycle (%)');
        legend(legStr, 'Location', 'best', 'Interpreter', 'none');
        title(sprintf('%s (%s)', colLabels{col}, grpNames{col}));
    end

    sgtitle(sprintf('Gait-Cycle-Normalised  |  %d cycles  |  Right leg  |  BW = %.0f N', ...
            nCycles, bodyWeight), 'FontSize', 12, 'FontWeight', 'bold');
end

%% ====================  LOCAL HELPERS  ====================
function plot_gc_shaded(x, data, color)
% plot mean ± shaded std for gait-cycle data (nCycles × nPts)
%   color (optional) – line and fill colour; default grey
if nargin < 3 || isempty(color)
    color = [0.5 0.5 0.5];
end
mu = mean(data, 1, 'omitnan');
sd = std(data, 0, 1, 'omitnan');
x2 = [x, fliplr(x)];
fill(x2, [mu - sd, fliplr(mu + sd)], color, ...
    'FaceAlpha', 0.15, 'EdgeColor', 'none');
hold on;
plot(x, mu, 'Color', color, 'LineWidth', 1.8);
hold off;
grid on;
end
