function play_simulation_video(model, modelInfo, playbackSpeed, frameStep, showGRF)
% Name: play_simulation_video
% Description: Play a 3D visual playback of a completed forward simulation
%   using the OpenSim visualiser.  Reads state history and time series
%   directly from modelInfo (no .sto file needed).
%   Optionally displays ground reaction force (GRF) arrows by adding
%   PrescribedForce objects to the model.
%
% Input:
%   model         - org.opensim.modeling.Model (will have visualiser enabled)
%   modelInfo     - ModelInfo object with dy.labelHistory and st.simInfo.timeSeries
%                   (and dy.grf fields if showGRF = true)
%   playbackSpeed - (optional) speed multiplier (>1 = faster, <1 = slower);
%                   default 1.0 (real-time)
%   frameStep     - (optional) render every N-th frame; default 1 (every frame)
%   showGRF       - (optional) logical; if true, add GRF arrows at feet.
%                   Default false.
%
% Usage:
%   play_simulation_video(model, modelInfo)
%   play_simulation_video(model, modelInfo, 2.0, 5)       % 2x speed, every 5th frame
%   play_simulation_video(model, modelInfo, 1.0, 1, true)  % real-time with GRF arrows
%
% Note: showGRF = true permanently adds PrescribedForce objects to the model.
%   Subsequent forward simulations that use this model should reinitialise it.
%
% Dependencies: OpenSim MATLAB API

if nargin < 3 || isempty(playbackSpeed)
    playbackSpeed = 1.0;
end
if nargin < 4 || isempty(frameStep)
    frameStep = 1;
end
if nargin < 5 || isempty(showGRF)
    showGRF = false;
end

%% ---- extract simulation data ----
labelHist = modelInfo.dy.labelHistory;          % label order
permLabelToInternal = modelInfo.st.model.permLabelToInternal;  % label→internal permutation
timeSeries = modelInfo.st.simInfo.timeSeries;
nValid = find(all(~isnan(labelHist)), 1, 'last');
nStates = size(labelHist, 1);

if nValid < 2
    warning('play_simulation_video: not enough valid frames (nValid = %d).', nValid);
    return
end

tVec = timeSeries(1:nValid);
frames = 1:frameStep:nValid;
nFrames = numel(frames);

fprintf('Playing 3D simulation video (%d frames, %.1fx speed, step=%d)...\n', ...
    nFrames, playbackSpeed, frameStep);

%% ---- set up visualiser ----
model.setUseVisualizer(true);
stateVis = model.initSystem();

% Set first-frame kinematics so initial body positions are correct
stateVis.setTime(tVec(1));
Y0 = stateVis.updY();
% Permute from label order → SimTK internal order for the visualiser
labelFrame = labelHist(:, 1);
internalFrame = labelFrame(permLabelToInternal);
for i = 0:nStates - 1
    Y0.set(i, internalFrame(i + 1));
end
model.realizePosition(stateVis);

%% ---- playback loop ----
nextGrfTime = 0;
grfInterval = 0.1;  % seconds between GRF snapshots
grfScale = 5e-4;

if showGRF
    vis = model.updVisualizer().updSimbodyVisualizer();
    bodyCalcnR = model.getBodySet().get('calcn_r');
    bodyCalcnL = model.getBodySet().get('calcn_l');
    grdIdx = 0;
end

for k = 1:nFrames
    fi = frames(k);
    tNow = tVec(fi);

    stateVis.setTime(tNow);
    Y = stateVis.updY();
    % Permute from label order → SimTK internal order for the visualiser
    labelFrame = labelHist(:, fi);
    internalFrame = labelFrame(permLabelToInternal);
    for i = 0:nStates - 1
        Y.set(i, internalFrame(i + 1));
    end

    model.realizeDynamics(stateVis);

    % Draw GRF snapshot every grfInterval seconds
    if showGRF && tNow >= nextGrfTime - 1e-9
        alpha = tNow / tVec(end);  % 0→1 over the simulation

        % Right foot — warm gradient: red → orange → yellow
        colR = org.opensim.modeling.Vec3(1, alpha, 0);
        posR = bodyCalcnR.getPositionInGround(stateVis);
        lineR = org.opensim.modeling.DecorativeLine();
        lineR.setColor(colR);
        lineR.setLineThickness(2);
        vis.addRubberBandLine(grdIdx, posR, grdIdx, ...
            org.opensim.modeling.Vec3(posR.get(0) + modelInfo.dy.grf.fxr(fi) * grfScale, ...
                 posR.get(1) + modelInfo.dy.grf.fyr(fi) * grfScale, posR.get(2)), lineR);

        % Left foot — cool gradient: dark blue → teal → green
        colL = org.opensim.modeling.Vec3(alpha * 0.3, alpha, 0.6 - alpha * 0.2);
        posL = bodyCalcnL.getPositionInGround(stateVis);
        lineL = org.opensim.modeling.DecorativeLine();
        lineL.setColor(colL);
        lineL.setLineThickness(2);
        vis.addRubberBandLine(grdIdx, posL, grdIdx, ...
            org.opensim.modeling.Vec3(posL.get(0) + modelInfo.dy.grf.fxl(fi) * grfScale, ...
                 posL.get(1) + modelInfo.dy.grf.fyl(fi) * grfScale, posL.get(2)), lineL);

        nextGrfTime = tNow + grfInterval;
    end

    model.updVisualizer().show(stateVis);

    if k < nFrames
        dt = (tVec(frames(k + 1)) - tNow) / playbackSpeed;
        pause(dt);
    end
end

fprintf('3D playback finished.\n');

end
