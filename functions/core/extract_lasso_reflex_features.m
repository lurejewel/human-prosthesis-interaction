function A = extract_lasso_reflex_features(modelInfo, side, frameIndex)
% Name: extract_lasso_reflex_features
% Description: Build the 1×22 feature row vector for one leg side.
%   Feature order (22 elements):
%     1-7   : normalized fiber length  of 7 same-side muscles
%     8-14  : normalized MTU force     of 7 same-side muscles
%     15    : pelvis tilt angle
%     16    : pelvis tilt angular velocity
%     17    : same-side hip angle
%     18    : same-side hip angular velocity
%     19    : same-side knee angle
%     20    : same-side knee angular velocity
%     21    : same-side ankle angle
%     22    : same-side ankle angular velocity
%
%   All joint angles/velocities reuse the same stateMap indexing and sign
%   conventions as the original cal_muscle_excitation / compute_leg_excitation.
%
% Input:
%   modelInfo  – ModelInfo handle object
%   side       – 'right' or 'left'
%   frameIndex – current frame index (1-based)

validatestring(side, {'right', 'left'});

nMusPerSide = 7;
if strcmp(side, 'right')
    muscleIdx = 1:nMusPerSide;
    legSuffix = 'r';
else
    muscleIdx = (nMusPerSide+1):(2*nMusPerSide);
    legSuffix = 'l';
end

% ---- muscle-level features (indices 1-14) ----
fiberLen = modelInfo.dy.muscle.lCEN(muscleIdx, frameIndex);
mtuForce = modelInfo.dy.muscle.fATN(muscleIdx, frameIndex);

% ---- joint kinematic features (indices 15-22) ----
% Reuse the exact same indexing from the old cal_muscle_excitation
map = modelInfo.st.model.map;

if frameIndex == 1
    initPose = modelInfo.st.model.initPose;
    pelvisTilt  = initPose(map('pelvis_tilt/value'));
    pelvisTiltV = initPose(map('pelvis_tilt/speed'));
else
    sh = modelInfo.dy.stateHistory;
    pelvisTilt  = sh(map('pelvis_tilt/value'), frameIndex - 1);
    pelvisTiltV = sh(map('pelvis_tilt/speed'), frameIndex - 1);
end

if frameIndex == 1
    hipAng  = initPose(map(['hip_flexion_' legSuffix '/value']));
    hipVel  = initPose(map(['hip_flexion_' legSuffix '/speed']));
    kneeAng = initPose(map(['knee_extension_' legSuffix '/value']));
    kneeVel = initPose(map(['knee_extension_' legSuffix '/speed']));
    ankleAng = initPose(map(['ankle_dorsiflexion_' legSuffix '/value']));
    ankleVel = initPose(map(['ankle_dorsiflexion_' legSuffix '/speed']));
else
    sh = modelInfo.dy.stateHistory;
    hipAng  = sh(map(['hip_flexion_' legSuffix '/value']), frameIndex - 1);
    hipVel  = sh(map(['hip_flexion_' legSuffix '/speed']), frameIndex - 1);
    kneeAng = sh(map(['knee_extension_' legSuffix '/value']), frameIndex - 1);
    kneeVel = sh(map(['knee_extension_' legSuffix '/speed']), frameIndex - 1);
    ankleAng = sh(map(['ankle_dorsiflexion_' legSuffix '/value']), frameIndex - 1);
    ankleVel = sh(map(['ankle_dorsiflexion_' legSuffix '/speed']), frameIndex - 1);
end

% ---- assemble 1×22 row vector ----
A = [fiberLen(:)', mtuForce(:)', ...
     pelvisTilt, pelvisTiltV, ...
     hipAng, hipVel, ...
     kneeAng, kneeVel, ...
     ankleAng, ankleVel];

% ---- validation ----
assert(isequal(size(A), [1, 22]), 'Feature vector must be 1×22.');
assert(all(isfinite(A)), 'Feature vector contains non-finite values.');

end
