function a_opt = iterative_static_optimization(projName, initPose, dofNames, ...
    tauTarget, tauLimit, soCoordNames)
% Name: iterative_static_optimization
% Description: Compute initial muscle activations for a static pose by
%   iteratively solving a convex QP that matches muscle-generated joint
%   moments to externally-supplied target moments.
%
%   The muscle force model accounts for force-length (f_l), force-velocity
%   (f_v), and pennation (cos α) at the current state.  Passive muscle
%   forces and knee coordinate-limit torques are subtracted from the net
%   joint moment before solving the QP.
%
% Input:
%   projName      - model base name (e.g. 'human0714')
%   initPose      - 18×1 double, initial kinematics [Q; U] for all 9 DOF
%   dofNames      - 9×1 cell of char, coordinate names matching initPose
%   tauTarget     - nCoords×1 double, net joint moments from UN.sto
%   tauLimit      - nCoords×1 double, knee limit-force torques (0 elsewhere)
%   soCoordNames  - nCoords×1 cell of char, coordinates to match
%
% Output:
%   a_opt         - nMusTotal×1 double, converged muscle activations
%
% Algorithm:
%   a ← 0.05;  R ← moment arms at a=0.05
%   repeat:
%       setActivation(a) → equilibrate → realise
%       Fcap ← Fopt * f_l * f_v * cos(α)
%       A_eq ← R * Fcap
%       b_eq ← tauTarget - tauLimit - Σ(Fpassive)
%       a_new ← quadprog(min Σa² s.t. A_eq·a = b_eq, a∈[0.01,1])
%   until max|a_new - a| < 1e-4 or 10 iterations

% ---- parameters ----
maxIter = 10;
tolAct  = 1e-4;

nDof      = numel(dofNames);
nSoCoords = numel(soCoordNames);

% ---- load model & set kinematics ----
import org.opensim.modeling.*
soModel = Model(['model/' projName '.osim']);
nMusTotal = soModel.getMuscles().getSize();
soState = soModel.initSystem();
soCoordSet = soModel.getCoordinateSet();
for j = 1:nDof
    c = soCoordSet.get(dofNames{j});
    c.setValue(soState, initPose(j));
    c.setSpeedValue(soState, initPose(nDof + j));
end

% ---- static data: muscle handles, Fopt, coord handles ----
muscleHandles = cell(nMusTotal, 1);
Fopt = zeros(nMusTotal, 1);
for i = 1:nMusTotal
    muscleHandles{i} = soModel.getMuscles().get(i - 1);
    Fopt(i) = muscleHandles{i}.getMaxIsometricForce();
end
soCoordHandles = cell(nSoCoords, 1);
for j = 1:nSoCoords
    soCoordHandles{j} = soCoordSet.get(soCoordNames{j});
end

% ---- moment-arm matrix R (once, purely kinematic) ----
a_opt = 0.05 * ones(nMusTotal, 1);
for i = 1:nMusTotal
    Muscle.safeDownCast(soModel.getMuscles().get(i - 1)).setActivation(soState, a_opt(i));
end
soModel.equilibrateMuscles(soState);
soModel.realizeVelocity(soState);
R = zeros(nSoCoords, nMusTotal);
for j = 1:nSoCoords
    for i = 1:nMusTotal
        R(j, i) = muscleHandles{i}.computeMomentArm(soState, soCoordHandles{j});
    end
end

% ---- QP parameters (constant) ----
H = 2 * eye(nMusTotal);
f_vec = zeros(nMusTotal, 1);
lb = 0.01 * ones(nMusTotal, 1);
ub = 1.00 * ones(nMusTotal, 1);
qpOpts = optimoptions('quadprog', 'Display', 'off', 'Algorithm', 'interior-point-convex');

% ---- iterative loop ----
for iter = 1:maxIter
    % Set activations on model
    for i = 1:nMusTotal
        Muscle.safeDownCast(soModel.getMuscles().get(i - 1)).setActivation(soState, a_opt(i));
    end
    soModel.equilibrateMuscles(soState);
    soModel.realizeVelocity(soState);

    % Recompute Fcap and passive moments
    Fcap = zeros(nMusTotal, 1);
    tauPassiveMuscle = zeros(nSoCoords, 1);
    for i = 1:nMusTotal
        fl  = muscleHandles{i}.getActiveForceLengthMultiplier(soState);
        fv  = muscleHandles{i}.getForceVelocityMultiplier(soState);
        cpa = muscleHandles{i}.getCosPennationAngle(soState);
        fp  = muscleHandles{i}.getPassiveForceMultiplier(soState);
        Fcap(i) = Fopt(i) * fl * fv * cpa;
        Fpassive_i = Fopt(i) * fp * cpa;
        for j = 1:nSoCoords
            tauPassiveMuscle(j) = tauPassiveMuscle(j) + Fpassive_i * R(j, i);
        end
    end

    % Build constraints
    A_eq = R .* repmat(Fcap', nSoCoords, 1);
    b_eq = tauTarget - tauLimit - tauPassiveMuscle;

    % Solve QP
    [a_new, ~, exitflag] = quadprog(H, f_vec, [], [], A_eq, b_eq, lb, ub, [], qpOpts);
    if exitflag <= 0
        nRes = nSoCoords;
        A_aug = [A_eq, eye(nSoCoords)];
        [a_aug, ~, exitflag] = quadprog(blkdiag(H, 2000 * eye(nRes)), ...
            zeros(nMusTotal + nRes, 1), [], [], A_aug, b_eq, ...
            [lb; -1000 * ones(nRes, 1)], [ub; 1000 * ones(nRes, 1)], [], qpOpts);
        if exitflag > 0
            a_new = a_aug(1:nMusTotal);
        end
    end

    a_new = a_new(:);
    deltaA = max(abs(a_new - a_opt));
    a_opt  = a_new;
    if deltaA < tolAct, break; end
end

fprintf('[%s] Iterative SO finished: %d iters, final max|da| = %.2e,  range [%.3f, %.3f]\n', ...
    char(datetime), iter, deltaA, min(a_opt), max(a_opt));

end
