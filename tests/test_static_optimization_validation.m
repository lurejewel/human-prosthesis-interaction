% -------------------------------------------------------------------------
% Name: test_static_optimization_validation.m
% Author(s): Jin Wei, Peking U. wjin24@stu.pku.edu.cn
% Description: Validate static-optimisation-derived initial muscle
%   activations against UN.sto reference data.  For the first frame of
%   UN.sto, the script:
%     1. Reads joint angles, joint moments, knee limit torques, and
%        muscle activations.
%     2. Sets the human0714.osim model to that static pose.
%     3. Builds and solves a convex QP:
%           min  Σ a_m²
%           s.t. Σ (a_m · F_cap_m · r_{m,j}) = τ_target(j) − τ_limit(j)
%        where F_cap = F_opt · f_l · cos(α)  (v = 0 ⇒ f_v = 1).
%     4. Compares QP activations with UN.sto activations via scatter and
%        bar plots.
%
% NOTE: The knee joint in human0714.osim has a CoordinateLimitForce.
%   Its contribution (knee_r.torque / knee_l.torque) is read from UN.sto
%   and subtracted from the net knee moment so that the QP only needs to
%   match the *muscle-generated* portion.
% -------------------------------------------------------------------------

function test_static_optimization_validation()
import org.opensim.modeling.*

%% ====== 0. Paths ======
[scriptDir, ~, ~] = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..');
cd(projectRoot);
addpath(genpath('assets'), genpath('model'), genpath('functions'));

stoPath = fullfile('Sparse Group LASSO Validation', 'UN.sto');
assert(isfile(stoPath), 'UN.sto not found: %s', stoPath);

fprintf('========== Static Optimisation Validation ==========\n');

%% ====== 1. Read UN.sto first frame ======
fprintf('[1/5] Reading UN.sto first frame ...\n');

% ---- manual tab-delimited read (no OpenSim Storage dependency) ----
fid = fopen(stoPath, 'r');
for k = 1:6
    fgetl(fid);  % skip metadata lines (line 1-6, line 6 = "endheader")
end
headerLine = fgetl(fid);           % line 7: column labels
dataLine   = str2double(strsplit(fgetl(fid), '\t'));  % line 8: first data row (t=0)
fclose(fid);

colNames = strsplit(headerLine, '\t');
nCols = numel(colNames);
colMap = containers.Map('KeyType', 'char', 'ValueType', 'int32');
for i = 1:nCols
    colMap(colNames{i}) = i;
end

% ---- load model early to extract muscle names & count ----
projName = 'human0714';
model = Model(fullfile('model', [projName '.osim']));
nMusTotal = model.getMuscles().getSize();
nMusPerLeg = nMusTotal / 2;

% Extract base muscle names from right-leg muscles (strip '_r' suffix)
muscleBaseNames = cell(nMusPerLeg, 1);
for i = 1:nMusPerLeg
    fullName = char(model.getMuscles().get(i - 1).getName());
    muscleBaseNames{i} = fullName(1:end-2);  % strip '_r'
end

sides = {'_r', '_l'};

% ---- 1a. Coordinate definitions (FIXED order, identical in UN.sto and .osim) ----
% All vectors & matrices below use this exact order.
% 1: hip R   2: hip L   3: knee R   4: knee L   5: ankle R   6: ankle L
coordNames = {'hip_flexion_r', 'hip_flexion_l', ...
              'knee_flexion_r', 'knee_flexion_l', ...
              'ankle_dorsiflexion_r', 'ankle_dorsiflexion_l'};
nCoords = numel(coordNames);

% ---- 1b. Joint angles and velocities (rad, rad/s, used to set model pose) ----
qTarget = zeros(nCoords, 1);
uTarget = zeros(nCoords, 1);
for j = 1:nCoords
    qTarget(j) = readCol(coordNames{j},         colMap, dataLine);
    uTarget(j) = readCol([coordNames{j}, '_u'],  colMap, dataLine);  % UN.sto velocity suffix
end
fprintf('  Joint kinematics loaded for %d coordinates.\n', nCoords);

% ---- 1c. Joint moments (target for QP) — SAME order as coordNames ----
momentCols = strcat(coordNames, '.moment');  % {pelvis_tilt.moment, hip_flexion_r.moment, ...}
tauTarget = zeros(nCoords, 1);
for j = 1:nCoords
    tauTarget(j) = readCol(momentCols{j}, colMap, dataLine);
end

% ---- 1d. Knee limit-force torques (from UN.sto) — SAME order ----
tauLimit = zeros(nCoords, 1);
kneeIdxR = find(strcmp(coordNames, 'knee_flexion_r'));
kneeIdxL = find(strcmp(coordNames, 'knee_flexion_l'));
tauLimit(kneeIdxR) = readCol('knee_r.torque', colMap, dataLine);
tauLimit(kneeIdxL) = readCol('knee_l.torque', colMap, dataLine);

% ---- 1e. Muscle activations (reference for comparison) ----
% Names derived from model (muscleBaseNames extracted above)
actUN = zeros(nMusTotal, 1);
muscleLabels = cell(nMusTotal, 1);  % for plotting

for s = 1:2  % 1=right, 2=left
    base = (s - 1) * nMusPerLeg;
    for m = 1:nMusPerLeg
        idx = base + m;
        colName = [muscleBaseNames{m}, sides{s}, '.activation'];
        actUN(idx) = readCol(colName, colMap, dataLine);
        muscleLabels{idx} = [muscleBaseNames{m}, sides{s}];
    end
end

fprintf('  Reference activations loaded for %d muscles.\n', nMusTotal);

%% ====== 2. Iterative Static Optimization ======
fprintf('[2/3] Iterative Static Optimization ...\n');

% ---- parameters ----
maxIter = 10;
tolAct  = 1e-4;

% ---- initialize state and set kinematics (once) ----
state = model.initSystem();
coordSet = model.getCoordinateSet();
for j = 1:nCoords
    cname = coordNames{j};
    coord = coordSet.get(cname);
    coord.setValue(state, qTarget(j));
    coord.setSpeedValue(state, uTarget(j));
end

% ---- build static data: muscle handles, Fopt, coordHandles ----
muscleHandles = cell(nMusTotal, 1);
Fopt  = zeros(nMusTotal, 1);
for i = 1:nMusTotal
    muscleHandles{i} = model.getMuscles().get(i - 1);
    Fopt(i) = muscleHandles{i}.getMaxIsometricForce();
end

coordHandles = cell(nCoords, 1);
for j = 1:nCoords
    coordHandles{j} = coordSet.get(coordNames{j});
end

% ---- compute moment-arm matrix R (once — purely kinematic) ----
a_current = 0.05 * ones(nMusTotal, 1);
for i = 1:nMusTotal
    Muscle.safeDownCast(model.getMuscles().get(i - 1)).setActivation(state, a_current(i));
end
model.equilibrateMuscles(state);
model.realizeVelocity(state);
R = zeros(nCoords, nMusTotal);
for j = 1:nCoords
    for i = 1:nMusTotal
        R(j, i) = muscleHandles{i}.computeMomentArm(state, coordHandles{j});
    end
end

% ---- QP parameters (constant across iterations) ----
H = 2 * eye(nMusTotal);
f_vec = zeros(nMusTotal, 1);
lb = 0.01 * ones(nMusTotal, 1);
ub = 1.00 * ones(nMusTotal, 1);
opts = optimoptions('quadprog', 'Display', 'off', 'Algorithm', 'interior-point-convex');

% ---- iterative loop ----
convHistory = nan(maxIter, 1);
reserveUsed = false;
qpIterations = nan(maxIter, 1);

for iter = 1:maxIter
    % (a) Set activations on model
    for i = 1:nMusTotal
        Muscle.safeDownCast(model.getMuscles().get(i - 1)).setActivation(state, a_current(i));
    end
    model.equilibrateMuscles(state);
    model.realizeVelocity(state);

    % (b) Recompute Fcap and passive moments
    Fcap = zeros(nMusTotal, 1);
    tauPassiveMuscle = zeros(nCoords, 1);
    for i = 1:nMusTotal
        fl  = muscleHandles{i}.getActiveForceLengthMultiplier(state);
        fv  = muscleHandles{i}.getForceVelocityMultiplier(state);
        cpa = muscleHandles{i}.getCosPennationAngle(state);
        fp  = muscleHandles{i}.getPassiveForceMultiplier(state);
        Fcap(i) = Fopt(i) * fl * fv * cpa;
        Fpassive_i = Fopt(i) * fp * cpa;
        for j = 1:nCoords
            tauPassiveMuscle(j) = tauPassiveMuscle(j) + Fpassive_i * R(j, i);
        end
    end

    % (c) Build constraints
    A_eq = R .* repmat(Fcap', nCoords, 1);
    b_eq = tauTarget - tauLimit - tauPassiveMuscle;

    % (d) Solve QP
    [a_new, ~, exitflag, qpOut] = quadprog(H, f_vec, [], [], A_eq, b_eq, lb, ub, [], opts);

    if exitflag <= 0
        nRes = nCoords;
        H_aug = blkdiag(H, 2000 * eye(nRes));
        A_aug = [A_eq, eye(nCoords)];
        lb_aug = [lb; -1000 * ones(nRes, 1)];
        ub_aug = [ub;  1000 * ones(nRes, 1)];
        [a_aug, ~, exitflag, qpOut] = quadprog(H_aug, f_vec, [], [], A_aug, b_eq, ...
            lb_aug, ub_aug, [], opts);
        if exitflag > 0
            a_new = a_aug(1:nMusTotal);
            reserveUsed = true;
        else
            warning('QP failed at iter %d (exitflag=%d). Using previous activations.', iter, exitflag);
            a_new = a_current;
        end
    end

    a_new = a_new(:);
    deltaA = max(abs(a_new - a_current));
    convHistory(iter) = deltaA;
    qpIterations(iter) = qpOut.iterations;

    fprintf('  Iter %d/%d: max|da| = %.2e, QP exitflag = %d', ...
        iter, maxIter, deltaA, exitflag);
    if exitflag > 0
        fprintf(' (%d QP iterations)', qpOut.iterations);
    end
    fprintf('\n');

    a_current = a_new;
    if deltaA < tolAct
        fprintf('  Converged at iteration %d.\n', iter);
        break
    end
end

a_opt = a_current;
convHistory = convHistory(1:iter);
qpIterations = qpIterations(1:iter);

% ---- fiber-length summary (final converged state) ----
lCEN_final = zeros(nMusTotal, 1);
for i = 1:nMusTotal
    lCEN_final(i) = muscleHandles{i}.getNormalizedFiberLength(state);
end
fprintf('  Final normalized fiber lengths: min=%.4f, max=%.4f, mean=%.4f\n', ...
    min(lCEN_final), max(lCEN_final), mean(lCEN_final));

% ---- iteration summary ----
fprintf('  Iterative SO finished: %d iterations, final max|da| = %.2e', iter, convHistory(end));
if reserveUsed, fprintf(', reserves needed'); end
fprintf('.\n');
fprintf('  Activation range: [%.4f, %.4f].\n', min(a_opt), max(a_opt));

%% ====== 3. Compare & visualise ======
fprintf('[3/3] Plotting comparisons ...\n');

% ---- 3a. Joint moment reconstruction check ----
tauRecon = A_eq * a_opt;
tauError = tauRecon - b_eq;
fprintf('\nJoint moment reconstruction error (should be ≈ 0):\n');
coordLabelsShort = {'hip R', 'hip L', 'knee R', 'knee L', ...
                    'ankle R', 'ankle L'};
for j = 1:nCoords
    fprintf('  %-12s  target = %8.2f  recon = %8.2f  err = %+.2e Nm\n', ...
        coordLabelsShort{j}, tauTarget(j), tauRecon(j) + tauLimit(j) + tauPassiveMuscle(j), tauError(j));
end

% ---- 3b. Scatter plot: QP vs UN.sto activations ----
actUN = actUN(:);   % ensure column vector
figure('Color', 'w', 'Position', [100, 100, 1100, 500]);

subplot(1, 2, 1);
scatter(actUN, a_opt, 50, 'filled', 'MarkerEdgeColor', 'k');
hold on;
mx = max(max(actUN), max(a_opt));
plot([0, mx], [0, mx], 'k--', 'LineWidth', 1);   % ideal line
xlabel('UN.sto activation');
ylabel('QP computed activation');
title('Static Optimisation: QP vs UN.sto');
grid on; axis equal; xlim([0, mx*1.05]); ylim([0, mx*1.05]);

% label each point
for i = 1:nMusTotal
    text(actUN(i) + 0.01, a_opt(i) + 0.01, muscleLabels{i}, ...
        'FontSize', 7, 'Interpreter', 'none');
end

% ---- stats ----
mask = ~isnan(actUN) & ~isnan(a_opt);
rmse = sqrt(mean((actUN(mask) - a_opt(mask)).^2));
rho  = corr(actUN(mask), a_opt(mask));
text(mx*0.02, mx*0.92, sprintf('RMSE = %.4f\nr = %.4f', rmse, rho), ...
    'FontSize', 10, 'BackgroundColor', [1 1 0.85], 'VerticalAlignment', 'top');

% ---- 3c. Grouped bar chart (right / left legs) ----
subplot(1, 2, 2);
nBars = nMusTotal;
barWidth = 0.35;
xPos = (1:nBars)';

b1 = bar(xPos - barWidth/2, actUN, barWidth, 'FaceColor', [0.4 0.6 0.9]);
hold on;
b2 = bar(xPos + barWidth/2, a_opt,   barWidth, 'FaceColor', [0.9 0.4 0.4]);

set(gca, 'XTick', 1:nBars, 'XTickLabel', muscleLabels, ...
    'XTickLabelRotation', 45, 'FontSize', 8);
ylabel('Activation');
legend([b1, b2], {'UN.sto', 'QP (Static Opt)'}, 'Location', 'northwest');
title(sprintf('Muscle Activation Comparison  (RMSE = %.4f)', rmse));
grid on; box on;

% add a vertical separator between right & left legs
xline(nMusPerLeg + 0.5, 'k--', 'LineWidth', 1.2);

sgtitle(sprintf('Static Optimisation Validation  |  Frame 0  |  %d muscles', nMusTotal), ...
    'FontSize', 11, 'FontWeight', 'bold');

% ---- 3d. Joint moment comparison (3-way) ----
tauReconUN = A_eq * actUN + tauLimit + tauPassiveMuscle;   % UN activations → moments
tauReconQP = A_eq * a_opt + tauLimit + tauPassiveMuscle;   % QP activations → moments

figure('Color', 'w', 'Position', [150, 150, 1000, 500]);
nJoints = nCoords;
barWidth = 0.22;
xPos = (1:nJoints)';

b1 = bar(xPos - barWidth, tauTarget,   barWidth, 'FaceColor', [0.4 0.6 0.9]);  hold on;
b2 = bar(xPos,           tauReconQP,   barWidth, 'FaceColor', [0.9 0.4 0.4]);
b3 = bar(xPos + barWidth, tauReconUN,   barWidth, 'FaceColor', [0.4 0.9 0.4]);

set(gca, 'XTick', 1:nJoints, 'XTickLabel', coordLabelsShort, 'FontSize', 10);
ylabel('Joint moment (Nm)');
legend([b1, b2, b3], {'UN.sto net moment', 'QP recon (A·a_{QP})', 'UN act. recon (A·a_{UN})'}, ...
    'Location', 'best');
title('Joint Moment Comparison: Target vs Reconstructed');
grid on; box on;

% ---- annotation: reconstruction error ----
fprintf('\nJoint moment reconstruction summary:\n');
fprintf('  %-12s  %10s  %10s  %10s\n', 'Coordinate', 'UN target', 'QP recon', 'UN act recon');
for j = 1:nJoints
    fprintf('  %-12s  %10.2f  %10.2f  %10.2f\n', ...
        coordLabelsShort{j}, tauTarget(j), tauReconQP(j), tauReconUN(j));
end

% ---- 3e. Convergence history plot ----
figure('Color', 'w', 'Position', [200, 200, 500, 400]);
semilogy(1:iter, convHistory, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Iteration'); ylabel('max|da| (inf-norm of activation change)');
title(sprintf('Iterative SO Convergence  (tol = %.0e)', tolAct));
grid on; box on;
yline(tolAct, 'r--', 'LineWidth', 1.2);
legend({'max|da|', 'tolerance'}, 'Location', 'best');

fprintf('\n========== Summary ==========\n');
fprintf('RMSE (QP vs UN.sto)     : %.4f\n', rmse);
fprintf('Correlation (r)          : %.4f\n', rho);
fprintf('Iterative SO converged in: %d iterations\n', iter);
fprintf('Final max|da|            : %.2e\n', convHistory(end));
if reserveUsed
    fprintf('Reserve actuators        : needed\n');
end
fprintf('==============================\n\n');

end

% ========================================================================
function val = readCol(name, colMap, dataLine)
% Read a single value from the first data row of UN.sto by column name.
    idx = colMap(name);
    val = dataLine(idx);
end
