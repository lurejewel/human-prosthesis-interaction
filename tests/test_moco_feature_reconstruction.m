% =========================================================================
% test_moco_feature_reconstruction.m
%
% 目的
% ----
% 验证 build_lasso_features_from_kinematics.m 的“由 q/qdot/activation 重建
% lCEN”是否可信。验证数据用 results/sim_result_*.sto：这是预测仿真自己的
% 前向结果，里面同时有 q、qdot、activation 和真实的 fiber_length state。
%
% 思路
% ----
%   1. 从 forward simulation 结果里取 q/qdot/activation + 真实 fiber_length；
%   2. 只把 q/qdot/activation 传给 build_lasso_features；
%   3. 重建 lCEN = fiber_length / lopt；
%   4. 与文件里的真实 fiber_length/lopt 比较。
%
% 为什么跳过前 0.5 s
% ------------------
% 接口没有第一帧的 fiber_length，只能用 lopt + equilibrateMuscles 求初值。
% 前 0.5 s 包含这个初值 transient；0.5 s 后积分已经“忘掉”初值。接口在
% LASSO 里也用 edgeMargin 排除首尾，所以这里跳过 0.5 s 与真实使用一致。
%
% 判据
% ----
% 子步积分（<=2 ms）对 18 块肌肉的 lCEN RMSE 应 ≲1e-3；这里放宽到 0.01，
% 以容忍不同 OpenSim 版本/equilibrate 的数值差异。若 RMSE 明显大于 0.01，
% 通常说明：
%   * 子步太大（substepMax 调大后）；
%   * 坐标/速度单位不是 rad/m；
%   * 模型不是 human0918（肌肉参数不同）。
%
% 运行
% ----
%   从仓库根目录：
%       test_moco_feature_reconstruction
% =========================================================================

clear; clc;

% -------------------------------------------------------------------------
% 0. 路径与测试数据
% -------------------------------------------------------------------------
% 测试脚本假定从仓库根目录运行；用 mfilename 推导 repoRoot，避免依赖 pwd。
here = fileparts(mfilename('fullpath'));
repoRoot = fileparts(here);
addpath(fullfile(repoRoot, 'Sparse Group LASSO Validation'));
addpath(fullfile(repoRoot, 'functions'));

simFile = fullfile(repoRoot, 'results', 'sim_result_2026-07-27_15-39-45.sto');
if ~isfile(simFile)
    error('TestDataMissing', 'Forward-simulation reference not found: %s', simFile);
end
T = read_sto_mot_table(simFile);
lab = T.labels;
% 列名 -> 列号；后面直接用 idx('pelvis_tilt/value') 这类写法。
idx = containers.Map(lab, num2cell(1:numel(lab)));

% -------------------------------------------------------------------------
% 1. 取出 9 个坐标的 value/speed
% -------------------------------------------------------------------------
% human0918 是 2D 模型，只有这 9 个 DOF；forward sim 结果也按这个顺序保存。
qNames = {'pelvis_tilt','pelvis_tx','pelvis_ty', ...
          'hip_flexion_r','knee_extension_r','ankle_dorsiflexion_r', ...
          'hip_flexion_l','knee_extension_l','ankle_dorsiflexion_l'};
q  = zeros(T.nRows, numel(qNames));
qd = zeros(T.nRows, numel(qNames));
for c = 1:numel(qNames)
    q(:, c)  = T.data(:, idx([qNames{c} '/value']));
    qd(:, c) = T.data(:, idx([qNames{c} '/speed']));
end

% -------------------------------------------------------------------------
% 2. 取出 18 块肌肉的 activation 与真实 fiber_length
% -------------------------------------------------------------------------
% muscleNames 顺序为 canonical right 9 + left 9，与 build_lasso_features
% 的列顺序假设一致。
muscleNames = {'hamstrings_r','bifemsh_r','glut_max_r','iliopsoas_r', ...
               'rect_fem_r','vasti_r','gastroc_r','soleus_r','tib_ant_r', ...
               'hamstrings_l','bifemsh_l','glut_max_l','iliopsoas_l', ...
               'rect_fem_l','vasti_l','gastroc_l','soleus_l','tib_ant_l'};
act = zeros(T.nRows, numel(muscleNames));
trueFL = zeros(T.nRows, numel(muscleNames));
for m = 1:numel(muscleNames)
    act(:, m)    = T.data(:, idx(['/forceset/' muscleNames{m} '/activation']));
    trueFL(:, m) = T.data(:, idx(['/forceset/' muscleNames{m} '/fiber_length']));
end

% -------------------------------------------------------------------------
% 3. 只给 q/qdot/activation，重建特征
% -------------------------------------------------------------------------
% substepMax=0.002 是接口默认值；verbose 打开便于看到子步数和特征范围。
F = build_lasso_features_from_kinematics( ...
    fullfile(repoRoot, 'model', 'human0918.osim'), T.time, q, qd, act, ...
    muscleNames, qNames, struct('substepMax', 0.002, 'verbose', true));

% -------------------------------------------------------------------------
% 4. 读取模型 lopt，把真实 fiber_length 转成 lCEN
% -------------------------------------------------------------------------
% 注意：这里必须用同一个模型（human0918）的 lopt，否则比较没有意义。
mdl = org.opensim.modeling.Model(fullfile(repoRoot, 'model', 'human0918.osim'));
lopt = zeros(1, numel(muscleNames));
for m = 1:numel(muscleNames)
    lopt(m) = mdl.getMuscles().get(muscleNames{m}).getOptimalFiberLength();
end

% -------------------------------------------------------------------------
% 5. 跳过初始 transient，逐肌肉比较
% -------------------------------------------------------------------------
skip = T.time < (T.time(1) + 0.5);
rmse = zeros(1, numel(muscleNames));
maxErr = zeros(1, numel(muscleNames));
for m = 1:numel(muscleNames)
    ref = trueFL(~skip, m) / lopt(m);
    got = F.lCEN(~skip, m);
    rmse(m) = sqrt(mean((got - ref).^2));
    maxErr(m) = max(abs(got - ref));
    fprintf('%-12s  RMSE=%.2e  max|err|=%.2e\n', muscleNames{m}, rmse(m), maxErr(m));
end

% -------------------------------------------------------------------------
% 6. 判据
% -------------------------------------------------------------------------
fprintf('\nmax RMSE over muscles = %.3e\n', max(rmse));
if max(rmse) > 0.01
    error('FeatureReconstructionFailed', ...
        'lCEN reconstruction RMSE %.3e exceeds 0.01.', max(rmse));
end
fprintf('Feature reconstruction test PASSED.\n');
