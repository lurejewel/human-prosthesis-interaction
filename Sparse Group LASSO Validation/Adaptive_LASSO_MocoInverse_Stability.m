% =========================================================================
% Adaptive_LASSO_MocoInverse_Stability.m
%
% 作用（多随机种子一致性结构辨识 / stability selection）
% ------------------------------------------------------
% 同一份实验数据、同一套辨识程序，换 K 个随机种子各跑一次 adaptive LASSO，
% 得到 K 个 26×9 的 0-1 结构（mask）。统计哪些 (feature, muscle) 条目反复
% 被选中（Pi = 入选频率），再用多数投票得到"共识结构"，最后在**全量数据**
% 上把 mask 冻结、只重新估计系数。
%
% 为什么要这样做
% --------------
% 现有链路（extractStoMuscleFeatures -> OLS -> adaptive 权重 -> lambda path
% -> Huber refit）是**完全确定性**的，换种子不会改变任何结果；真正让结构
% 抖动的是 lambda 选择那条规则（"在 min RMSE 的 5% 容差内取最稀疏"）——
% 数据动一点，选中的 lambda 就换一档，support 整块变化。因此这里用
% **按步态片段(run)分块重采样**注入随机性：每个种子 = 一次"换一批等效实验
% 样本"的合法抽样，其余超参（gamma/rho/c/tolSelect/lambda 网格/Huber）全部
% 保持不变。
%
% 输出
% ----
%   results/lasso_controller_result_moco_consensus.mat
%       lasso 结构（grouped v2），可直接被 load_lasso_reflex_controller.m
%       读取；额外带 .identification / .stability 诊断字段（loader 忽略）。
%   results/lasso_stability_diagnostics.mat
%       R：每个种子的 mask/RMSE、Pi、Jaccard 矩阵、聚类、全量重拟合结果。
%   results/lasso_stability_checkpoint.mat
%       断点续跑（options.resumeFromCheckpoint = true）。
%   results/lasso_stability_*.png
%       诊断图（Pi 热图、mask 频次、Jaccard、support size、RMSE、Pi 直方图）。
%
% 运行顺序
% --------
%   1) moco_inverse_to_lasso_dataset.m      （生成 moco_lasso_dataset.*）
%   2) 本脚本                                （得到共识 controller）
%   3) demo / prepare_fresh_start：把 lassoFile 指向
%      results/lasso_controller_result_moco_consensus.mat
%
% 说明
% ----
%   * K=20、16 核：并行大约 2 波；单个种子约 2 分钟（基线实测 126 s）。
%   * CMA-ES 只能改"有值条目"的数值，不能改 0-1 图（unpack 每帧强制
%     beta(~mask)=0，向量长度在 load 时就固定）。所以这里冻结的是结构，
%     系数只是 initial guess。
% =========================================================================

close all; clear; clc;

scriptFullPath = mfilename('fullpath');
if isempty(scriptFullPath)
    scriptDir = pwd;
else
    scriptDir = fileparts(scriptFullPath);
end
repoRoot = fileparts(scriptDir);
addpath(scriptDir);
% functions/core is only needed for the loader cross-check at the end, but the
% check is part of the contract, so put it on the path up front.
coreDir = fullfile(repoRoot, 'functions', 'core');
if exist(coreDir, 'dir')
    addpath(coreDir);
end
resultsDir = fullfile(repoRoot, 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

logFilePath = fullfile(scriptDir, 'Adaptive_LASSO_MocoInverse_Stability.log');
if isfile(logFilePath)
    delete(logFilePath);
end
diary(logFilePath);
diary on;
logCleanup = onCleanup(@() diary('off'));

fprintf('=====================================================================\n');
fprintf(' Multi-seed stability identification of the LASSO reflex structure\n');
fprintf(' started : %s\n', char(datetime));
fprintf('=====================================================================\n');

%% ======================= CONFIG ========================================
datasetDir  = fullfile(scriptDir, 'moco_dataset');
datasetSto  = fullfile(datasetDir, 'moco_lasso_dataset.sto');
datasetMat  = fullfile(datasetDir, 'moco_lasso_dataset.mat');

outMat      = fullfile(resultsDir, 'lasso_controller_result_moco_consensus.mat');
diagMat     = fullfile(resultsDir, 'lasso_stability_diagnostics.mat');
ckptMat     = fullfile(resultsDir, 'lasso_stability_checkpoint.mat');

% ---- identification settings: SHARED by every seed --------------------
% 必须与 Adaptive_LASSO_MocoInverse.m 保持一致，否则不同种子的结果不能直接
% 与单次拟合的基线比较。
fitOptions = struct();
fitOptions.gamma            = 1;
fitOptions.epsilonOLS       = 1e-8;
fitOptions.nzTol            = 1e-5;
fitOptions.lambda1          = 0;
fitOptions.c                = 0.01;   % must match cleanOpts.floor
fitOptions.rho              = 50;
fitOptions.delaySteps       = 1;      % 10 ms at the 100 Hz LASSO grid
fitOptions.excitationThreshold = 0.01;
fitOptions.robustInitIter   = 2;
fitOptions.robustIter       = 2;
fitOptions.huberDeltaFactor = 1.345;
fitOptions.source           = 'MocoInverse';
fitOptions.groupLabels      = {'stance', 'swing'};
fitOptions.groupPhases      = {[0 1], [2 3 4]};
fitOptions.lambdaPathConfig = struct( ...
    'nLambda', 60, ...
    'lambdaMinRatio', 1e-4, ...
    'tolSelect', 0.05, ...
    'selectionMetric', 'post_refit_train_censored_rmse');
% 每个 seed 只用各自的 mask，不需要逐次打印 lambda path
fitOptions.verbose          = false;

% ---- ensemble settings ------------------------------------------------
options = struct();
options.nSeeds             = 20;          % 种子数 K
options.baseSeed           = 2026;        % 第 k 个种子 = baseSeed + k
options.resampleMode       = 'run_bootstrap';   % 按步态片段有放回重采样
options.subsampleFraction  = 0.8;         % 仅 run_subsample 模式使用
options.maxBlockLen        = 0;           % 0 = 整段 run 作为一个 block

% 共识规则：默认逐条目多数投票（Pi >= 0.5），等价于"13/20 个种子给出同一
% 框架就选它"（因为 13/20 = 0.65 > 0.5，每个条目都被多数覆盖）。
% 'mode' / 'medoid' / 'stability' 只作为诊断对比。
options.consensusRule      = 'majority';
options.piThreshold        = 0.7;         % 仅 'stability' 规则使用
options.tieBreak           = 'include';   % Pi 恰好 = 0.5 时取"有值"
options.stableOnPi         = 0.9;         % 诊断: Pi >= 0.9 视为稳定激活
options.stableOffPi        = 0.1;         % 诊断: Pi <= 0.1 视为稳定不激活
options.clusterJaccard     = 0.9;         % mask 相似图连边阈值（严格聚类诊断）
options.conformJaccard     = 0.8;         % 某 seed 与共识 mask 的 Jaccard >= 此值算"支持"

% 质量门限（都不是 Pi 阈值）：
%   minStabilityJaccard : 种子之间 mask 的两两 Jaccard 均值，低于它 => 警告
%                         "0-1 结构可重复性弱"（这是真正的稳定性指标）
%   maxRMSEIncrease     : 冻结 mask 后全量重拟合 RMSE 相对**未冻结**的全量
%                         单次拟合的涨幅（种子 RMSE 是 in-bag 的，不能当参照）
%   minDominantFraction : 仅作为报告值（严格"同一框架"的占比门槛）
options.minStabilityJaccard = 0.35;
options.maxRMSEIncrease     = 0.15;
options.computeReferenceFit = true;   % 多跑一次普通全量拟合作为 RMSE 参照

options.useParallel        = true;
options.workerThreads      = 1;    % 每个 worker 只开 1 个 BLAS 线程（矩阵只有 26×26）
options.nRetry             = 3;
options.saveResults        = true;
options.resumeFromCheckpoint = false;     % 中断后改成 true 可续跑
options.verbose            = true;
options.outFile            = outMat;
options.diagFile           = diagMat;
options.checkpointFile     = ckptMat;
options.fitOptions         = fitOptions;

makeFitPlots = true;
%% =======================================================================

if ~isfile(datasetSto)
    error('DatasetNotFound', ...
        ['LASSO dataset not found:\n  %s\n' ...
         'Run moco_inverse_to_lasso_dataset.m first.'], datasetSto);
end
if ~isfile(datasetMat)
    error('DatasetMatNotFound', ...
        ['Dataset .mat not found:\n  %s\n' ...
         'Run moco_inverse_to_lasso_dataset.m first.'], datasetMat);
end

% sampleWeight 是"接口时间网格"长度（delay 之前）；loader 会按 delaySteps 裁剪。
S = load(datasetMat, 'ds');
ds = S.ds;
fitOptions.sampleWeight  = ds.sampleWeight(:);
fitOptions.preprocessing = ds.meta;
options.fitOptions       = fitOptions;

fprintf('[input] dataset   : %s\n', datasetSto);
fprintf('[input] time      : %.3f .. %.3f s (%d samples, %.1f Hz)\n', ...
    ds.time(1), ds.time(end), numel(ds.time), 1 / median(diff(ds.time)));
fprintf('[input] valid     : %d / %d samples (%.1f %%)\n', ...
    nnz(ds.validMask), numel(ds.validMask), 100 * mean(ds.validMask));
fprintf('[config] seeds    : %d (baseSeed=%d)\n', options.nSeeds, options.baseSeed);
fprintf('[config] consensus: %s (Pi >= 0.5 majority; tie=%s)\n', ...
    options.consensusRule, options.tieBreak);

%% ---- run ---------------------------------------------------------------
t0 = tic;
[R, D] = lasso_stability_identification(datasetSto, options);
elapsed = toc(t0);

%% ---- report ------------------------------------------------------------
lassoFinal = R.lasso;
C = R.consensus;
muscleNames = lassoFinal.muscleNames;

fprintf('\n===== FROZEN 0-1 STRUCTURE (consensus) =====\n');
fprintf('  %-4s %-10s %6s %6s %8s %8s\n', 'grp', 'label', 'nnz', 'of', 'median', 'pi>=.9');
for g = 1:numel(C.mask)
    fprintf('  %-4d %-10s %6d %6d %8g %8d\n', g, ...
        lassoFinal.identification.groupLabels{g}, nnz(C.mask{g}), C.nElem, ...
        median(C.supportSizes{g}), C.nStableOn(g));
end
fprintf('\n  %-12s %8s %8s %8s\n', 'muscle', 'stance', 'swing', 'total');
totalNnz = 0;
for m = 1:numel(muscleNames)
    n1 = nnz(C.mask{1}(:, m));
    n2 = nnz(C.mask{2}(:, m));
    totalNnz = totalNnz + n1 + n2;
    fprintf('  %-12s %8d %8d %8d\n', muscleNames{m}, n1, n2, n1 + n2);
end
fprintf('  %-12s %8d %8d %8d\n', '(total)', ...
    nnz(C.mask{1}), nnz(C.mask{2}), totalNnz);
fprintf('\n  CMA-ES vector length = %d  (nnz_stance + 9 + nnz_swing + 9)\n', ...
    R.totalLen);
fprintf('  elapsed = %.1f s (%.1f min)\n', elapsed, elapsed / 60);

%% ---- verify the emitted controller loads -------------------------------
try
    [initPara, reflexParamMap, reflexTemplate] = ...
        load_lasso_reflex_controller(outMat);
    fprintf('\n[check] load_lasso_reflex_controller OK: totalLen=%d, %d groups, %d phases\n', ...
        reflexParamMap.totalLen, numel(reflexParamMap.groups), reflexTemplate.nPhases);
    assert(reflexParamMap.totalLen == R.totalLen, ...
        'Total parameter count mismatch: loader %d vs expected %d.', ...
        reflexParamMap.totalLen, R.totalLen);
    fprintf('[check] totalLen matches the consensus structure.\n');
catch ME
    warning('Stability:LoaderCheckFailed', ...
        'Emitted controller did not load as expected: %s', ME.message);
end

%% ---- plots -------------------------------------------------------------
if makeFitPlots
    plot_lasso_stability_report(R, resultsDir);
    fprintf('[save] figures: %s\\lasso_stability_*.png\n', resultsDir);
end

fprintf('\nDone.\n');
fprintf('Next: set lassoFile = ''%s'' in prepare_fresh_start.m / the demo.\n', outMat);
fprintf('      The 0-1 map is frozen; CMA-ES only tunes the values inside it.\n');
fprintf('=====================================================================\n');
diary off;
