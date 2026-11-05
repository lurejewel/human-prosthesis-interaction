% =========================================================================
% Adaptive_LASSO_MocoInverse.m
%
% 作用
% ----
% 这是“MocoInverse 实测数据”分支的 LASSO main。它不修改旧的
% Adaptive_LASSO_StanceSwing.m（UN.sto 路径），而是读取接口生成的
% moco_lasso_dataset.sto/.mat，拟合 stance/swing 两组 adaptive LASSO，
% 输出一个新的控制器文件：
%
%     results/lasso_controller_result_moco.mat
%
% 输入
% ----
%   Sparse Group LASSO Validation/moco_dataset/moco_lasso_dataset.sto
%       —— 由 moco_inverse_to_lasso_dataset.m 生成；列名与 SCONE UN.sto 兼容。
%   Sparse Group LASSO Validation/moco_dataset/moco_lasso_dataset.mat
%       —— 提供 sampleWeight（排除 Moco 分段接缝/无效样本）、meta、QC。
%
% 输出
% ----
%   results/lasso_controller_result_moco.mat
%       lasso 结构，grouped v2 格式；可直接被 load_lasso_reflex_controller.m
%       读取，作为 demo 的 reflexParamMap/reflexTemplate/initPara 来源。
%   results/lasso_moco_fit_diagnostics.mat
%       fit 结构：系数、预测、lambda path、每块肌肉/相位的 RMSE、权重等。
%   Adaptive_LASSO_MocoInverse.log
%       运行日志。
%
% 与旧 main 的差异
% ----------------
%   1. 数据集已经由接口清洗/相位标注，main 不再做滤波/相位判断；
%   2. 用 sampleWeight 排除接缝/坏样本，而不是删行——lambda path 依赖
%      时间连续性，删行会让样本索引和相位边界错位；
%   3. 选定 support 后的 post-selection refit 用 Huber IRLS，避免少数
%      excitation burst 主导系数；
%   4. 输出文件名带 _moco，避免覆盖旧的 UN.sto 控制器。
%
% 运行顺序
% --------
%   1) moco_inverse_to_lasso_dataset.m
%   2) Adaptive_LASSO_MocoInverse.m
%   3) 把 demo / prepare_fresh_start 的 lassoFile 指向
%      results/lasso_controller_result_moco.mat
% =========================================================================

close all; clear; clc;

scriptFullPath = mfilename('fullpath');
if isempty(scriptFullPath)
    scriptDir = pwd;
else
    scriptDir = fileparts(scriptFullPath);
end
addpath(scriptDir);
repoRoot = fileparts(scriptDir);
resultsDir = fullfile(repoRoot, 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

logFilePath = fullfile(scriptDir, 'Adaptive_LASSO_MocoInverse.log');
if isfile(logFilePath)
    delete(logFilePath);
end
diary(logFilePath);
diary on;
logCleanup = onCleanup(@() diary('off'));

fprintf('=====================================================================\n');
fprintf(' Adaptive LASSO for MocoInverse-derived excitations\n');
fprintf(' started : %s\n', char(datetime));
fprintf('=====================================================================\n');

%% ======================= CONFIG ========================================
% 这里只放“拟合策略”参数；数据清洗参数已经在
% moco_inverse_to_lasso_dataset.m 的 cleanOpts 里，二者不要重复设置。
%
% 关键参数说明：
%   c = 0.01               censoring floor，必须与 cleanOpts.floor 一致；
%   rho = 50               censored hinge 强度；越大越强制
%                          A_C*beta+b <= c（不是 >= c）；
%   delaySteps = 1          100 Hz 网格下 10 ms 反射延迟；
%   robustIter = 2         选定 support 后的 Huber IRLS 次数；0 可做消融；
%   huberDeltaFactor=1.345 正态分布下 95% 效率的经典 Huber 常数；
%   lambdaPathConfig       nLambda/lambdaMinRatio/tolSelect 与旧 main 对齐。
datasetDir  = fullfile(scriptDir, 'moco_dataset');
datasetSto  = fullfile(datasetDir, 'moco_lasso_dataset.sto');
datasetMat  = fullfile(datasetDir, 'moco_lasso_dataset.mat');

outMat      = fullfile(resultsDir, 'lasso_controller_result_moco.mat');
diagMat     = fullfile(resultsDir, 'lasso_moco_fit_diagnostics.mat');

options = struct();
options.gamma            = 1;
options.epsilonOLS       = 1e-8;
options.nzTol            = 1e-5;
options.lambda1          = 0;
options.c                = 0.01;   % must match cleanOpts.floor
options.rho              = 50;
options.delaySteps       = 1;      % 10 ms at the 100 Hz LASSO grid
options.excitationThreshold = 0.01;
options.robustInitIter   = 2;      % Huber IRLS for initial adaptive-weight OLS
options.robustIter       = 2;      % Huber IRLS on selected support
options.huberDeltaFactor = 1.345;
options.source           = 'MocoInverse';
options.groupLabels      = {'stance', 'swing'};
options.groupPhases      = {[0 1], [2 3 4]};
options.lambdaPathConfig = struct( ...
    'nLambda', 60, ...
    'lambdaMinRatio', 1e-4, ...
    'tolSelect', 0.05, ...
    'selectionMetric', 'post_refit_train_censored_rmse');

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

% sampleWeight 的长度是“接口时间网格”的长度（在 delay 之前）。
% fit_adaptive_lasso_phase_model 会在内部按 delaySteps 自动裁剪，和
% preprocessExcitationAndFeatures 的 target 对齐，所以这里直接传整列。
S = load(datasetMat, 'ds');
ds = S.ds;
options.sampleWeight = ds.sampleWeight(:);
options.preprocessing = ds.meta;

fprintf('[input] dataset : %s\n', datasetSto);
fprintf('[input] time    : %.3f .. %.3f s (%d samples, %.1f Hz)\n', ...
    ds.time(1), ds.time(end), numel(ds.time), 1 / median(diff(ds.time)));
fprintf('[input] target  : %s\n', ds.qc.cleaning.targetSource);
fprintf('[input] valid   : %d / %d samples (%.1f %%)\n', ...
    nnz(ds.validMask), numel(ds.validMask), 100 * mean(ds.validMask));
fprintf('[input] model   : %s\n', ds.meta.modelFile);

%% ---- fit ---------------------------------------------------------------
[fit, lasso] = fit_adaptive_lasso_phase_model(datasetSto, options);

%% ---- save --------------------------------------------------------------
% 给 lasso 加来源/预处理/QC 字段，便于以后追溯“这个控制器是用哪版数据、
% 哪套清洗参数、哪个模型拟合的”。load_lasso_reflex_controller.m 会忽略
% 这些额外字段，不影响旧 loader。
lasso.source        = 'MocoInverse';
lasso.datasetFile   = datasetSto;
lasso.created       = char(datetime);
lasso.preprocessing = ds.meta;
lasso.qc            = ds.qc;
save(outMat, 'lasso');
save(diagMat, 'fit', '-v7.3');

fprintf('\n[save] controller : %s\n', outMat);
fprintf('[save] diagnostics: %s\n', diagMat);

%% ---- fit-quality plots -------------------------------------------------
if makeFitPlots
    muscleNames = fit.muscleNames;
    t = ds.time;
    stanceMask = fit.phase.stanceMask;  % after delay mapping
    Y = fit.Y;
    nY = size(Y, 1);

    % fit.Y 已经按 delaySteps 去掉第一帧，而 ds.time 还是接口全长；
    % tFit 必须与 fit.Y 逐行对应，否则相位阴影/曲线会错位。
    %   * delaySteps=0：长度相同，直接用 t；
    %   * delaySteps>0 且长度差正好是 delaySteps：t(delay+1:end)；
    %   * 其他情况（理论上不该发生）：用 linspace 兜底，只用于画图。
    if numel(stanceMask) == numel(t)
        tFit = t;
    elseif numel(t) - options.delaySteps == nY
        tFit = t((options.delaySteps + 1):end);
    else
        tFit = linspace(t(1), t(end), nY).';
    end

    for m = 1:numel(muscleNames)
        figure('Name', sprintf('LASSO fit: %s', muscleNames{m}), ...
            'Color', 'w', 'Position', [80 80 1100 420]);
        hold on; grid on;
        shade_mask(tFit, stanceMask, [0.80 0.88 0.98]);
        plot(tFit, Y(:, m), 'k-', 'LineWidth', 1.4);
        plot(tFit, fit.YPred_adaptive(:, m), '-', 'Color', [0.20 0.45 0.75], ...
            'LineWidth', 1.0);
        plot(tFit, fit.YPred_refit(:, m), '-', 'Color', [0.85 0.33 0.10], ...
            'LineWidth', 1.2);
        xlim([tFit(1) tFit(end)]);
        ylim([0 1.05]);
        xlabel('t (s)'); ylabel('excitation');
        title(sprintf('%s  |  RMSE refit = %.4g  |  stance (blue) / swing (white)', ...
            muscleNames{m}, fit.RMSE_refit_full(m)), 'Interpreter', 'none');
        legend({'target (cleaned)', 'adaptive LASSO', 'robust refit'}, ...
            'Location', 'best');
    end

    figure('Name', 'LASSO fit summary', 'Color', 'w', ...
        'Position', [100 100 1000 360]);
    bar([fit.RMSE_adaptive_full, fit.RMSE_refit_full]);
    grid on;
    set(gca, 'XTick', 1:numel(muscleNames), 'XTickLabel', muscleNames, ...
        'XTickLabelRotation', 35);
    ylabel('weighted RMSE');
    legend({'adaptive LASSO', 'robust refit'}, 'Location', 'best');
    title('Per-muscle fit error (MocoInverse target)');
end

fprintf('\nDone. Next: set lassoFile = ''%s'' in the predictive-simulation demo.\n', outMat);
fprintf('=====================================================================\n');
diary off;

%% ======================= local helpers ==================================
% 在时间轴上给 stance 区间画浅色背景。
% 只用于可视化，不参与拟合；mask 与 t 长度不一致时直接返回，避免画错。
function shade_mask(t, mask, color)
mask = logical(mask(:));
if numel(mask) ~= numel(t)
    return;
end
d = diff([false; mask; false]);
starts = find(d == 1);
finish = find(d == -1) - 1;
yl = ylim();
for r = 1:numel(starts)
    a = starts(r); b = finish(r);
    if a >= 1 && b <= numel(t)
        patch([t(a), t(b), t(b), t(a)], ...
            [yl(1), yl(1), yl(2), yl(2)], color, ...
            'EdgeColor', 'none', 'FaceAlpha', 0.45, 'HandleVisibility', 'off');
    end
end
end
