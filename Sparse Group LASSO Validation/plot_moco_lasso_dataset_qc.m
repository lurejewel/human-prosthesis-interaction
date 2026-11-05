function plot_moco_lasso_dataset_qc(ds, outDir, outStem)
%PLOT_MOCO_LASSO_DATASET_QC 保存 Moco->LASSO 数据集的 QC 总览图。
%
% 生成三张图：
%   1) *_qc_drives.png   ：右腿 9 块肌肉的 activation、raw excitation（若有）
%                          与清洗后 target 的对比，并用浅蓝色标出 stance；
%   2) *_qc_phase.png    ：用于相位判断的竖直 GRF、stance mask、0..4 标签；
%   3) *_qc_features.png ：lCEN / fATN 的肌肉-时间热图。
%
% 这些图不参与拟合，只用于快速判断：
%   * target 是否还有明显尖刺/饱和；
%   * stance 是否被切得太碎；
%   * 特征是否出现异常大/负值或某块肌肉整段不变。
%
% 输入：
%   ds      接口生成的 dataset struct（含 time/activation/excitation/lCEN/
%           fATN/phase/qc/muscleNames/rawExcitation 等）；
%   outDir  图片输出目录；
%   outStem 文件名前缀。
%
% 注意：所有绘图都做了 try/catch 的外层保护（调用方），单张图失败不会
% 影响 .sto/.mat 已经写出的结果。

if nargin < 2 || isempty(outDir)
    outDir = pwd;
end
if nargin < 3 || isempty(outStem)
    outStem = 'moco_lasso_dataset';
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

t = ds.time;
nT = numel(t);
rightIdx = 1:min(9, numel(ds.muscleNames));
muscleNames = ds.muscleNames(rightIdx);

% ---- figure 1: raw vs cleaned drives ------------------------------------
% 只画右腿 9 块（LASSO 实际拟合的对象）；左腿列存在但不画，避免图太满。
% 三条线：
%   activation（蓝）     —— 清洗后的激活，用来判断 target 是否跟得上
%   raw x（灰点线）      —— 仅当接口有 raw excitation 时
%   target x（橙）       —— LASSO 真正拟合的目标
% 浅蓝背景 = stance；空白 = swing。
f1 = figure('Name', 'Moco -> LASSO drives', 'Color', 'w', ...
    'Position', [80 80 1400 850]);
tl = tiledlayout(f1, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:numel(rightIdx)
    j = rightIdx(k);
    nexttile; hold on; grid on;
    shade_stance(t, ds.phase.stanceR);
    if ds.qc.cleaning.hasRawActivation
        plot(t, ds.activation(:, j), '-', 'Color', [0.20 0.45 0.75], ...
            'LineWidth', 1.1);
    end
    if isfield(ds, 'rawExcitation') && ~isempty(ds.rawExcitation)
        plot(t, ds.rawExcitation(:, j), ':', 'Color', [0.55 0.55 0.55], ...
            'LineWidth', 0.9);
    end
    plot(t, ds.excitation(:, j), '-', 'Color', [0.85 0.33 0.10], ...
        'LineWidth', 1.3);
    ylim([0 1.05]);
    xlim([t(1) t(end)]);
    title(muscleNames{j}, 'Interpreter', 'none');
    if k == 1
        leg = {'activation'};
        if isfield(ds, 'rawExcitation') && ~isempty(ds.rawExcitation)
            leg{end+1} = 'raw x';
        end
        leg{end+1} = 'LASSO target x';
        legend(leg, 'Location', 'southoutside', 'FontSize', 7);
    end
    if k > numel(rightIdx) - 3, xlabel('t (s)'); end
    if mod(k - 1, 3) == 0, ylabel('(-)'); end
end
sgtitle(sprintf('Right-leg drives: target = %s', ds.qc.cleaning.targetSource), ...
    'Interpreter', 'none');
save_fig(f1, fullfile(outDir, [outStem '_qc_drives.png']));

% ---- figure 2: GRF and state labels -------------------------------------
% 第一行：实际使用的左右竖直 GRF，检查 0.02*BW 阈值是否合理；
% 第二行：LASSO stance mask（state 0/1；左腿整体上移 1.2）；
% 第三行：严格按 cal_gait_phase.m 计算的 0..4（左腿整体上移 5.2）。
f2 = figure('Name', 'GRF and phase labels', 'Color', 'w', ...
    'Position', [120 120 1100 650]);
tiledlayout(f2, 3, 1, 'TileSpacing', 'compact');
nexttile; hold on; grid on;
plot(t, ds.phase.forceR, '-', 'LineWidth', 1.1);
plot(t, ds.phase.forceL, '--', 'LineWidth', 1.0);
ylabel('F_y (N)'); legend({'right','left'}, 'Location', 'best');
title('Vertical GRF used for stance detection');
nexttile; hold on; grid on;
plot(t, ds.phase.stanceR, '-', 'LineWidth', 1.4);
plot(t, ds.phase.stanceL + 1.2, '--', 'LineWidth', 1.0);
ylim([-0.1 2.3]); yticks([0 1 1.2 2.2]);
yticklabels({'swing','stance','swing','stance'});
title('LASSO stance masks = state 0/1 (right solid, left dashed)');
nexttile; hold on; grid on;
plot(t, ds.phase.stateR, '-', 'LineWidth', 1.2);
plot(t, ds.phase.stateL + 5.2, '--', 'LineWidth', 1.0);
ylim([-0.2 9.6]); yticks(0:4);
title('cal_gait_phase 0..4 labels (right solid, left dashed offset)', ...
    'Interpreter', 'none');
xlabel('t (s)');
save_fig(f2, fullfile(outDir, [outStem '_qc_phase.png']));

% ---- figure 3: features ------------------------------------------------
% 用 imagesc 看 9 块肌肉 × 时间的 lCEN/fATN。正常情况：
%   * lCEN 大致在 0.5–1.4，随步态周期上下波动；
%   * fATN 在 0–0.8，stance 期明显高于 swing；
%   * 不应出现某一块肌肉整段贴着 0 或 2.5（后者是积分发散被 clip 的迹象）。
f3 = figure('Name', 'LASSO features', 'Color', 'w', ...
    'Position', [100 100 1200 700]);
tiledlayout(f3, 2, 1, 'TileSpacing', 'compact');
nexttile;
imagesc(t, 1:numel(rightIdx), ds.lCEN(:, rightIdx).');
axis xy; colorbar; yticks(1:numel(rightIdx));
yticklabels(regexprep(muscleNames, '_r$', ''));
title('Normalized fiber length (lCEN)'); xlabel('t (s)');
nexttile;
imagesc(t, 1:numel(rightIdx), ds.fATN(:, rightIdx).');
axis xy; colorbar; yticks(1:numel(rightIdx));
yticklabels(regexprep(muscleNames, '_r$', ''));
title('Normalized MTU force (fATN)'); xlabel('t (s)');
save_fig(f3, fullfile(outDir, [outStem '_qc_features.png']));
end

% ========================================================================
% 局部函数：在坐标轴上给 stance 区间画浅色背景
% ========================================================================
% 用 patch 而不是 area，避免影响坐标轴范围；HandleVisibility='off'
% 保证 legend 不会把背景当作数据系列。
function shade_stance(t, mask)
mask = logical(mask(:));
d = diff([false; mask; false]);
starts = find(d == 1);
ends   = find(d == -1) - 1;
yl = ylim();
for r = 1:numel(starts)
    if starts(r) <= numel(t) && ends(r) <= numel(t)
        patch([t(starts(r)), t(ends(r)), t(ends(r)), t(starts(r))], ...
            [yl(1), yl(1), yl(2), yl(2)], [0.85 0.90 0.98], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.45, 'HandleVisibility', 'off');
    end
end
end

% ========================================================================
% 局部函数：保存图片（优先 exportgraphics，旧 MATLAB 用 print）
% ========================================================================
function save_fig(fig, file)
try
    exportgraphics(fig, file, 'Resolution', 150);
catch
    print(fig, file, '-dpng', '-r150');
end
end
