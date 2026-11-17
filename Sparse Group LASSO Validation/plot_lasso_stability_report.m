function plot_lasso_stability_report(R, figDir)
%PLOT_LASSO_STABILITY_REPORT  Diagnostic figures for the consensus structure.
%
%   plot_lasso_stability_report(R, figDir)
%
%   R is the struct returned by lasso_stability_from_design /
%   lasso_stability_identification (or loaded from
%   results/lasso_stability_diagnostics.mat).  Figures are written as PNGs to
%   figDir (default: current folder).
%
%   Figures
%   -------
%   1. Pi heatmaps per phase group (feature x muscle selection frequency);
%   2. exact-mask frequency spectrum ("how many seeds agree on the same mask");
%   3. Jaccard similarity matrix per group (with the detected clusters);
%   4. support-size distribution of the replicates vs the consensus mask;
%   5. per-muscle RMSE: replicate spread vs the frozen-mask full-data refit;
%   6. Pi histogram (stable-on / ambiguous / stable-off entry counts).
%
%   See also lasso_stability_from_design, consensus_reflex_structure.

if nargin < 2 || isempty(figDir)
    figDir = pwd;
end
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

if ~isfield(R, 'consensus') || ~isfield(R, 'lasso')
    error('plot_lasso_stability_report:BadInput', ...
        'R must contain .consensus and .lasso (see lasso_stability_from_design).');
end
C = R.consensus;
lasso = R.lasso;
nGroups = numel(C.mask);
groupLabels = get_group_labels(R, nGroups);
muscleNames = cellstr(lasso.muscleNames);
featureNames = cellstr(lasso.featureNames);
if numel(featureNames) ~= size(C.Pi{1}, 1)
    featureNames = arrayfun(@(i) sprintf('f%d', i), 1:size(C.Pi{1}, 1), ...
        'UniformOutput', false).';
end
featureLabels = short_feature_names(featureNames);
muscleLabels = short_muscle_names(muscleNames);

%% ---- Figure 1: Pi heatmaps ---------------------------------------------
fig = figure('Name', 'Pi (selection frequency)', 'Color', 'w', ...
    'Position', [40 40 1500 760]);
for g = 1:nGroups
    subplot(1, nGroups, g);
    imagesc(C.Pi{g});
    colormap(gca, flipud(gray(256)));
    caxis([0 1]); colorbar;
    set(gca, 'XTick', 1:numel(muscleLabels), 'XTickLabel', muscleLabels, ...
        'XTickLabelRotation', 40, 'YTick', 1:numel(featureLabels), ...
        'YTickLabel', featureLabels, 'FontSize', 7, ...
        'TickLabelInterpreter', 'none');   % keep glut_max literal, not TeX
    title(sprintf('%s: Pi  (nnz consensus = %d, ambiguous = %d)', ...
        groupLabels{g}, nnz(C.mask{g}), C.nAmbiguous(g)));
    xlabel('muscle'); ylabel('feature');
end
save_fig(fig, figDir, 'lasso_stability_pi_heatmap');

%% ---- Figure 2: exact-mask frequency spectrum ---------------------------
fig = figure('Name', 'Exact-mask frequencies', 'Color', 'w', ...
    'Position', [80 80 1000 420]);
for g = 1:nGroups
    subplot(1, nGroups, g);
    cnt = C.exactMaskCounts{g}(:);
    kk = numel(cnt);
    bar(1:kk, cnt, 'FaceColor', [0.35 0.55 0.80]);
    hold on;
    plot([0.5 kk + 0.5], [C.modeCount(g) C.modeCount(g)], 'r--', 'LineWidth', 1.0);
    grid on;
    xlim([0.5 max(kk, 1) + 0.5]);
    xlabel('distinct exact masks (ranked)');
    ylabel('number of seeds');
    title(sprintf('%s: modal mask %d/%d = %.2f, mean Jaccard %.2f', ...
        groupLabels{g}, C.modeCount(g), C.nValid(g), C.modeFraction(g), ...
        C.jaccardMean(g)));
end
save_fig(fig, figDir, 'lasso_stability_mask_frequencies');

%% ---- Figure 3: Jaccard matrices ----------------------------------------
fig = figure('Name', 'Mask similarity', 'Color', 'w', ...
    'Position', [100 100 1000 460]);
for g = 1:nGroups
    subplot(1, nGroups, g);
    J = C.jaccardMatrix{g};
    imagesc(J, [0 1]);
    axis square; colorbar;
    colormap(gca, parula(256));
    title(sprintf('%s: Jaccard (dominant cluster %.2f, %d cluster(s))', ...
        groupLabels{g}, C.dominantClusterFraction(g), numel(C.clusterSizes{g})));
    xlabel('replicate'); ylabel('replicate');
end
save_fig(fig, figDir, 'lasso_stability_jaccard');

%% ---- Figure 4: support-size distribution -------------------------------
fig = figure('Name', 'Support sizes', 'Color', 'w', ...
    'Position', [120 120 900 420]);
for g = 1:nGroups
    subplot(1, nGroups, g);
    sz = C.supportSizes{g}(:);
    histogram(sz, 'BinMethod', 'integers', 'FaceColor', [0.55 0.70 0.85]);
    hold on;
    yl = ylim();
    plot([nnz(C.mask{g}) nnz(C.mask{g})], yl, 'r-', 'LineWidth', 1.6);
    grid on;
    xlabel('nnz of the replicate mask');
    ylabel('count');
    title(sprintf('%s: consensus nnz = %d (median replicate %g)', ...
        groupLabels{g}, nnz(C.mask{g}), median(sz)));
    legend({'replicates', 'consensus'}, 'Location', 'best');
end
save_fig(fig, figDir, 'lasso_stability_support_sizes');

%% ---- Figure 5: per-muscle RMSE, replicates vs frozen-mask refit --------
fig = figure('Name', 'RMSE comparison', 'Color', 'w', ...
    'Position', [140 140 1050 440]);
for g = 1:nGroups
    subplot(1, nGroups, g);
    [Rmat, valid] = collect_rmse(R, g);
    if isempty(Rmat)
        title(sprintf('%s: no replicate RMSE', groupLabels{g}));
        continue;
    end
    nM = size(Rmat, 2);
    hold on;
    for i = 1:size(Rmat, 1)
        plot(1:nM, Rmat(i, :), '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.6);
    end
    plot(1:nM, median(Rmat, 1), 'k-', 'LineWidth', 1.8);
    leg = {'replicates', 'replicate median'};
    if isfield(R, 'fitReference') && ~isempty(R.fitReference) && ...
            isfield(R.fitReference, 'RMSE_refit_full')
        plot(1:nM, R.fitReference.RMSE_refit_full(:).', 's-', ...
            'Color', [0.20 0.55 0.25], 'LineWidth', 1.4, 'MarkerSize', 4);
        leg{end + 1} = 'unfrozen full-data fit';
    end
    if isfield(R, 'fitFinal') && isfield(R.fitFinal, 'RMSE_refit_full')
        plot(1:nM, R.fitFinal.RMSE_refit_full(:).', 'r-o', 'LineWidth', 1.6, ...
            'MarkerSize', 4);
        leg{end + 1} = 'frozen-mask full-data refit';
    end
    legend(leg, 'Location', 'best');
    grid on;
    set(gca, 'XTick', 1:nM, 'XTickLabel', muscleLabels, 'XTickLabelRotation', 35, ...
        'TickLabelInterpreter', 'none');
    xlabel('muscle'); ylabel('weighted RMSE');
    title(sprintf('%s: fit error (%d valid replicates)', groupLabels{g}, nnz(valid)));
end
save_fig(fig, figDir, 'lasso_stability_rmse');

%% ---- Figure 6: Pi histogram -------------------------------------------
fig = figure('Name', 'Pi distribution', 'Color', 'w', ...
    'Position', [160 160 900 380]);
for g = 1:nGroups
    subplot(1, nGroups, g);
    edges = 0:0.05:1;
    histogram(C.Pi{g}(:), edges, 'FaceColor', [0.45 0.65 0.55]);
    hold on;
    yl = ylim();
    plot([0.5 0.5], yl, 'k--', 'LineWidth', 1.0);
    plot([0.7 0.7], yl, 'r:', 'LineWidth', 1.0);
    grid on;
    xlabel('\pi (selection frequency)'); ylabel('entries');
    title(sprintf('%s: stable-on %d | ambiguous %d | stable-off %d', ...
        groupLabels{g}, C.nStableOn(g), C.nAmbiguous(g), C.nStableOff(g)));
end
save_fig(fig, figDir, 'lasso_stability_pi_histogram');
end

%% ========================================================================
% Local helpers.
% ========================================================================
function labels = get_group_labels(R, nGroups)
labels = cell(1, nGroups);
if isfield(R.lasso, 'identification') && isfield(R.lasso.identification, 'groupLabels')
    labels = R.lasso.identification.groupLabels;
elseif isfield(R.options, 'fitOptions') && isfield(R.options.fitOptions, 'groupLabels')
    labels = R.options.fitOptions.groupLabels;
end
if numel(labels) ~= nGroups
    labels = arrayfun(@(i) sprintf('group%d', i), 1:nGroups, 'UniformOutput', false);
end
end

function [Rmat, valid] = collect_rmse(R, g)
valid = R.valid(:, g);
idx = find(valid);
Rmat = [];
if isempty(idx) || ~isfield(R, 'replicates')
    return;
end
first = R.replicates(idx(1)).rmse{g};
if isempty(first)
    return;
end
Rmat = nan(numel(idx), numel(first));
for i = 1:numel(idx)
    v = R.replicates(idx(i)).rmse{g};
    if numel(v) == size(Rmat, 2)
        Rmat(i, :) = v;
    end
end
end

function out = short_feature_names(names)
out = names(:);
for i = 1:numel(out)
    s = out{i};
    s = regexprep(s, '_r\.fiber_length_norm$', ' lCEN');
    s = regexprep(s, '_r\.mtu_force_norm$', ' fATN');
    s = regexprep(s, '_l\.fiber_length_norm$', ' lCEN');
    s = regexprep(s, '_l\.mtu_force_norm$', ' fATN');
    s = strrep(s, 'pelvis_tilt_u', 'pelvis tilt v');
    s = strrep(s, 'pelvis_tilt', 'pelvis tilt');
    s = strrep(s, 'hip_flexion_r', 'hip angle');
    s = strrep(s, 'hip_flexion_l', 'hip angle');
    s = strrep(s, 'knee_extension_r', 'knee angle');
    s = strrep(s, 'knee_extension_l', 'knee angle');
    s = strrep(s, 'ankle_dorsiflexion_r', 'ankle angle');
    s = strrep(s, 'ankle_dorsiflexion_l', 'ankle angle');
    s = strrep(s, '_u', ' v');
    out{i} = s;
end
end

function out = short_muscle_names(names)
out = names(:);
for i = 1:numel(out)
    out{i} = regexprep(out{i}, '_r$', '');
end
end

function save_fig(fig, figDir, baseName)
f = fullfile(figDir, [baseName '.png']);
try
    exportgraphics(fig, f, 'Resolution', 150);
catch
    print(fig, f, '-dpng', '-r150');
end
end
