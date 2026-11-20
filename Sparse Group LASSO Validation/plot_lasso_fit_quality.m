function plot_lasso_fit_quality(diagFile, dsMatFile, outDir, zoomWin, tag)
%PLOT_LASSO_FIT_QUALITY  Per-muscle measured-vs-predicted excitation curves.
%
%   plot_lasso_fit_quality(diagFile, dsMatFile, outDir, zoomWin)
%
%   RMSE numbers are hard to interpret; this renders the actual excitation
%   traces for the 9 right-leg muscles so the fit quality can be judged by eye.
%
%   Outputs (PNG, in outDir)
%   ------------------------
%   lasso_fit_quality.png        full 30 s record, measured vs the majority
%                                consensus and (if available) the pi>=0.7 one
%   lasso_fit_quality_zoom.png   a 3 s window (default centred in the record)
%                                so single gait cycles are visible
%   lasso_fit_rmse_comparison.png  per-muscle RMSE bars: unfrozen full fit vs
%                                majority consensus vs pi>=0.7 consensus
%
%   Inputs
%   ------
%   diagFile  : results/lasso_stability_diagnostics.mat (majority consensus)
%   dsMatFile : .../moco_dataset/moco_lasso_dataset.mat (time vector)
%   outDir    : output folder (default results/)
%   zoomWin   : [t0 t1] zoom window in seconds (default: 3 s, 1/3 into record)
%
%   See also lasso_stability_from_design, analyze_lasso_sparsity.

if nargin < 1 || isempty(diagFile)
    here = fileparts(mfilename('fullpath'));
    diagFile = fullfile(fileparts(here), 'results', 'lasso_stability_diagnostics.mat');
end
if nargin < 2 || isempty(dsMatFile)
    here = fileparts(mfilename('fullpath'));
    dsMatFile = fullfile(here, 'moco_dataset', 'moco_lasso_dataset.mat');
end
here = fileparts(mfilename('fullpath'));
repoRoot = fileparts(here);
if nargin < 3 || isempty(outDir)
    outDir = fullfile(repoRoot, 'results');
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

D = load(diagFile, 'R');
R = D.R;
fit = R.fitFinal;
ds = load(dsMatFile, 'ds');
ds = ds.ds;

Y = fit.Y;
YP = fit.YPred_refit;
stanceMask = fit.phase.stanceMask;
nY = size(Y, 1);
delaySteps = fit.options.delaySteps;
t = ds.time(:);
if numel(t) == nY
    tFit = t;
elseif numel(t) - delaySteps == nY
    tFit = t((delaySteps + 1):end);
else
    tFit = linspace(t(1), t(end), nY).';
end

if nargin < 4 || isempty(zoomWin)
    T = tFit(end) - tFit(1);
    zoomWin = [tFit(1) + 0.33 * T, tFit(1) + 0.33 * T + 3];
end
if nargin < 5 || isempty(tag)
    tag = '';
end
sfx = tag;

muscles = fit.muscleNames;
nMus = numel(muscles);

% optional second (sparser) consensus for a side-by-side visual comparison
YP70 = [];
len70 = NaN;
pi70File = fullfile(fileparts(diagFile), 'lasso_stability_diagnostics_pi70.mat');
if isfile(pi70File)
    R70 = load(pi70File, 'R');
    if isfield(R70.R, 'totalLen') && R70.R.totalLen ~= R.totalLen && ...
            isfield(R70.R, 'fitFinal') && isfield(R70.R.fitFinal, 'YPred_refit')
        YP70 = R70.R.fitFinal.YPred_refit;
        len70 = R70.R.totalLen;
    end
end

%% ---- Figure 1: full record ---------------------------------------------
fig = figure('Name', 'LASSO fit quality (full record)', 'Color', 'w', ...
    'Position', [40 40 1500 900]);
draw_grid(tFit, Y, YP, YP70, muscles, R.totalLen, len70, [tFit(1) tFit(end)], stanceMask);
sgtitle(sprintf('Measured vs predicted excitation, full 30 s (black = measured, %d vs %d parameters; bg yellow = stance, blue = swing)', ...
    R.totalLen, len70), 'FontSize', 11);
f1 = fullfile(outDir, ['lasso_fit_quality' sfx '.png']);
save_fig(fig, f1);

%% ---- Figure 2: 3 s zoom -------------------------------------------------
fig = figure('Name', 'LASSO fit quality (zoom)', 'Color', 'w', ...
    'Position', [40 40 1500 900]);
draw_grid(tFit, Y, YP, YP70, muscles, R.totalLen, len70, zoomWin, stanceMask);
if isempty(YP70)
    sgtitle(sprintf('Zoom %.1f-%.1f s: one gait cycle structure (black = measured, %d parameters; bg yellow = stance, blue = swing)', ...
        zoomWin(1), zoomWin(2), R.totalLen), 'FontSize', 11);
else
    sgtitle(sprintf('Zoom %.1f-%.1f s: one gait cycle structure (black = measured, %d vs %d parameters; bg yellow = stance, blue = swing)', ...
        zoomWin(1), zoomWin(2), R.totalLen, len70), 'FontSize', 11);
end
f2 = fullfile(outDir, ['lasso_fit_quality_zoom' sfx '.png']);
save_fig(fig, f2);

%% ---- Figure 3: per-muscle RMSE, three structures ----------------------
fig = figure('Name', 'Per-muscle fit error', 'Color', 'w', ...
    'Position', [80 80 1200 420]);
hold on; grid on;
M = [];
leg = {};
if isfield(R, 'fitReference') && ~isempty(R.fitReference)
    M(:, end + 1) = R.fitReference.RMSE_refit_full(:);
    leg{end + 1} = sprintf('unfrozen full fit (%d)', R.referenceTotalLen);
end
M(:, end + 1) = fit.RMSE_refit_full(:);
leg{end + 1} = sprintf('frozen mask (%d)', R.totalLen);
if ~isempty(YP70)
    M(:, end + 1) = R70.R.fitFinal.RMSE_refit_full(:);
    leg{end + 1} = sprintf('pi>=0.7 consensus (%d)', len70);
end
bar(M);
set(gca, 'XTick', 1:nMus, 'XTickLabel', strrep(muscles, '_', '\_'), ...
    'XTickLabelRotation', 35);
ylabel('weighted RMSE (excitation units)');
legend(leg, 'Location', 'northwest');
title('Per-muscle fit error: does freezing / sparsifying cost much?');
f3 = fullfile(outDir, ['lasso_fit_rmse_comparison' sfx '.png']);
save_fig(fig, f3);

fprintf('[fit-quality] wrote:\n  %s\n  %s\n  %s\n', f1, f2, f3);
end

%% ========================================================================
function draw_grid(tFit, Y, YP, YP70, muscles, len1, len70, xl, stanceMask)
nMus = numel(muscles);
for m = 1:nMus
    subplot(3, 3, m);
    hold on; grid on;
    ylim([0 1.05]);
    if ~isempty(stanceMask)
        shade_phase(tFit, stanceMask, [1.00 0.98 0.84]);   % very pale yellow
        shade_phase(tFit, ~stanceMask, [0.87 0.94 1.00]);  % very pale blue
    end
    hY = plot(tFit, Y(:, m), 'k-', 'LineWidth', 1.2);
    hP = plot(tFit, min(YP(:, m), 1), '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.1);
    hP70 = gobjects(1, 1);
    r1 = rmse_simple(Y(:, m), YP(:, m));
    if ~isempty(YP70)
        hP70 = plot(tFit, min(YP70(:, m), 1), '--', 'Color', [0.15 0.35 0.75], 'LineWidth', 1.0);
        r2 = rmse_simple(Y(:, m), YP70(:, m));
        title(sprintf('%s   RMSE %.4f (%d) / %.4f (%d)', ...
            strrep(muscles{m}, '_', '\_'), r1, len1, r2, len70), 'FontSize', 8);
    else
        title(sprintf('%s   RMSE=%.4f   R^2=%.2f', ...
            strrep(muscles{m}, '_', '\_'), r1, rsq(Y(:, m), YP(:, m))), 'FontSize', 9);
    end
    xlim(xl);
    ylim([0 1.05]);
    if m > 6
        xlabel('t (s)');
    end
    if mod(m - 1, 3) == 0
        ylabel('excitation');
    end
    if m == 1
        if isempty(YP70)
            legend([hY hP], {'measured', 'prediction'}, 'Location', 'northwest', ...
                'FontSize', 7);
        else
            legend([hY hP hP70], {'measured', sprintf('consensus %d', len1), ...
                sprintf('consensus %d', len70)}, 'Location', 'northwest', 'FontSize', 7);
        end
    end
end
end

%% ========================================================================
function shade_phase(t, mask, color)
%SHADE_PHASE  Paint contiguous TRUE runs of mask with a pale background.
mask = logical(mask(:));
if numel(mask) ~= numel(t)
    return;
end
yl = ylim();
d = diff([false; mask; false]);
starts = find(d == 1);
finish = find(d == -1) - 1;
for r = 1:numel(starts)
    a = starts(r); b = finish(r);
    if a >= 1 && b <= numel(t)
        patch([t(a), t(b), t(b), t(a)], ...
            [yl(1), yl(1), yl(2), yl(2)], color, ...
            'EdgeColor', 'none', 'FaceAlpha', 0.55, 'HandleVisibility', 'off');
    end
end
end

%% ========================================================================
function r = rmse_simple(y, yh)
r = sqrt(mean((y(:) - yh(:)).^2));
end

function r2 = rsq(y, yh)
y = y(:); yh = yh(:);
ssRes = sum((y - yh).^2);
ssTot = sum((y - mean(y)).^2);
r2 = 1 - ssRes / max(ssTot, eps);
end

function save_fig(fig, f)
try
    exportgraphics(fig, f, 'Resolution', 150);
catch
    print(fig, f, '-dpng', '-r150');
end
end
