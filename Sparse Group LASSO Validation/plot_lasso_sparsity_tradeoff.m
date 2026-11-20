function plot_lasso_sparsity_tradeoff(A, R, figFile)
%PLOT_LASSO_SPARSITY_TRADEOFF  Per-muscle lambda paths + summary sweeps.
%
%   plot_lasso_sparsity_tradeoff(A, R, figFile)
%
%   A  : struct from analyze_lasso_sparsity (results/lasso_sparsity_analysis.mat)
%   R  : OPTIONAL diagnostics struct (results/lasso_stability_diagnostics.mat
%        -> .R).  Only used for muscle names.
%   figFile : output PNG (default results/lasso_sparsity_tradeoff.png)
%
%   Layout (4 x 3)
%   ---------------
%   Rows 1-3 : one panel per muscle (9 panels).  Each panel shows that
%              muscle's stance (solid) and swing (dashed) lambda paths:
%              x = support size of that muscle-group,
%              y = RMSE / that muscle-group's best RMSE.
%              A pentagram marks each curve's own Kneedle knee.
%   Row 4    : consensus pi sweep, tolSelect sweep, and a legend/summary tile.
%
%   Rationale: lambda is selected independently for every muscle-group, so a
%   group-level median curve (the previous panel 1) is not a realizable
%   configuration; per-muscle panels are the honest view.
%
%   See also analyze_lasso_sparsity.

if nargin < 3 || isempty(figFile)
    here = fileparts(mfilename('fullpath'));
    figFile = fullfile(fileparts(here), 'results', 'lasso_sparsity_tradeoff.png');
end
if nargin < 2
    R = [];
end

nMus = size(A.path.nnzMat, 2);
if ~isempty(R) && isfield(R, 'fitReference') && isfield(R.fitReference, 'muscleNames')
    muscleNames = R.fitReference.muscleNames;
else
    muscleNames = {'hamstrings_r', 'bifemsh_r', 'glut_max_r', 'iliopsoas_r', ...
        'rect_fem_r', 'vasti_r', 'gastroc_r', 'soleus_r', 'tib_ant_r'};
end
shortName = strrep(muscleNames, '_r', '');
cols = {[0.90 0.55 0.10], [0.15 0.45 0.75]};
styles = {'-', '--'};

fig = figure('Name', 'LASSO sparsity (per muscle)', 'Color', 'w', ...
    'Position', [30 30 1500 1480]);
tl = tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

%% ---- 9 per-muscle panels ----------------------------------------------
kneeNnz = zeros(2, nMus);
for m = 1:nMus
    nexttile(m); hold on; grid on;
    for g = 1:2
        x = squeeze(A.path.nnzMat(g, m, :));
        y = squeeze(A.path.errMat(g, m, :));
        ok = isfinite(x) & isfinite(y) & x >= 1;   % nnz=0 is not a structure
        x = x(ok); y = y(ok);
        if numel(x) < 2
            continue;
        end
        % Pareto front first (line underneath), then every evaluated lambda
        % as a filled dot on top so the cloud is visible.
        [xf, yf] = pareto_front(x, y);
        plot(xf, yf, styles{g}, 'Color', cols{g}, 'LineWidth', 1.6);
        scatter(x, y, 24, cols{g}, 'filled', 'MarkerFaceAlpha', 0.55, ...
            'MarkerEdgeColor', 'none', 'HandleVisibility', 'off');
        [kx, ky] = local_knee(x, y);
        plot(kx, ky, 'p', 'MarkerSize', 10, 'MarkerFaceColor', cols{g}, ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
        kneeNnz(g, m) = kx;
    end
    xlim([0 25]); ylim([0.98 3.0]);
    title(sprintf('%s   (knee S%d/W%d)', shortName{m}, ...
        round(kneeNnz(1, m)), round(kneeNnz(2, m))), 'FontSize', 9, ...
        'Interpreter', 'none');
    if m > 6
        xlabel('support size (nonzeros)');
    end
    if mod(m - 1, 3) == 0
        ylabel('RMSE / best');
    end
end

%% ---- pi sweep ----------------------------------------------------------
nexttile(10); hold on; grid on;
piS = A.piSweep;
xK = piS.params(:);
h1 = plot(xK, piS.trainCV(:), '-o', 'LineWidth', 1.6);
h2 = plot(xK, piS.cvRMSE(:), '-s', 'LineWidth', 1.6);
kn = piS.knee;
plot(kn.params, kn.cvRMSE, 'p', 'MarkerSize', 14, 'LineWidth', 1.2, ...
    'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
xlabel('CMA-ES parameters (total)'); ylabel('mean RMSE');
title(sprintf('pi sweep (knee: %d params, pi>=%.2f)', ...
    round(kn.params), kn.pi), 'FontSize', 9);
legend([h1 h2], {'in-fold train', 'out-of-fold'}, 'Location', 'northeast', ...
    'FontSize', 7);

%% ---- tolSelect sweep ---------------------------------------------------
nexttile(11); hold on; grid on;
tS = A.tolSweep;
plot(tS.params(:), tS.trainRMSE(:), '-o', 'LineWidth', 1.6);
xlabel('CMA-ES parameters (total)'); ylabel('mean RMSE');
title('tolSelect sweep (full-data fit)', 'FontSize', 9);

%% ---- legend / summary tile --------------------------------------------
nexttile(12); axis off; hold on;
hS = plot(nan, nan, '-o', 'Color', cols{1}, 'LineWidth', 1.6, 'MarkerSize', 3);
hW = plot(nan, nan, '--o', 'Color', cols{2}, 'LineWidth', 1.6, 'MarkerSize', 3);
hK = plot(nan, nan, 'p', 'MarkerSize', 10, 'MarkerFaceColor', [0.6 0.6 0.6], ...
    'MarkerEdgeColor', 'k');
legend([hS hW hK], {'stance: dots = lambdas, solid = Pareto front', ...
    'swing: dots = lambdas, dashed = Pareto front', 'knee (per muscle-group, nnz>=1)'}, ...
    'Location', 'northwest', 'FontSize', 7);
totalKnee = sum(kneeNnz(:)) + 18;
text(0.02, 0.45, sprintf(['knee support: stance %d + swing %d nonzeros\n' ...
    '=> %d CMA-ES parameters (incl. 18 biases)'], ...
    sum(kneeNnz(1, :)), sum(kneeNnz(2, :)), totalKnee), ...
    'Units', 'normalized', 'FontSize', 9, 'VerticalAlignment', 'top');

title(tl, 'LASSO sparsity: per-muscle lambda paths (rows 1-3) and summary sweeps (row 4)', ...
    'FontSize', 11);

save_fig(fig, figFile);
fprintf('[sparsity-figure] wrote %s\n', figFile);
end

%% ========================================================================
function [xf, yf] = pareto_front(x, y)
%PARETO_FRONT  Non-dominated points (minimise both), sorted by x.
x = x(:); y = y(:);
n = numel(x);
nd = true(n, 1);
for i = 1:n
    for j = 1:n
        if j ~= i && x(j) <= x(i) && y(j) <= y(i) && (x(j) < x(i) || y(j) < y(i))
            nd(i) = false;
            break;
        end
    end
end
xf = x(nd); yf = y(nd);
[xf, is] = sort(xf); yf = yf(is);
end

%% ========================================================================
function [kx, ky] = local_knee(x, y)
%LOCAL_KNEE  Kneedle knee on one muscle-group's path (Pareto-filtered first).
ok = isfinite(x) & isfinite(y);
x = x(ok); y = y(ok);
n = numel(x);
if n < 3
    kx = x(end); ky = y(end);
    return;
end
nd = true(n, 1);
for i = 1:n
    for j = 1:n
        if j ~= i && x(j) <= x(i) && y(j) <= y(i) && (x(j) < x(i) || y(j) < y(i))
            nd(i) = false;
            break;
        end
    end
end
x = x(nd); y = y(nd);
[x, is] = sort(x); y = y(is);
xn = (x - min(x)) / max(max(x) - min(x), eps);
yn = (y - min(y)) / max(max(y) - min(y), eps);
d = (1 - xn) - yn;
[~, i] = max(d);
kx = x(i); ky = y(i);
end

%% ========================================================================
function save_fig(fig, f)
try
    exportgraphics(fig, f, 'Resolution', 150);
catch
    print(fig, f, '-dpng', '-r150');
end
end
