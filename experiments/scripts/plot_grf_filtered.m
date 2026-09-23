%% plot_grf_filtered.m -- raw vs processed ground reaction data
%
% One subplot per GRF component: three rows (force / torque / COP) and two
% columns (belt1 = left plate, belt2 = right plate), 12 curves in all, each
% drawn twice -- raw underneath, processed on top.
%
%   raw       : level_walking_grf.mot
%   processed : level_walking_grf_filtered.mot
%
% What the processing did, per belt (see run_ik_id.py -> filter_grf()):
%   * vx, vy, vz, ty : 15 Hz zero-phase low-pass, applied to every stance phase
%                      separately with 0.3 s of mirror padding so the filter
%                      never sees the load/unload step; exactly zero during
%                      swing; 25 ms cosine taper at each stance end
%   * tx, tz         : exactly zero (the raw columns carry only ~1e-5 Nm
%                      rounding noise; only ty, the free moment about the
%                      vertical axis, holds signal)
%   * px, py, pz     : UNFILTERED, bit-identical to the raw recording.  The
%                      swing-to-stance COP step is a real measurement (the
%                      foot landing where it lands) and filtering it injected
%                      up to 31 Nm of spurious moment into the loading
%                      response, so it is deliberately left alone.
%
% Usage:
%   plot_grf_filtered                    % 170-240 s, figure only
%   plot_grf_filtered([180 200])         % custom window
%   plot_grf_filtered([170 240], true)   % also save a PNG next to this file
%
% The numeric companion is experiments/scripts/check_grf_cop.py, which
% prints the same comparison as pass/fail checks.

function plot_grf_filtered(tlim, save_figs)

if nargin < 1 || isempty(tlim)
    tlim = [170 240];           % s
end
if nargin < 2 || isempty(save_figs)
    save_figs = false;
end
tlim = tlim(:)';
assert(numel(tlim) == 2 && tlim(1) < tlim(2), 'tlim must be [t0 t1]');

here     = fileparts(mfilename('fullpath'));
data_dir = fullfile(here, '..', 'data', 'SQR_walking');
f_raw    = fullfile(data_dir, 'level_walking_grf.mot');
f_proc   = fullfile(data_dir, 'level_walking_grf_filtered.mot');

assert(exist(f_raw, 'file') == 2, 'raw GRF file not found: %s', f_raw);
assert(exist(f_proc, 'file') == 2, ...
    ['processed GRF file not found: %s\n' ...
     'Run:  run_ik_id.py   (it always regenerates the full pipeline)'], f_proc);

[tr, raw,  nr] = read_opensim_mot(f_raw);
[tp, proc, np] = read_opensim_mot(f_proc);

fprintf('raw       : %s  (%d rows, %.1f-%.1f s)\n', f_raw,  numel(tr), tr(1), tr(end));
fprintf('processed : %s  (%d rows, %.1f-%.1f s)\n', f_proc, numel(tp), tp(1), tp(end));
assert(isequal(nr, np), 'column names differ between raw and processed');
assert(numel(tr) == numel(tp) && max(abs(tr(:) - tp(:))) < 1e-9, ...
    'raw and processed are not on the same time base');

% --- column lookup -------------------------------------------------------
% header: time | ground_force1_vx vy vz | ground_force1_px py pz |
%               ground_torque1_x y z | ground_force2_... (same 9)
expect = {};
for b = 1:2
    for s = {'vx', 'vy', 'vz', 'px', 'py', 'pz'}
        expect{end+1} = sprintf('ground_force%d_%s', b, s{1});   %#ok<AGROW>
    end
    for s = {'x', 'y', 'z'}
        expect{end+1} = sprintf('ground_torque%d_%s', b, s{1});  %#ok<AGROW>
    end
end

idx = zeros(1, numel(expect));
for k = 1:numel(expect)
    j = find(strcmp(nr, expect{k}), 1);
    assert(~isempty(j), 'column "%s" not found in %s', expect{k}, f_raw);
    idx(k) = j;
end

% --- window --------------------------------------------------------------
selr = tr >= tlim(1) & tr <= tlim(2);
selp = tp >= tlim(1) & tp <= tlim(2);
assert(any(selr), 'no raw samples inside [%g %g] s', tlim(1), tlim(2));

% --- layout: row = quantity, column = belt -------------------------------
units  = {'N', 'N\cdotm', 'm'};
labels = {{'vx','vy','vz'}, {'tx','ty','tz'}, {'px','py','pz'}};
rowname = {'Force', 'Torque', 'COP'};
% offsets of vx,vy,vz / tx,ty,tz / px,py,pz inside a 9-column belt block
offs   = {[1 2 3], [7 8 9], [4 5 6]};

fig = figure('Name', sprintf('GRF raw vs processed  (%.0f-%.0f s)', tlim(1), tlim(2)), ...
             'Color', 'w', 'Position', [60 40 1500 900]);
if exist('tiledlayout', 'file') == 2 || exist('tiledlayout', 'builtin') == 5
    tl = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf(['Ground reaction data -- raw vs processed, %.0f-%.0f s' ...
          '   (processed: 15 Hz per-stance, COP untouched)'], tlim(1), tlim(2)), ...
          'FontWeight', 'bold', 'FontSize', 13);
    use_tiles = true;
else
    use_tiles = false;                      % pre-R2019b fallback
end

fprintf('\n%-34s %11s %11s %9s\n', 'channel', 'raw peak', 'proc peak', 'kept %');
fprintf('%s\n', repmat('-', 1, 68));

for r = 1:3
    for b = 1:2
        if use_tiles
            ax = nexttile(tl, (r - 1) * 2 + b);
        else
            ax = subplot(3, 2, (r - 1) * 2 + b, 'Parent', fig);
        end
        hold(ax, 'on');

        hs = gobjects(1, 3);
        lg = cell(1, 3);
        for k = 1:3
            j = idx((b - 1) * 9 + offs{r}(k));

            plot(ax, tr(selr), raw(selr, j), '-', 'Color', [0.74 0.74 0.74], ...
                 'LineWidth', 1.1, 'HandleVisibility', 'off');
            hs(k) = plot(ax, tp(selp), proc(selp, j), '-', 'LineWidth', 1.1);
            lg{k} = sprintf('%s (%s)', labels{r}{k}, units{r});

            a = max(abs(raw(selr, j)));
            c = max(abs(proc(selp, j)));
            if a > 0
                kept = sprintf('%9.1f', 100 * c / a);
            else
                kept = '      n/a';
            end
            fprintf('belt%d %-27s %11.4g %11.4g %s\n', ...
                    b, sprintf('%s.%s', rowname{r}, labels{r}{k}), a, c, kept);
        end
        hold(ax, 'off');

        grid(ax, 'on');
        box(ax, 'on');
        xlim(ax, tlim);
        xlabel(ax, 'time (s)');
        ylabel(ax, sprintf('%s (%s)', rowname{r}, units{r}));
        title(ax, sprintf('belt%d -- %s', b, rowname{r}), 'FontWeight', 'normal');
        if r == 1
            legend(ax, hs, lg, 'Location', 'best');
        end
    end
end

if save_figs
    out = fullfile(here, sprintf('grf_raw_vs_processed_%g_%g.png', tlim(1), tlim(2)));
    try
        exportgraphics(fig, out, 'Resolution', 150);
    catch
        print(fig, out, '-dpng', '-r150');   % pre-R2020a fallback
    end
    fprintf('\nsaved figure -> %s\n', out);
    close(fig);
end

end % plot_grf_filtered


% -------------------------------------------------------------------------
function [t, data, names] = read_opensim_mot(path)
%READ_OPENSIM_MOT Read an OpenSim .mot/.sto file (header, labels, numbers).

fid = fopen(path, 'r');
assert(fid > 0, 'cannot open %s', path);
cleanup = onCleanup(@() fclose(fid));

while true
    ln = fgetl(fid);
    assert(ischar(ln), 'no endheader in %s', path);
    if strncmpi(strtrim(ln), 'endheader', 9)
        break;
    end
end

% the first non-empty line after endheader is the column-label line
names = {};
while true
    ln = fgetl(fid);
    assert(ischar(ln), 'no column labels in %s', path);
    if ~isempty(strtrim(ln))
        names = strsplit(strtrim(ln));
        break;
    end
end

data = textscan(fid, repmat('%f', 1, numel(names)), 'CollectOutput', true);
data = data{1};
t    = data(:, 1);

assert(size(data, 1) == numel(t), ...
    'truncated or ragged data in %s: read %d rows', path, size(data, 1));
assert(size(data, 2) == numel(names), ...
    'expected %d columns in %s, found %d', numel(names), path, size(data, 2));
end
