function visualize_preprocessing()
% Plot raw versus preprocessed GRF and sEMG over the 30 s analysis window.
% Run this file inside a JW_* trial folder; all paths and PNG outputs are local.

here = fileparts(mfilename('fullpath'));
summary = jsondecode(fileread(fullfile(here, 'steady_state_summary.json')));
win = double(summary.analysis_window(:).');

rawGrfPath = pickFile(here, '*_grf.mot');
[~, stem] = fileparts(rawGrfPath);
stem = erase(stem, '_grf');
[tRaw, grfRaw, labels] = readMot(rawGrfPath);
[tProc, grfProc, ~] = readMot(fullfile(here, [stem '_grf_filtered.mot']));

plotGrfRawFiltered(tRaw, grfRaw, tProc, grfProc, labels, win, stem, ...
    fullfile(here, [stem '_grf_raw_vs_filtered.png']));

semgPath = pickFile(here, '*_sEMG.csv');
[fs, ts, rawSemg] = readSemg(semgPath, size(grfRaw, 1));
sel = ts >= max(ts(1), win(1) - 1) & ts <= min(ts(end), win(2) + 1);
t = ts(sel);
x = rawSemg(sel, :);
[clean, replaced] = despikeSemg(x);
[hpB, hpA] = butter(2, 10 / (fs / 2), 'high');
[lpB, lpA] = butter(2, 5 / (fs / 2), 'low');
env = abs(filtfilt(hpB, hpA, clean));
env = filtfilt(lpB, lpA, env);
env = max(env, 0);
inner = t >= win(1) & t <= win(2);
scale = max(env(inner, :), [], 1);
scale(scale <= 0) = 1;
envNorm = env ./ scale;
plotSemg(t, x, clean, replaced, envNorm, win, stem, ...
    fullfile(here, [stem '_semg_raw_spikes_filtered.png']));

fprintf('[%s] wrote two preprocessing figures for %.2f-%.2f s\n', ...
    stem, win(1), win(2));
end

function plotGrfRawFiltered(tRaw, raw, tProc, proc, labels, win, stem, outPath)
componentLabels = {{'vx', 'vy', 'vz'}, {'tx', 'ty', 'tz'}, {'px', 'py', 'pz'}};
rowNames = {'Force', 'Torque', 'COP'};
units = {'N', 'N m', 'm'};
offsets = {[1 2 3], [7 8 9], [4 5 6]};
expected = {};
for b = 1:2
    for suffix = {'vx', 'vy', 'vz', 'px', 'py', 'pz'}
        expected{end + 1} = sprintf('ground_force%d_%s', b, suffix{1}); %#ok<AGROW>
    end
    for suffix = {'x', 'y', 'z'}
        expected{end + 1} = sprintf('ground_torque%d_%s', b, suffix{1}); %#ok<AGROW>
    end
end
columns = zeros(1, numel(expected));
for k = 1:numel(expected)
    columns(k) = find(strcmp(labels, expected{k}), 1);
    if columns(k) == 0
        error('Missing GRF column %s', expected{k});
    end
end
selRaw = tRaw >= win(1) & tRaw <= win(2);
selProc = tProc >= win(1) & tProc <= win(2);
fig = figure('Color', 'w', 'Position', [60 40 1500 900]);
for r = 1:3
    for b = 1:2
        ax = subplot(3, 2, (r - 1) * 2 + b, 'Parent', fig);
        hold(ax, 'on');
        handles = gobjects(1, 3);
        legendLabels = cell(1, 3);
        for k = 1:3
            col = columns((b - 1) * 9 + offsets{r}(k));
            plot(ax, tRaw(selRaw), raw(selRaw, col), '-', ...
                'Color', [0.74 0.74 0.74], 'LineWidth', 1.1, ...
                'HandleVisibility', 'off');
            handles(k) = plot(ax, tProc(selProc), proc(selProc, col), '-', ...
                'LineWidth', 1.1);
            legendLabels{k} = sprintf('%s (%s)', componentLabels{r}{k}, units{r});
        end
        hold(ax, 'off');
        grid(ax, 'on'); box(ax, 'on'); xlim(ax, win);
        xlabel(ax, 'time (s)');
        ylabel(ax, sprintf('%s (%s)', rowNames{r}, units{r}));
        title(ax, sprintf('belt%d -- %s', b, rowNames{r}), 'FontWeight', 'normal');
        if r == 1
            legend(ax, handles, legendLabels, 'Location', 'best');
        end
    end
end
sgtitle(sprintf('%s | GRF raw vs filtered | %.2f-%.2f s', stem, win(1), win(2)), ...
    'Interpreter', 'none');
drawnow;
print(fig, outPath, '-dpng', '-r170');
set(fig, 'Visible', 'on');
fprintf('wrote %s\n', outPath);
end

function plotSemg(t, raw, clean, replaced, envNorm, win, stem, outPath)
muscles = {'rect_fem_r', 'vasti_r', 'hamstrings_r', 'tib_ant_r', 'gastroc_r', ...
    'rect_fem_l', 'vasti_l', 'hamstrings_l', 'tib_ant_l', 'gastroc_l'};
shortNames = {'RF_R', 'VL_R', 'BF_R', 'TA_R', 'GAS_R', ...
    'RF_L', 'VL_L', 'BF_L', 'TA_L', 'GAS_L'};
plotMask = t >= win(1) & t <= win(2);
idx = find(plotMask);
idx = idx(1:5:end);
fig = figure('Color', 'w', 'Position', [40 20 1600 1800]);
for i = 1:10
    subplot(10, 2, 2 * i - 1);
    hRaw = plot(t(idx), raw(idx, i), 'Color', [0.60 0.60 0.60]); hold on;
    hClean = plot(t(idx), clean(idx, i), 'Color', [0.10 0.35 0.75]);
    legendHandles = [hRaw hClean];
    legendLabels = {'raw', 'despiked'};
    rep = find(plotMask & replaced(:, i));
    if ~isempty(rep)
        hSpikes = plot(t(rep), raw(rep, i), 'r.', 'MarkerSize', 8);
        legendHandles(end + 1) = hSpikes;
        legendLabels{end + 1} = 'flagged spikes';
    end
    grid on; box off; xlim(win);
    title(sprintf('EMG%d %s (%d flagged samples)', i + 6, muscles{i}, ...
        numel(rep)), 'Interpreter', 'none', 'FontWeight', 'normal');
    ylabel({shortNames{i}, 'raw (V)'});
    if i == 1
        legend(legendHandles, legendLabels, 'Location', 'northeast');
    end

    subplot(10, 2, 2 * i);
    plot(t(idx), envNorm(idx, i), 'Color', [0.85 0.25 0.20], 'LineWidth', 0.8);
    grid on; box off; xlim(win); ylim([0 1.05]);
    title('filtered envelope: 10 Hz HP, rectify, 5 Hz LP, max-normalised', ...
        'FontWeight', 'normal');
    ylabel({shortNames{i}, 'norm.'});
    if i == 10
        subplot(10, 2, 2 * i - 1); xlabel('time (s)');
        subplot(10, 2, 2 * i); xlabel('time (s)');
    end
end
sgtitle(sprintf('%s | sEMG raw, flagged spikes and filtered envelope | %.2f-%.2f s', ...
    stem, win(1), win(2)), 'Interpreter', 'none');
drawnow;
print(fig, outPath, '-dpng', '-r170');
set(fig, 'Visible', 'on');
fprintf('wrote %s\n', outPath);
end

function [t, data, labels] = readMot(path)
fid = fopen(path, 'r');
if fid < 0
    error('Cannot open %s', path);
end
cleanup = onCleanup(@() fclose(fid));
while ~feof(fid)
    line = fgetl(fid);
    if strcmp(strtrim(line), 'endheader')
        break;
    end
end
labels = strsplit(strtrim(fgetl(fid)));
fmt = repmat('%f', 1, numel(labels));
values = textscan(fid, fmt, 'CollectOutput', true);
data = values{1};
t = data(:, 1);
end

function [fs, t, x] = readSemg(path, nRows)
fid = fopen(path, 'r', 'n', 'UTF-8');
if fid < 0
    error('Cannot open %s', path);
end
cleanup = onCleanup(@() fclose(fid));
fgetl(fid);
fs = str2double(fgetl(fid));
fgetl(fid);
header = strsplit(strtrim(fgetl(fid)), ',');
fgetl(fid);
indices = zeros(1, 10);
for i = 1:10
    indices(i) = find(strcmp(header, sprintf('EMG%d', i + 6)), 1);
    if indices(i) == 0
        error('EMG%d not found in %s', i + 6, path);
    end
end
fmt = '';
for k = 1:max(indices)
    if k <= 2 || any(indices == k)
        fmt = [fmt '%f']; %#ok<AGROW>
    else
        fmt = [fmt '%*f']; %#ok<AGROW>
    end
end
fmt = [fmt '%*[^\n]'];
values = textscan(fid, fmt, nRows, 'Delimiter', ',', 'EmptyValue', NaN);
frame = values{1};
subFrame = values{2};
if numel(frame) ~= nRows
    error('Expected %d samples in the 2 kHz sEMG block, got %d', ...
        nRows, numel(frame));
end
nSub = max(subFrame) + 1;
t = ((frame - 1) * nSub + subFrame) / fs;
if any(diff(t) <= 0)
    error('sEMG time base is not strictly increasing');
end
x = zeros(nRows, 10);
for i = 1:10
    x(:, i) = values{i + 2};
end
end

function [clean, replaced] = despikeSemg(x)
[n, nCh] = size(x);
clean = x;
replaced = false(n, nCh);
for i = 1:nCh
    ax = abs(x(:, i));
    threshold = 2.5 * sortedPercentile(ax, 99.9);
    flagged = ax > threshold;
    if ~any(flagged)
        continue;
    end
    grown = conv(double(flagged), ones(21, 1), 'same') > 0;
    replaced(:, i) = grown;
    d = diff([0; grown; 0]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;
    for k = 1:numel(starts)
        s = starts(k);
        e = ends(k);
        sup = (max(1, s - 50):min(n, e + 50)).';
        sup = sup(~grown(sup));
        if numel(sup) < 4
            sup = find(~grown);
        end
        if numel(sup) < 2
            clean(s:e, i) = 0;
        else
            clean(s:e, i) = interp1(sup, x(sup, i), (s:e).', 'pchip');
        end
    end
end
end

function value = sortedPercentile(x, p)
xs = sort(x(:));
k = min(numel(xs), max(1, ceil(p / 100 * numel(xs))));
value = xs(k);
end

function path = pickFile(folder, pattern)
hits = dir(fullfile(folder, pattern));
if isempty(hits)
    error('No file matching %s in %s', pattern, folder);
end
path = fullfile(folder, hits(1).name);
end