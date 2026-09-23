%% process_semg.m
% =========================================================================
% sEMG 信号处理程序 (SQR level walking)
%
% 处理流程:
%   0) 去伪迹(despike): 剔除异常尖峰并插值修复
%      判据: |x| > kP * P99.9(该通道本身), kP = 2.5
%      修复: 被标记区间(外扩 ±5 ms)用保形三次插值(PCHIP)重建
%      (详见 diagnose_semg_artifacts.m 的诊断结论; 本步骤可置 enableDespike=false 关闭)
%   1) 10 Hz 二阶 Butterworth 高通滤波  (去除直流/基线漂移)
%   2) 全波整流 abs(x)
%   3) 5 Hz 二阶 Butterworth 低通滤波   (得到平滑包络)
%   4) 按各肌肉包络的最大值归一化 (0 ~ 1)
%
% 数据来源: experiments\data\SQR_walking\level_walking_semg.csv
%   (Delsys EMGworks 导出的多设备 CSV, 本程序只读取第一段 EMG 数据块,
%    采样率 2000 Hz, 每帧 20 个子样本)
%
% 通道对应关系:
%   EMG7  - 右侧股直肌   Right Rectus Femoris
%   EMG8  - 右侧股外侧肌 Right Vastus Lateralis (Vasti Externus)
%   EMG9  - 右侧股二头肌 Right Biceps Femoris
%   EMG10 - 右侧胫骨前肌 Right Tibialis Anterior
%   EMG11 - 右侧腓肠肌   Right Gastrocnemius
%   EMG12 - 左侧股直肌   Left Rectus Femoris
%   EMG13 - 左侧股外侧肌 Left Vastus Lateralis (Vasti Externus)
%   EMG14 - 左侧股二头肌 Left Biceps Femoris
%   EMG15 - 左侧胫骨前肌 Left Tibialis Anterior
%   EMG16 - 左侧腓肠肌   Left Gastrocnemius
%
% 输出 (experiments\moco\outputs\semg\):
%   level_walking_semg_raw_and_processed.png/.fig  各肌肉 raw(去伪迹后) 与处理后对比图
%   level_walking_semg_despike_effect.png          去伪迹前后归一化包络对比图
%   level_walking_semg_processed.csv               归一化包络 (Time_s, EMG7~EMG16, 全时程)
%   level_walking_semg_processed.mat               数据 (t, emgRaw, emgClean, emgEnvNorm ...)
%   semg_artifact_report.csv                       每个被替换区间的清单(可核查)
%   注: 图只显示 plotWin 区间(默认 170~240 s); CSV/MAT 仍是全时程 0~249.6 s 数据。
%
% 依赖: Signal Processing Toolbox (butter, filtfilt); 建议 MATLAB R2020b+
% =========================================================================

clear; close all; clc;

%% 1. 参数设置 -----------------------------------------------------------
scriptDir = fileparts(mfilename('fullpath'));   % 脚本所在目录 experiments\moco
semgFile  = fullfile(scriptDir, '..', 'data', 'SQR_walking', 'level_walking_semg.csv');
outDir    = fullfile(scriptDir, 'outputs', 'semg');
if ~exist(outDir, 'dir'); mkdir(outDir); end

% --- 去伪迹参数 (见 diagnose_semg_artifacts.m) ---
enableDespike    = true;  % 是否执行去伪迹
kP               = 2.5;   % 判据: |x| > kP * P99.9(该通道)
despikeMargin    = 10;    % 标记区间向两侧外扩的样本数 (10 = 5 ms @2000 Hz)
despikeNeighbors = 50;    % 插值支撑点: 区间两侧各取多少个干净样本 (50 = 25 ms)

% --- 滤波与归一化参数 ---
hpFc       = 10;        % 高通截止频率 (Hz)
lpFc       = 5;         % 低通截止频率 (Hz)
filtOrder  = 2;         % Butterworth 阶数
zeroPhase  = true;      % true : filtfilt 零相位滤波(正反各2阶, 等效4阶, 无相位延迟, EMG 处理常用)
                        % false: filter() 严格单次2阶前向滤波(有相位延迟)
plotWin      = [170 240]; % 显示(绘图)区间 (s): 只呈现这一段时程, 便于看清细节
plotDecimate = 1;         % 绘图降采样因子(仅加快绘图, 不影响处理结果; 1 = 全部点)

% 通道与肌肉对应关系 (与文件列 EMG7~EMG16 一一对应)
channels  = 7:16;
muscleNames = { ...
    'Right Rectus Femoris (右侧股直肌)', ...
    'Right Vastus Lateralis (右侧股外侧肌)', ...
    'Right Biceps Femoris (右侧股二头肌)', ...
    'Right Tibialis Anterior (右侧胫骨前肌)', ...
    'Right Gastrocnemius (右侧腓肠肌)', ...
    'Left Rectus Femoris (左侧股直肌)', ...
    'Left Vastus Lateralis (左侧股外侧肌)', ...
    'Left Biceps Femoris (左侧股二头肌)', ...
    'Left Tibialis Anterior (左侧胫骨前肌)', ...
    'Left Gastrocnemius (左侧腓肠肌)'};
shortNames = {'RF_R','VL_R','BF_R','TA_R','GAS_R', ...
              'RF_L','VL_L','BF_L','TA_L','GAS_L'};

%% 2. 读取数据 (只取 CSV 第一段: EMG 数据块) -----------------------------
fprintf('读取数据: %s\n', semgFile);
raw = fileread(semgFile);
lines = regexp(raw, '\r?\n', 'split')';         % 全部行
lines = lines(:);

fs = str2double(strtrim(lines{2}));             % 第2行为采样率 (2000 Hz)
if ~contains(lines{4}, 'EMG')
    warning('第4行列名中未找到 EMG, 请确认数据文件是否为 EMG 导出');
end

% 第一段数据块在第一个空行之前结束 (其后是测力台/动捕等其他设备的数据)
blankIdx  = find(cellfun(@(s) isempty(strtrim(s)), lines(6:end)), 1, 'first') + 5;
dataLines = lines(6:blankIdx-1);
block = strjoin(dataLines(:)', newline);
C = textscan(block, '%f%f%f%f%f%f%f%f%f%f%f%f', 'Delimiter', ',', 'CollectOutput', true);
data = C{1};

frame    = data(:,1);       % 帧号
subFrame = data(:,2);       % 帧内子样本号
emgRaw   = data(:,3:12);    % EMG7 ~ EMG16, 单位 V

nSamples = size(emgRaw, 1);
nCh      = numel(channels);
nSub     = max(subFrame) + 1;                    % 每帧子样本数 (本数据为 20)
t        = ((frame - 1) * nSub + subFrame) / fs; % 时间轴 (s)
fprintf('采样率: %d Hz | 样本数: %d | 时长: %.2f s\n', fs, nSamples, t(end));

%% 3. 去伪迹 (剔除异常尖峰 + PCHIP 插值修复) -----------------------------
if enableDespike
    fprintf('\n去伪迹: 判据 |x| > %.2f * P99.9(该通道), 标记区间外扩 ±%d 样本 (%.1f ms)\n', ...
            kP, despikeMargin, despikeMargin/fs*1000);
    [emgClean, replacedMask, dRows] = despikeEmg(emgRaw, t, fs, kP, despikeMargin, despikeNeighbors);

    if isempty(dRows)
        artTab = table();
        fprintf('未检出异常尖峰, 数据未被修改。\n');
    else
        artTab = cell2table(dRows, 'VariableNames', ...
            {'ChIdx','sIdx','eIdx','t_start_s','t_end_s','SamplesFlagged','Interval_ms', ...
             'PeakRemoved_V','PeakOverP999','PeakOverSigma','SupportPoints'});
        nEv = height(artTab);
        artTab.EMG    = reshape(channels(artTab.ChIdx),   nEv, 1);   % 显式列向量, 避免方向歧义
        artTab.Muscle = reshape(shortNames(artTab.ChIdx), nEv, 1);
        artTab.ChIdx  = [];
        artTab = movevars(artTab, {'EMG','Muscle'}, 'Before', 't_start_s');
        artTab = sortrows(artTab, 'PeakOverSigma', 'descend');

        fprintf('共替换 %d 个区间, %d 个样本 (等效占单通道记录长度的 %.3f%%)\n', ...
                height(artTab), nnz(replacedMask), ...
                100 * nnz(replacedMask) / nSamples);
        nShowTab = min(20, height(artTab));
        disp(artTab(1:nShowTab, :));
        if height(artTab) > nShowTab
            fprintf('... 其余 %d 个区间见 semg_artifact_report.csv\n', height(artTab)-nShowTab);
        end
    end
else
    emgClean     = emgRaw;
    replacedMask = false(nSamples, nCh);
    artTab       = table();
    fprintf('\n已按参数设置跳过去伪迹 (enableDespike = false)。\n');
end

%% 4. sEMG 处理: 高通 -> 整流 -> 低通 -> 最大值归一化 ----------------------
assert(hpFc < fs/2, '高通截止频率必须小于 Nyquist 频率 %.1f Hz', fs/2);
[hpB, hpA] = butter(filtOrder, hpFc/(fs/2), 'high');   % 10 Hz 二阶高通
[lpB, lpA] = butter(filtOrder, lpFc/(fs/2), 'low');    % 5 Hz 二阶低通

emgEnv     = computeEnvelope(emgClean, hpB, hpA, lpB, lpA, zeroPhase);  % 去伪迹后
emgEnvNorm = emgEnv ./ max(emgEnv, [], 1);                             % 4) 最大值归一化

% 对照: 未去伪迹的包络(仅用于验证去伪迹效果, 不参与输出数据)
emgEnvBef     = computeEnvelope(emgRaw, hpB, hpA, lpB, lpA, zeroPhase);
emgEnvNormBef = emgEnvBef ./ max(emgEnvBef, [], 1);

% 修复前后对比 (max/P99.9 是关键指标: 接近 1~2 说明最大值已由真实 EMG 决定)
p999Raw = zeros(1, nCh);  p999Cln = zeros(1, nCh);
for i = 1:nCh
    p999Raw(i) = pctlSorted(abs(emgRaw(:,i)),   0.999);
    p999Cln(i) = pctlSorted(abs(emgClean(:,i)), 0.999);
end
fprintf('\n================ 去伪迹前后对比 ================\n');
fprintf('%-5s %-8s %13s %13s %13s %13s %12s %12s\n', ...
    'EMG','muscle','max/P99.9前','max/P99.9后','峰值时刻前(s)','峰值时刻后(s)','env>0.1前','env>0.1后');
for i = 1:nCh
    [~, ipkB] = max(emgEnvBef(:,i));
    [~, ipkA] = max(emgEnv(:,i));
    fprintf('%-5d %-8s %13.1f %13.1f %13.2f %13.2f %12.3f %12.3f\n', ...
        channels(i), shortNames{i}, ...
        max(abs(emgRaw(:,i)))/p999Raw(i), max(abs(emgClean(:,i)))/p999Cln(i), ...
        t(ipkB), t(ipkA), mean(emgEnvNormBef(:,i) > 0.1), mean(emgEnvNorm(:,i) > 0.1));
end
fprintf('(env>0.1 = 归一化包络中大于 0.1 的时间占比; 修复前该值极小即"被压小")\n');

% 各肌肉峰值出现时刻(修复后): 归一化按全时程最大值, 另列出显示区间内的峰值位置
wIdxShow = find(t >= plotWin(1) & t <= plotWin(2));
fprintf('\n%-6s %-8s %-16s %-16s %-18s %-16s\n', ...
    'EMG', 'Muscle', '全时程峰值(s)', 'Peak (norm)', '显示区间内峰值(s)', 'Peak (norm)');
for i = 1:nCh
    [pv, ipk] = max(emgEnv(:,i));
    [~, iw]   = max(emgEnv(wIdxShow, i));
    fprintf('%-6d %-8s %-16.2f %-16.3f %-18.2f %-16.3f\n', ...
        channels(i), shortNames{i}, t(ipk), emgEnvNorm(ipk,i), ...
        t(wIdxShow(iw)), emgEnvNorm(wIdxShow(iw), i));
end

%% 5. 图 1: 每块肌肉一行, 左列 raw(去伪迹后), 右列处理后 ------------------
% 只显示 plotWin 区间 (s): 全时程 249.6 s 压缩到一张图里看不清细节
winIdx = wIdxShow;                     % 只显示 plotWin 区间
idx    = winIdx(1:plotDecimate:end);   % 仅绘图时降采样
fig  = figure('Name', 'sEMG raw vs processed', 'Color', 'w', ...
              'Position', [60 30 1500 1650]);
for i = 1:nCh
    % 左列: 原始 sEMG (去伪迹后, 红点 = 被替换的样本)
    subplot(nCh, 2, 2*i-1);
    plot(t(idx), emgClean(idx, i), 'Color', [0.15 0.35 0.65]); hold on;
    idxRep = find(replacedMask(:, i) & t >= plotWin(1) & t <= plotWin(2));
    plot(t(idxRep), emgClean(idxRep, i), 'r.', 'MarkerSize', 9);
    grid on; box off;
    ylabel({shortNames{i}, 'Raw (V)'});
    title(muscleNames{i}, 'FontWeight', 'normal');
    xlim(plotWin);

    % 右列: 处理后 sEMG (归一化包络)
    subplot(nCh, 2, 2*i);
    plot(t(idx), emgEnvNorm(idx, i), 'Color', [0.85 0.25 0.20], 'LineWidth', 0.8);
    grid on; box off;
    ylabel({shortNames{i}, 'Processed (norm.)'});
    title('Processed', 'FontWeight', 'normal');
    xlim(plotWin);
    ylim([0, 1.05]);
end
subplot(nCh, 2, 1);
title({sprintf('RAW sEMG (去伪迹后, 红点=被替换样本)  显示 %.0f-%.0f s', plotWin(1), plotWin(2)); ...
       muscleNames{1}}, 'FontWeight', 'bold');
subplot(nCh, 2, 2);
title({'PROCESSED sEMG (归一化仍按全时程最大值)'; ...
       '10 Hz HP | rectify | 5 Hz LP | max-norm'}, 'FontWeight', 'bold');
subplot(nCh, 2, 2*nCh-1); xlabel('Time (s)');
subplot(nCh, 2, 2*nCh);   xlabel('Time (s)');
sgtitle(sprintf('sEMG Raw and Processed (Level Walking, fs = %d Hz, 显示区间 %.0f-%.0f s / 全长 %.1f s)', ...
        fs, plotWin(1), plotWin(2), t(end)), 'FontSize', 13);

%% 6. 图 2: 去伪迹前后归一化包络对比 (同样只显示 plotWin 区间) -----------
fig2 = figure('Name', 'despike effect', 'Color', 'w', 'Position', [80 30 1400 1650]);
for i = 1:nCh
    subplot(nCh, 1, i);
    plot(t(idx), emgEnvNormBef(idx, i), 'Color', [0.65 0.65 0.65]); hold on;
    plot(t(idx), emgEnvNorm(idx, i), 'Color', [0.10 0.30 0.75], 'LineWidth', 0.8);
    grid on; box off; xlim(plotWin); ylim([0, 1.05]);
    ylabel(shortNames{i});
    if i == 1
        title('去伪迹前(灰) vs 去伪迹后(蓝) 的归一化包络');
        legend({'despike 前', 'despike 后'}, 'Location', 'northeast');
    end
end
xlabel('Time (s)');
sgtitle(sprintf('去伪迹对归一化包络的影响 (灰: 未去伪迹, 最大值被伪迹绑架 → 整体被压小; 显示 %.0f-%.0f s)', ...
        plotWin(1), plotWin(2)), 'FontSize', 13);

%% 7. 保存结果 -----------------------------------------------------------
figName = fullfile(outDir, 'level_walking_semg_raw_and_processed');
saveas(fig,  [figName '.png']);
savefig(fig, figName);
saveas(fig2, fullfile(outDir, 'level_walking_semg_despike_effect.png'));

T = array2table([t(:), emgEnvNorm], 'VariableNames', ...
    [{'Time_s'}, arrayfun(@(c) sprintf('EMG%d', c), channels, 'UniformOutput', false)]);
writetable(T, fullfile(outDir, 'level_walking_semg_processed.csv'));

if ~isempty(artTab)
    writetable(artTab, fullfile(outDir, 'semg_artifact_report.csv'));
end

% 只保存必要数据(验证用的 emgEnvBef/emgEnvNormBef 属中间结果, 不入库;
%  需要时把 enableDespike 置 false 重跑即可得到)
save(fullfile(outDir, 'level_walking_semg_processed.mat'), ...
     't', 'fs', 'emgRaw', 'emgClean', 'replacedMask', ...
     'emgEnv', 'emgEnvNorm', 'artTab', ...
     'channels', 'muscleNames', 'shortNames', 'hpFc', 'lpFc', 'filtOrder', ...
     'zeroPhase', 'kP', 'plotWin', '-v7.3');

fprintf('\n已保存输出到: %s\n', outDir);

%% 局部函数 -------------------------------------------------------------
function env = computeEnvelope(x, hpB, hpA, lpB, lpA, zeroPhase)
% 高通 -> 全波整流 -> 低通, 得到线性包络
n   = size(x, 1);
env = zeros(size(x));
for i = 1:size(x, 2)
    if zeroPhase
        y = filtfilt(hpB, hpA, x(:, i));
    else
        y = filter(hpB, hpA, x(:, i));
    end
    y = abs(y);
    if zeroPhase
        y = filtfilt(lpB, lpA, y);
    else
        y = filter(lpB, lpA, y);
    end
    env(:, i) = max(y, 0);      % 零相位滤波边缘可能有极小负振铃
end
end

function [xClean, maskReplaced, rows] = despikeEmg(x, t, fs, kP, marginSamples, nbrSamples)
% 去伪迹: 检出异常尖峰(|x| > kP*P99.9), 用 PCHIP 从两侧干净样本插值重建
[n, nCh]      = size(x);
xClean        = x;
maskReplaced  = false(n, nCh);
rows          = {};
for i = 1:nCh
    xi   = x(:, i);
    ax   = abs(xi);
    p999 = pctlSorted(ax, 0.999);
    mask = ax > kP * p999;                      % 主判据
    if ~any(mask), continue; end
    sigma = 1.4826 * pctlSorted(abs(xi - pctlSorted(xi, 0.5)), 0.5);
    % 标记区间向两侧外扩, 覆盖脉冲邻域; 相邻标记合并为同一区间
    maskD = conv(double(mask), ones(2*marginSamples+1, 1), 'same') > 0;
    maskReplaced(:, i) = maskD;
    [s0, e0] = getIntervals(maskD);
    for k = 1:numel(s0)
        s = s0(k);
        e = e0(k);
        nFlagged = nnz(mask(s:e));
        peakVal  = max(ax(s:e));
        % 支撑点: 区间两侧各 nbrSamples 个"未被替换"的原始样本
        lo  = max(1, s - nbrSamples);
        hi  = min(n, e + nbrSamples);
        sup = (lo:hi)';
        sup = sup(~maskD(sup));
        if numel(sup) < 4
            sup = find(~maskD);                 % 退化情形: 用全记录干净点
        end
        if numel(sup) < 2
            xi(s:e) = 0;                        % 极端退化: 置零
        else
            xi(s:e) = interp1(sup, xi(sup), (s:e)', 'pchip');
        end
        rows(end+1, :) = {i, s, e, t(s), t(e), nFlagged, (e-s+1)/fs*1000, ...
                          peakVal, peakVal/p999, peakVal/sigma, numel(sup)}; %#ok<SAGROW>
    end
    xClean(:, i) = xi;
end
end

function v = pctlSorted(x, p)
% 排序法求分位数(不依赖 Statistics Toolbox)
xs = sort(x(:));
n  = numel(xs);
if n == 0, v = NaN; return; end
k  = min(n, max(1, ceil(p * n)));
v  = xs(k);
end

function [s0, e0] = getIntervals(mask)
% 把逻辑掩码转换成若干连续区间的起止下标
mask = mask(:) > 0;
d  = diff([0; mask; 0]);
s0 = find(d == 1);
e0 = find(d == -1) - 1;
end
