%% diagnose_semg_artifacts.m
% =========================================================================
% sEMG 异常尖峰(伪迹)诊断程序 —— 只做统计与画图, 不滤波 / 不归一化 / 不修改任何数据
%
% 诊断内容:
%   1) 每通道稳健统计: max|x|、P99.9、MAD 稳健标准差 sigma_robust、max/P99.9、max/sigma
%      max/P99.9 是"伪迹主导程度"指标: 接近 1~2 = 正常; > 5 = 最大值被伪迹绑架
%   2) 局部稳健基线 baseline(n) = P75( |x| , n 前后各 T/2 ), 相对幅度比
%         ratio(n) = |x(n)| / baseline(n)                       <-- 判据 A 的输入
%      用两种窗口 T = 200 ms 与 600 ms 各算一遍并对比:
%        T 太短时, 相位性肌肉(如腓肠肌)在爆发边缘会被误判; T 长一些可明显减少误检
%   3) ratio 分布直方图(对数横轴, 两种窗口叠加):
%      真实 EMG 聚在左侧, 伪迹孤立在右侧, 中间的"断口"就是阈值取值位置
%   4) 阈值扫描: 不同 k1(相对倍数)、k2(绝对倍数) 下被判为可疑的样本数
%   5) 候选伪迹事件清单(时间/真实长度/峰值/相对倍数/绝对倍数/跨通道共现)
%   6) 每通道"可疑样本 Top-N"(按相对倍数 与 按绝对幅度 各取 N 个, 取并集),
%      避免漏掉"较宽但相对倍数不高"的瞬变
%   7) 图: ratio 直方图 / 原始信号标红 / 可疑事件的 ±0.5 s 放大窗
%
% 判据说明(本数据集实测结论):
%   判据B(绝对, 推荐主判据) |x| > kP * P99.9(该通道本身), 推荐 kP ≈ 2.5
%       P99.9 已包含该通道真实爆发的幅度水平, 因此对强相位性肌肉天然适应;
%       实测: 伪迹峰值 = 该通道 P99.9 的 4~30 倍, 而干净通道(RF_R)的最大值仅 1.4 倍,
%       且 RF_R 在整个记录中没有任何样本超过 2.5*P99.9 —— 说明该判据不会误伤真实爆发。
%   判据A(相对) ratio = |x| / P75(|x|, ±100 ms) > k1
%       单用不够: 伪迹若发生在真实爆发期内, 局部 P75 也高, 相对倍数会降到 15 以下,
%       于是漏检 (如 RF_L 52.44 s、BF_R 184.80 s 的伪迹);
%       在强相位性通道(如 GAS_R)上又会把真实爆发峰误标。故它只作为辅助说明量。
%       (注意: 也不要用全局 sigma_robust 作基准 —— 强相位性肌肉的真实爆发峰可达 45~48 sigma,
%        用 10*sigma 会把真实爆发峰误判为伪迹, 程序里保留了该对比表以便查看。)
%   最终判据 = 判据B (可选用 requireRelative 加严), 再配合"图3放大窗"人工确认。
%
% 说明:
%   - 全部计算都在"原始域"进行(去伪迹必须在滤波之前, 否则脉冲能量已被扩散)
%   - sigma_robust = 1.4826 * median(|x - median(x)|), 对少量尖峰稳健
%   - 分位数用排序法自行实现, 不依赖 Statistics Toolbox
%
% 输出 (experiments\moco\outputs\semg\diagnostics\):
%   diag_01_ratio_histogram.png   ratio 分布直方图 (10 通道, 两种基线窗口叠加)
%   diag_02_raw_marked.png        原始信号 + 候选伪迹样本标红
%   diag_03_event_zoom.png        可疑事件放大窗 (覆盖所有通道, 每通道取最大的几个)
%   diag_artifact_candidates.csv  候选伪迹事件清单
%   diag_top_samples.csv          每通道可疑样本 Top-N
%   diag_channel_summary.csv      每通道稳健统计汇总
%   diag_semg_data.mat            ratio / baseline / 统计量 (供后续定阈值用)
% =========================================================================

clear; close all; clc;

%% 1. 参数 ---------------------------------------------------------------
scriptDir = fileparts(mfilename('fullpath'));   % experiments\moco
semgFile  = fullfile(scriptDir, '..', 'data', 'SQR_walking', 'level_walking_semg.csv');
outDir    = fullfile(scriptDir, 'outputs', 'semg', 'diagnostics');
if ~exist(outDir, 'dir'); mkdir(outDir); end

baselineWinList = [0.20 0.60];  % 局部基线窗口长度 (s), 各算一遍做对比
baselineStep    = 0.05;         % 基线滑动步长 (s)
pctl            = 0.75;         % 基线所用的分位数 (P75)
wRel            = 1;            % 相对判据(判据A)所用基线窗口序号 (1 = 200 ms)
growSamples     = 10;           % 候选事件向两侧扩展的样本数 (5 ms @2000 Hz)
kPList          = [2 2.5 3 4];  % 绝对幅度判据(判据B): |x| > kP * P99.9(该通道本身)
kPUse           = 2.5;          % 提取候选事件所用 kP
requireRelative = false;        % true: 还要求同时满足相对判据 ratio > k1Use (更保守, 但会漏掉
                                %       "发生在真实爆发期内"的伪迹, 见文末判读要点)
k1List          = [8 10 15 20 30 50 100];   % 相对倍数阈值扫描
k1Use           = 15;           % 提取候选事件所用相对阈值 (判据A)
k2List          = [8 10];       % 旧判据(×sigma_robust)阈值扫描, 仅作对比说明
topN            = 10;           % 每通道按"相对倍数""绝对幅度"各取前 N 个可疑样本
maxPerChannel   = 4;            % 放大窗图中每通道最多展示的事件数
maxZoomEvents   = 40;           % 放大窗图总事件数上限
zoomHalfWin     = 0.5;          % 放大窗半宽 (s)
plotDecimate    = 5;            % 绘图降采样(仅影响显示速度)
dominatedRatio  = 5;            % max/P99.9 超过该值视为"最大值被伪迹主导"

% 通道与肌肉对应关系
channels = 7:16;
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
shortNames = {'RF_R','VL_R','BF_R','TA_R','GAS_R','RF_L','VL_L','BF_L','TA_L','GAS_L'};

%% 2. 读取数据 (只取 CSV 第一段 EMG 数据块) -----------------------------
fprintf('读取数据: %s\n', semgFile);
raw   = fileread(semgFile);
lines = regexp(raw, '\r?\n', 'split')';
lines = lines(:);
fs = str2double(strtrim(lines{2}));             % 采样率 (2000 Hz)
blankIdx  = find(cellfun(@(s) isempty(strtrim(s)), lines(6:end)), 1, 'first') + 5;
dataLines = lines(6:blankIdx-1);
C = textscan(strjoin(dataLines(:)', newline), '%f%f%f%f%f%f%f%f%f%f%f%f', ...
             'Delimiter', ',', 'CollectOutput', true);
data = C{1};
frame    = data(:,1);
subFrame = data(:,2);
emgRaw   = data(:,3:12);                        % EMG7~EMG16 (V)

nSamples = size(emgRaw, 1);
nCh      = numel(channels);
nSub     = max(subFrame) + 1;
t        = ((frame - 1) * nSub + subFrame) / fs;
fprintf('采样率: %d Hz | 样本数: %d | 时长: %.2f s\n', fs, nSamples, t(end));

%% 3. 稳健统计 + 局部基线 + 相对幅度比 -----------------------------------
axAll    = abs(emgRaw);
nWin     = numel(baselineWinList);
baseAll  = cell(1, nWin);       % 每种窗口长度的基线
ratioAll = cell(1, nWin);       % 每种窗口长度的 ratio
sigmaRob = zeros(1, nCh);
maxVal   = zeros(1, nCh);
tMax     = zeros(1, nCh);
p999     = zeros(1, nCh);
nAb90    = zeros(1, nCh);       % 超过 90% 峰值的样本数(判断单点还是平台)

for w = 1:nWin
    W    = round(baselineWinList(w) * fs);
    step = max(1, round(baselineStep * fs));
    starts  = (1:step:(nSamples - W + 1))';
    centers = starts + (W - 1) / 2;
    idxMat  = starts + (0:W-1);                 % nB × W 索引矩阵
    pPos    = max(1, round(pctl * W));          % P75 在窗内的排序位置
    baseAll{w}  = zeros(nSamples, nCh);
    ratioAll{w} = zeros(nSamples, nCh);
    fprintf('计算局部 P75 基线: 窗口 %.0f ms, 步长 %.0f ms (%d 个窗口) ...\n', ...
            baselineWinList(w)*1000, baselineStep*1000, size(idxMat,1));
    for i = 1:nCh
        % 注意: A(idxMat, i) 会被线性化, 必须 reshape 回 nB × W
        vals = reshape(axAll(idxMat(:), i), size(idxMat));
        lp   = sort(vals, 2);
        lp   = lp(:, pPos);
        base = max(interp1(centers, lp, (1:nSamples)', 'linear', 'extrap'), realmin);
        baseAll{w}(:, i)  = base;
        ratioAll{w}(:, i) = axAll(:, i) ./ base;
    end
end

% 稳健统计(与基线窗口无关)
for i = 1:nCh
    x = emgRaw(:, i);
    sigmaRob(i) = 1.4826 * pctlSorted(abs(x - pctlSorted(x, 0.5)), 0.5);
    [maxVal(i), iM] = max(axAll(:, i));
    tMax(i)  = t(iM);
    p999(i)  = pctlSorted(axAll(:, i), 0.999);
    nAb90(i) = nnz(axAll(:, i) > 0.9 * maxVal(i));
end

%% 4. 报告 1: 每通道稳健统计 --------------------------------------------
fprintf('\n================ 1) 每通道稳健统计 ================\n');
fprintf('%-5s %-8s %12s %10s %12s %11s %13s %10s %8s %-22s\n', ...
    'EMG','muscle','max|x|(V)','t_max(s)','P99.9(V)','max/P99.9','sigma_rob(V)','max/sigma','n>0.9max','verdict');
for i = 1:nCh
    if maxVal(i)/p999(i) > dominatedRatio
        verdict = '最大值被伪迹主导';
    else
        verdict = '未见孤立伪迹主导';
    end
    fprintf('%-5d %-8s %12.3e %10.2f %12.3e %11.1f %13.3e %10.1f %8d %-22s\n', ...
        channels(i), shortNames{i}, maxVal(i), tMax(i), p999(i), ...
        maxVal(i)/p999(i), sigmaRob(i), maxVal(i)/sigmaRob(i), nAb90(i), verdict);
end

%% 5. 报告 2: 绝对幅度判据 |x| > kP * P99.9(该通道) ----------------------
% P99.9 已经包含了该通道"真实爆发"的幅度水平, 因此以它为基准的判据
% 天然适应强相位性肌肉, 不会把真实爆发峰当成伪迹(这是本数据集的推荐主判据)
fprintf('\n================ 2) 绝对幅度判据: |x| > kP * P99.9(该通道) ================\n');
fprintf('%-5s %-8s %12s %11s', 'EMG','muscle','P99.9(V)','max/P99.9');
for k = kPList, fprintf('%14s', sprintf('n(>%gxP99.9)', k)); end
fprintf('\n');
for i = 1:nCh
    fprintf('%-5d %-8s %12.3e %11.1f', channels(i), shortNames{i}, p999(i), maxVal(i)/p999(i));
    for k = kPList
        fprintf('%14d', nnz(axAll(:,i) > k * p999(i)));
    end
    fprintf('\n');
end
fprintf('判读: 干净通道(如 RF_R)在各 kP 下计数均为 0; 有伪迹的通道只在 kP 较小时出现个位数计数。\n');

%% 6. 报告 3: ratio 分布与阈值扫描 ---------------------------------------
bands = [1 8; 8 15; 15 30; 30 100; 100 inf];
for w = 1:nWin
    fprintf('\n================ 2) ratio 分布 (基线窗口 %.0f ms) ================\n', baselineWinList(w)*1000);
    fprintf('(ratio = |x| / P75(|x|, ±%.0f ms))\n', baselineWinList(w)*1000/2);
    fprintf('%-5s %-8s %10s %10s %10s %10s %10s\n', ...
        'EMG','muscle','[1,8)','[8,15)','[15,30)','[30,100)','[100,inf)');
    for i = 1:nCh
        r   = ratioAll{w}(:, i);
        cnt = zeros(1, size(bands,1));
        for k = 1:size(bands,1)
            cnt(k) = nnz(r >= bands(k,1) & r < bands(k,2));
        end
        fprintf('%-5d %-8s %10d %10d %10d %10d %10d\n', channels(i), shortNames{i}, cnt);
    end
end

% 两种窗口的"分布干净程度"对比: 右侧区间计数越少, 断口越清晰
fprintf('\n窗口长度对比 (全部通道合计):\n');
fprintf('%-12s %10s %10s %10s\n', 'baseline', 'n[8,15)', 'n[15,30)', 'n[30,inf)');
for w = 1:nWin
    rr = ratioAll{w};
    fprintf('%-12s %10d %10d %10d\n', sprintf('%.0f ms', baselineWinList(w)*1000), ...
        nnz(rr >= 8 & rr < 15), nnz(rr >= 15 & rr < 30), nnz(rr >= 30));
end

% 候选事件用"较短窗口"(断口更清晰) + "P99.9 绝对判据"
wUse = wRel;
rUse = ratioAll{wUse};
fprintf('\n================ 3) 阈值扫描 (基线窗口 %.0f ms) ================\n', baselineWinList(wUse)*1000);
fprintf('\n[A] 仅判据A: ratio > k1\n');
fprintf('%-5s %-8s', 'EMG','muscle');
for k = k1List, fprintf('%10s', sprintf('>%g', k)); end
fprintf('\n');
for i = 1:nCh
    fprintf('%-5d %-8s', channels(i), shortNames{i});
    for k = k1List, fprintf('%10d', nnz(rUse(:,i) > k)); end
    fprintf('\n');
end
fprintf('\n[B-推荐] 判据B 单独使用: |x| > kP*P99.9 (见表 2, 这是本数据集的推荐主判据)\n');
fprintf('\n[A∧B] 判据A ∧ 判据B: ratio > k1 且 |x| > %g*P99.9 (过于保守: 会漏掉爆发期内的伪迹)\n', kPUse);
fprintf('%-5s %-8s', 'EMG','muscle');
for k = k1List, fprintf('%10s', sprintf('>%g', k)); end
fprintf('\n');
for i = 1:nCh
    big = axAll(:,i) > kPUse * p999(i);
    fprintf('%-5d %-8s', channels(i), shortNames{i});
    for k = k1List, fprintf('%10d', nnz(rUse(:,i) > k & big)); end
    fprintf('\n');
end
for k2 = k2List
    fprintf('\n[B-对比] 判据A ∧ |x| > %g*sigma_robust (sigma 是全段稳健标准差, 会误检强相位性肌肉的真实爆发峰)\n', k2);
    fprintf('%-5s %-8s', 'EMG','muscle');
    for k = k1List, fprintf('%10s', sprintf('>%g', k)); end
    fprintf('\n');
    for i = 1:nCh
        big = axAll(:,i) > k2 * sigmaRob(i);
        fprintf('%-5d %-8s', channels(i), shortNames{i});
        for k = k1List, fprintf('%10d', nnz(rUse(:,i) > k & big)); end
        fprintf('\n');
    end
end

%% 7. 候选伪迹事件 (判据B, 可选 A∧B) + 真实长度 + 精确 P75 复算 ----------
mask = false(nSamples, nCh);
for i = 1:nCh
    mask(:, i) = axAll(:, i) > kPUse * p999(i);
    if requireRelative
        mask(:, i) = mask(:, i) & rUse(:, i) > k1Use;
    end
end
maskDil = false(nSamples, nCh);
for i = 1:nCh
    maskDil(:, i) = conv(double(mask(:, i)), ones(2*growSamples+1, 1), 'same') > 0;
end

Wuse = round(baselineWinList(wUse) * fs);
evCh = []; evMus = {}; evT0 = []; evT1 = []; evNs = []; evInt = []; evPk = []; evPkT = [];
evP75 = []; evR = []; evRs = []; evOp = [];
for i = 1:nCh
    [s0, e0] = getIntervals(maskDil(:, i));
    for k = 1:numel(s0)
        seg = s0(k):e0(k);
        [pv, ip] = max(axAll(seg, i));
        pk = seg(ip);
        % 该点周围 ±T/2 内"未被标红"的干净样本, 精确复算 P75
        w0 = max(1, pk - round(Wuse/2));
        w1 = min(nSamples, pk + round(Wuse/2));
        wn = (w0:w1)';
        clean = wn(~maskDil(wn, i));
        if numel(clean) < 50, clean = wn; end
        evCh(end+1,1)  = channels(i);
        evMus{end+1,1} = shortNames{i};
        evT0(end+1,1)  = t(seg(1));
        evT1(end+1,1)  = t(seg(end));
        evNs(end+1,1)  = nnz(mask(seg, i));             % 真实被判定的样本数
        evInt(end+1,1) = numel(seg) / fs * 1000;        % 含扩展后的区间长度 (ms)
        evPk(end+1,1)  = pv;
        evPkT(end+1,1) = t(pk);
        evP75(end+1,1) = pctlSorted(axAll(clean, i), pctl);
        evR(end+1,1)   = pv / evP75(end);
        evRs(end+1,1)  = pv / sigmaRob(i);
        evOp(end+1,1)  = pv / p999(i);                  % 峰值 / 该通道 P99.9
    end
end

% 跨通道共现(时间重叠, 容差 ±10 ms)
nEv = numel(evCh);
coinStr = repmat("-", nEv, 1);
for a = 1:nEv
    others = [];
    for b = 1:nEv
        if b == a || evCh(b) == evCh(a), continue; end
        if evT0(b) <= evT1(a) + 0.01 && evT1(b) >= evT0(a) - 0.01
            others(end+1) = evCh(b); %#ok<SAGROW>
        end
    end
    if ~isempty(others)
        coinStr(a) = strjoin(string(unique(others)), ",");
    end
end

[~, ord] = sort(evRs, 'descend');       % 按"绝对倍数"排序, 最极端的排在最前
evTab = table(evCh(ord), evMus(ord), (1:nEv)', evT0(ord), evT1(ord), evNs(ord), evInt(ord), ...
              evPk(ord), evPkT(ord), evP75(ord), evR(ord), evOp(ord), evRs(ord), coinStr(ord), ...
    'VariableNames', {'EMG','Muscle','EventIdx','t_start_s','t_end_s','ArtifactSamples','Interval_ms', ...
                      'Peak_V','PeakTime_s','LocalP75_V','Ratio_local','PeakOverP999','Ratio_sigmaRobust', ...
                      'CoincidentChannels'});

if requireRelative
    critStr = sprintf(' 且 ratio>%g', k1Use);
else
    critStr = '';
end
fprintf('\n================ 4) 候选伪迹事件 (|x|>%g*P99.9%s, 按绝对倍数排序) ================\n', kPUse, critStr);
if isempty(evTab)
    fprintf('未检出候选伪迹事件。\n');
else
    fprintf('共 %d 个事件 (ArtifactSamples = 真实被判定的样本数, Interval_ms = 含扩展的区间):\n', height(evTab));
    nShowTab = min(40, height(evTab));
    disp(evTab(1:nShowTab, :));
    if height(evTab) > nShowTab
        fprintf('... 其余 %d 个事件见 CSV 文件。\n', height(evTab)-nShowTab);
    end
end

%% 7. 每通道可疑样本 Top-N (按相对倍数 与 按绝对幅度 取并集) --------------
suspect = false(nSamples, nCh);
topRows = {};
for i = 1:nCh
    [~, o1] = maxk(rUse(:, i), topN);           % 相对倍数最大的 N 个
    [~, o2] = maxk(axAll(:, i), topN);          % 绝对幅度最大的 N 个
    o = unique([o1; o2]);
    suspect(o, i) = true;
    for k = 1:numel(o)
        topRows(end+1, :) = {channels(i), shortNames{i}, t(o(k)), axAll(o(k), i), ...
            baseAll{wUse}(o(k), i), rUse(o(k), i), axAll(o(k), i)/p999(i), ...
            axAll(o(k), i)/sigmaRob(i), logical(maskDil(o(k), i))}; %#ok<SAGROW>
    end
end
topTab = cell2table(topRows, 'VariableNames', {'EMG','Muscle','Time_s','AbsAmplitude_V', ...
    'LocalP75_V','Ratio_local','PeakOverP999','Ratio_sigmaRobust','FlaggedByThresholds'});
topTab = sortrows(topTab, 'Ratio_local', 'descend');

% 取"每个独立事件"的代表点用于放大窗: 按区间聚类, 每通道取幅度最大的几个
zoomRows = {};
for i = 1:nCh
    [s0, e0] = getIntervals(conv(double(suspect(:,i)), ones(2*growSamples+1,1), 'same') > 0);
    pkInfo = zeros(numel(s0), 3);   % [峰值幅度, 峰值索引, 相对倍数]
    for k = 1:numel(s0)
        seg = s0(k):e0(k);
        [pv, ip] = max(axAll(seg, i));
        pkInfo(k, :) = [pv, seg(ip), rUse(seg(ip), i)];
    end
    if isempty(pkInfo), continue; end
    [~, oo] = sort(pkInfo(:,1), 'descend');
    oo = oo(1:min(maxPerChannel, numel(oo)));
    for k = 1:numel(oo)
        zoomRows(end+1, :) = {channels(i), shortNames{i}, pkInfo(oo(k),2), ...
            pkInfo(oo(k),1), pkInfo(oo(k),3)}; %#ok<SAGROW>
    end
end
zoomTab = cell2table(zoomRows, 'VariableNames', {'EMG','Muscle','PeakIdx','Peak_V','Ratio_local'});
zoomTab = sortrows(zoomTab, 'Ratio_local', 'descend');
if height(zoomTab) > maxZoomEvents
    zoomTab = zoomTab(1:maxZoomEvents, :);
end

%% 8. 图 1: ratio 分布直方图 (两种基线窗口叠加) --------------------------
edges = -0.5:0.05:3.5;
ctr   = edges(1:end-1) + 0.025;
f1 = figure('Name','ratio histogram','Color','w','Position',[80 40 1400 1250]);
for i = 1:nCh
    subplot(5, 2, i);
    cnt1 = histcounts(log10(max(ratioAll{1}(:,i), 1e-3)), edges);
    cnt2 = histcounts(log10(max(ratioAll{nWin}(:,i), 1e-3)), edges);
    bar(ctr, max(cnt2, 0.5), 1, 'FaceColor', [0.85 0.45 0.30], 'EdgeColor', 'none'); hold on;
    bar(ctr, max(cnt1, 0.5), 1, 'FaceColor', [0.30 0.50 0.80], 'EdgeColor', 'none');
    set(gca, 'YScale', 'log'); grid on; box off;
    yl = ylim;
    plot([log10(k1Use) log10(k1Use)], yl, 'k--', 'LineWidth', 1.2);
    xlim([-0.5 3.5]);
    title(sprintf('%s: max/P99.9 = %.1f, n(ratio>%g) = %d / %d', shortNames{i}, ...
        maxVal(i)/p999(i), k1Use, ...
        nnz(ratioAll{1}(:,i) > k1Use), nnz(ratioAll{nWin}(:,i) > k1Use)), 'FontWeight', 'normal');
    if i >= nCh-1, xlabel('log_{10}( |x| / P75 local )'); end
    if mod(i,2) == 1, ylabel('count'); end
end
sgtitle(sprintf('相对幅度比分布  蓝: 基线窗口 %.0f ms | 橙: %.0f ms   (黑虚线 = 阈值 %g)', ...
    baselineWinList(1)*1000, baselineWinList(nWin)*1000, k1Use), 'FontSize', 13);

%% 9. 图 2: 原始信号 + 候选样本标红 --------------------------------------
di = 1:plotDecimate:nSamples;
f2 = figure('Name','raw with candidates','Color','w','Position',[60 30 1500 1650]);
for i = 1:nCh
    subplot(nCh, 1, i);
    plot(t(di), emgRaw(di, i), 'Color', [0.20 0.35 0.60]); hold on;
    idxMark = find(maskDil(:, i));
    plot(t(idxMark), emgRaw(idxMark, i), 'r.', 'MarkerSize', 9);
    grid on; box off; xlim([0, t(end)]);
    ylabel(shortNames{i});
    if i == 1
        title(sprintf('原始 sEMG (蓝) 与候选伪迹样本 (红点), 判据: ratio>%g 且 |x|>%g*P99.9', k1Use, kPUse));
    end
end
xlabel('Time (s)');

%% 10. 图 3: 可疑事件放大窗 ---------------------------------------------
nShow = height(zoomTab);
if nShow > 0
    nCol = 4;
    nRow = ceil(nShow / nCol);
    f3 = figure('Name','event zoom','Color','w','Position',[40 20 1700 min(1600, 260*nRow+100)]);
    for k = 1:nShow
        chIdx = find(channels == zoomTab.EMG(k));
        pk = zoomTab.PeakIdx(k);
        w0 = max(1, pk - round(zoomHalfWin*fs));
        w1 = min(nSamples, pk + round(zoomHalfWin*fs));
        subplot(nRow, nCol, k);
        plot(t(w0:w1), emgRaw(w0:w1, chIdx), 'Color', [0.20 0.35 0.60]); hold on;
        plot(t(pk), emgRaw(pk, chIdx), 'ro', 'MarkerSize', 6, 'LineWidth', 1.2);
        grid on; box off;
        title(sprintf('%s  t=%.2f s\nratio=%.0f, |x|/sigma=%.0f', ...
            zoomTab.Muscle{k}, t(pk), zoomTab.Ratio_local(k), ...
            zoomTab.Peak_V(k)/sigmaRob(chIdx)), 'FontWeight', 'normal', 'FontSize', 8);
    end
    sgtitle(sprintf('可疑事件放大窗 (±%.1f s, 每通道最多 %d 个, 按相对倍数排序)', ...
        zoomHalfWin, maxPerChannel), 'FontSize', 12);
end

%% 11. 保存结果 ----------------------------------------------------------
sumTab = table(channels(:), shortNames(:), muscleNames(:), maxVal(:), tMax(:), p999(:), ...
               (maxVal./p999)', sigmaRob(:), (maxVal./sigmaRob)', nAb90(:), ...
    'VariableNames', {'EMG','Muscle','MuscleFull','MaxAbs_V','MaxTime_s','P999_V', ...
                      'MaxOverP999','SigmaRobust_V','MaxOverSigma','SamplesAbove90pctMax'});
writetable(sumTab, fullfile(outDir, 'diag_channel_summary.csv'));
writetable(evTab,   fullfile(outDir, 'diag_artifact_candidates.csv'));
writetable(topTab,  fullfile(outDir, 'diag_top_samples.csv'));
saveas(f1, fullfile(outDir, 'diag_01_ratio_histogram.png'));
saveas(f2, fullfile(outDir, 'diag_02_raw_marked.png'));
if exist('f3', 'var'), saveas(f3, fullfile(outDir, 'diag_03_event_zoom.png')); end
save(fullfile(outDir, 'diag_semg_data.mat'), 't', 'fs', 'channels', 'shortNames', ...
     'baselineWinList', 'sigmaRob', 'maxVal', 'tMax', 'p999', 'maskDil', 'suspect', ...
     'evTab', 'topTab', 'ratioAll', 'baseAll', '-v7.3');

fprintf('\n诊断结果已保存到: %s\n', outDir);
fprintf('判读要点:\n');
fprintf('  * max/P99.9 >> 1 的通道, 最大値由孤立伪迹决定, 归一化必然被压小 (本数据: 7/10 个通道);\n');
fprintf('  * 推荐主判据 |x| > %.1f*P99.9 (表 2), 它不会误伤真实爆发; 相对判据单独用会漏检;\n', kPUse);
fprintf('  * 图 3 的放大窗用于人工确认: 单点/短尖峰 = 伪迹; 有明确爆发形态 = 真实 EMG, 不应剔除。\n');

%% 局部函数 -------------------------------------------------------------
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
