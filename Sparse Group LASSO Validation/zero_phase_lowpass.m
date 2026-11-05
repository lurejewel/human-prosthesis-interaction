function y = zero_phase_lowpass(x, fs, fc, order)
%ZERO_PHASE_LOWPASS 对列信号做零相位 Butterworth 低通。
%
% 为什么必须零相位
% ----------------
% LASSO 的输入特征和 target 都带有明确的时间先后关系（反射延迟只有
% 10 ms 左右）。普通因果低通会引入与截止频率相关的相位延迟，相当于
% 人为改变反射延迟；filtfilt() 先正向、再反向滤波，净相位为 0，
% 只改变幅频特性，不移动峰值时刻。
%
% 为什么放在这里而不是直接用 filtfilt
% -----------------------------------
% * 需要处理多列信号（18 块肌肉 / 26 个特征）并保持列独立；
% * 需要先按列插值 NaN，否则 filtfilt 会报错或把 NaN 扩散到整列；
% * 需要处理“整列全 NaN”的可选通道（例如只有 5 块肌肉的 sEMG
%   excitation 文件），不能因为一个可选通道缺失就让整个接口失败；
% * 需要在没有 Signal Processing Toolbox 的机器上仍能工作。
%
% 输入
% ----
%   x     : nT×nC 或 nT×1 数值矩阵；
%   fs    : 采样频率 [Hz]；
%   fc    : 截止频率 [Hz]；fc<=0 或 fc>=fs/2 时原样返回；
%   order : Butterworth 阶数，默认 4。
%
% 输出
% ----
%   y     : 与 x 同尺寸（若 x 是行向量，输出为列向量，便于后续矩阵运算）。
%
% NaN 处理约定
% ------------
% * 部分 NaN：按列线性插值（含端点外推）后再滤波，输出不含 NaN；
% * 整列 NaN：保持 NaN，不参与滤波；调用者自己决定这是错误（activation）
%   还是可接受（可选 excitation 通道）。

% -------------------------------------------------------------------------
% 0. 参数校验
% -------------------------------------------------------------------------
if nargin < 3 || isempty(fs) || isempty(fc)
    error('zero_phase_lowpass:badInput', 'x, fs and fc are required.');
end
if nargin < 4 || isempty(order)
    order = 4;
end

y = x;
% 这些情况下不做任何滤波，直接返回原输入：
% * 空矩阵；
% * fc<=0：用户显式关闭滤波；
% * fc/fs 非有限或 fs<=0：参数无效，交给调用者检查；
% * fc>=Nyquist：没有可滤的频率。
if isempty(x) || fc <= 0 || ~isfinite(fc) || ~isfinite(fs) || fs <= 0
    return;
end
if fc >= fs / 2
    return;
end

% 记录输入是否是向量：向量输入时先转成 nT×1，输出再转回列向量。
wasVector = isvector(x);
if wasVector
    x = x(:);
end

% -------------------------------------------------------------------------
% 1. 逐列插值 NaN
% -------------------------------------------------------------------------
% filtfilt 要求输入有限。对每个通道单独插值，避免不同通道之间互相影响。
% 全 NaN 列跳过；后面 allBadCols 会把它们排除在滤波之外。
nT = size(x, 1);
for c = 1:size(x, 2)
    v = x(:, c);
    bad = ~isfinite(v);
    if all(bad)
        continue;
    end
    if any(bad)
        good = find(~bad);
        % 'extrap' 处理首尾 NaN；线性插值对短缺口足够，且不会引入振铃。
        v(bad) = interp1(good, v(good), find(bad), 'linear', 'extrap');
        x(:, c) = v;
    end
end

% -------------------------------------------------------------------------
% 2. 只对“非全 NaN”的列滤波
% -------------------------------------------------------------------------
% 这样可选通道缺失时不会让 filtfilt 报错；缺失列原样保留 NaN。
allBadCols = all(~isfinite(x), 1);
y = x;
if all(allBadCols)
    if wasVector
        y = y(:);
    end
    return;
end
xGood = x(:, ~allBadCols);

% -------------------------------------------------------------------------
% 3. 零相位低通
% -------------------------------------------------------------------------
% 正常路径：4 阶 Butterworth + filtfilt。
% nT > 3*order 是 filtfilt 的保守长度要求；太短的数据走移动平均 fallback。
%
% fallback 不是 Butterworth，而是两次居中移动平均：
%   * 仍然零相位（居中）；
%   * 不需要任何工具箱；
%   * 频率响应不如 Butterworth 陡，但足以去掉 Moco 的高频 chatter；
%   * 会给出一次 warning，提示用户结果不是理想 Butterworth。
if exist('butter', 'file') == 2 && exist('filtfilt', 'file') == 2 && nT > 3 * order
    [b, a] = butter(order, fc / (fs / 2), 'low');
    yGood = filtfilt(b, a, xGood);
else
    % 窗口长度按截止频率估计：fs/(2*fc) 个样本约对应半个周期。
    % 强制为奇数，保证窗口中心落在样本上（零相位）。
    win = max(3, 2 * round(fs / (2 * fc)) + 1);
    if win >= nT
        win = max(3, 2 * floor((nT - 1) / 2) + 1);
    end
    yGood = movmean(movmean(xGood, win, 1), win, 1);
    warning('zero_phase_lowpass:fallback', ...
        ['Signal Processing Toolbox (butter/filtfilt) not found; used a ' ...
         'moving-average fallback (window = %d samples).'], win);
end

% 把滤波结果写回非全 NaN 的列；全 NaN 列保持原样。
y(:, ~allBadCols) = yGood;

if wasVector
    y = y(:);
end
end
