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
% 实现约定
% --------
% * 使用 Signal Processing Toolbox 的 butter + filtfilt；
% * 逐列调用 filtfilt，保持通道之间互不影响；
% * 输入必须是有限值：上游 clean_excitation_signals / 特征重建已经
%   保证这一点；如果仍有 NaN/Inf，直接报错，不再插值或静默跳过；
% * 不再提供“无 Signal Processing Toolbox”的移动平均 fallback。
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

% -------------------------------------------------------------------------
% 0. 参数校验
% -------------------------------------------------------------------------
if nargin < 3 || isempty(fs) || isempty(fc)
    error('zero_phase_lowpass:badInput', 'x, fs and fc are required.');
end
if nargin < 4 || isempty(order)
    order = 4;
end

% 记录输入是否是向量：向量输入时先转成 nT×1，输出再转回列向量。
% 放在 early return 之前，保证 fc<=0 / fc>=Nyquist 时也遵守这个约定。
wasVector = isvector(x);
if wasVector
    x = x(:);
end

% 输入必须是有限值；不管后面是否真的滤波，都先检查，避免 fc<=0 或
% fc>=Nyquist 时把 NaN 原样返回。上游已经保证有限，这里只是显式防线。
if any(~isfinite(x(:)))
    error('zero_phase_lowpass:nonFiniteInput', ...
        ['x contains NaN/Inf. Upstream cleaning/reconstruction must ' ...
         'provide finite signals before low-pass filtering.']);
end

% 这些情况下不做任何滤波，直接返回原输入：
% * 空矩阵；
% * fc<=0：用户显式关闭滤波；
% * fc/fs 非有限或 fs<=0：参数无效，交给调用者检查；
% * fc>=Nyquist：没有可滤的频率。
if isempty(x) || fc <= 0 || ~isfinite(fc) || ~isfinite(fs) || fs <= 0
    y = x;
    return;
end
if fc >= fs / 2
    y = x;
    return;
end

% filtfilt 对短信号有最小长度要求（约 3*order）；太短就明确报错，
% 不再用移动平均 fallback 掩盖问题。
nT = size(x, 1);
if nT <= 3 * order
    error('zero_phase_lowpass:tooShort', ...
        'Need more than 3*order = %d samples for filtfilt, got %d.', ...
        3 * order, nT);
end

% -------------------------------------------------------------------------
% 1. 设计滤波器
% -------------------------------------------------------------------------
% fc 归一化到 Nyquist；fc < fs/2 已在上面的检查中保证。
[b, a] = butter(order, fc / (fs / 2), 'low');

% -------------------------------------------------------------------------
% 2. 逐列零相位滤波
% -------------------------------------------------------------------------
% 逐列而不是一次 filtfilt(b,a,x)：通道数不多，循环更直观，也避免不同
% 通道的边界条件互相影响。
y = zeros(size(x));
for c = 1:size(x, 2)
    y(:, c) = filtfilt(b, a, x(:, c));
end

if wasVector
    y = y(:);
end
end
