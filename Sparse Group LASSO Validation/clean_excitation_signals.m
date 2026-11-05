function [xClean, aClean, qc] = clean_excitation_signals(time, xRaw, aRaw, opts)
%CLEAN_EXCITATION_SIGNALS 把 MocoInverse 的肌肉驱动清洗成 LASSO 可用目标。
%
% 背景：为什么实测数据接入 MocoInverse 后需要这一步
% ------------------------------------------------
% MocoInverse 求的是“在满足动力学/实测运动学/GRF 的前提下”的肌肉控制与
% 激活轨迹。由于实测数据本身有噪声、优化网格有限、目标函数只惩罚 excitation
% 幅值而不惩罚变化率，解里常见三类问题：
%
%   (i)  高频 chatter：collocation 网格尺度的交替分量；
%   (ii) 短 burst/spike：优化器利用 activation 快速通道（tau_act=10 ms）
%        产生的非生理尖峰；
%   (iii) excitation x 与 activation a 只是“近似”满足
%        da/dt = (x-a)/tau_a，而不是严格满足。
%
% LASSO 使用的是平方误差（外加 censoring hinge），对少数大幅值残差非常
% 敏感：几个 burst 就能把系数向量拉偏，导致学到的反射律在预测仿真里
% 产生抖动。因此本函数提供一套透明、可配置的清洗流程，并把每一步的
% 统计量放进 qc，方便做消融实验。
%
% 清洗流程（按执行顺序）
% ----------------------
%   1. 非有限值按列线性插值（全缺失列保留 NaN，由调用者判断是否致命）；
%   2. 可选 Hampel/中值 despike：只替换长度 <= maxBurstSamples 的短 run，
%      持续的真实峰值保留；用 despikeMinSigma 防止“0.01 平台导致鲁棒
%      尺度退化为 0，从而把所有正常波动都判成 spike”；
%   3. 零相位 Butterworth 低通（默认 6 Hz），避免引入相位延迟；
%   4. 用清洗后的 raw excitation 作为 LASSO 回归 target；
%   5. 裁剪到 [floor, ceiling]，默认 [0.01, 1]，其中 floor 必须与 LASSO
%      的 censoring floor c 一致，否则 censored 样本的统计假设会错位。
%
% 重要：target 只来自 excitation，不再从 activation 反演
% ------------------------------------------------------
% 正式 MocoInverse 结果一定包含 /forceset/<muscle> 的 excitation control。
% 因此本函数**不再**提供 activation-consistent 之类的回退方案：
%   * 如果 xRaw 为空，直接报错；
%   * 如果某个肌肉的 excitation 整列缺失，直接报错；
%   * aRaw 只用于计算 fATN/lCEN 等特征，不参与构造 target。
% 如果 aRaw 缺失，可以用 xRaw 积分出 activation 供特征使用，但这只是
% “补特征所需的 state”，不是“自己求解 excitation”。
%
% 输入
% ----
%   time  : nT×1 [s]，严格递增；内部用 median(diff(time)) 估计采样率；
%   xRaw  : nT×nMus raw excitation，**必须提供**（空矩阵会报错）；
%   aRaw  : nT×nMus raw activation；只用于特征计算，缺失时可由 xRaw 积分；
%   opts  : struct，所有字段可选：
%       .lowpassHz         默认 6；0 表示关闭
%       .filterOrder       默认 4
%       .doDespike         默认 true
%       .medianWindow      默认 5 个样本（100 Hz 下 50 ms）
%       .hampelSigma       默认 4（4×MAD）
%       .despikeMinSigma   默认 0.005（绝对噪声地板）
%       .maxBurstSamples   默认 5（只替换不超过 5 个样本的 run）
%       .floor             默认 0.01
%       .ceiling           默认 1.0
%       .tauAct            默认 0.01 s（DGF/Millard 激活时间常数）
%       .tauDeact          默认 0.04 s（去激活时间常数）
%
% 输出
% ----
%   xClean : nT×nMus 清洗后的 excitation target；必须全部有限，且落在
%            [floor, ceiling]；
%   aClean : nT×nMus 清洗后的 activation，用于 build_lasso_features 计算
%            lCEN / fATN（不参与 target 构造）；
%   qc     : 每块肌肉的统计量、flag 矩阵；qc.targetSource 固定为
%            'excitation'。
%
% 重要假设
% --------
% * time 是均匀网格；Moco/IK 插值到 targetDt 后满足；
% * x 和 a 的列顺序必须与 muscleNames 一致（本仓库为 canonical right 9 +
%   left 9）；本函数只按列处理，不解释肌肉名字；
% * xRaw 必须包含所有肌肉的 excitation；缺失任何一个通道都会报错，
%   不再尝试从 activation 反演。

% -------------------------------------------------------------------------
% 0. 默认参数、时间轴与采样率
% -------------------------------------------------------------------------
% 至少 4 个样本：filtfilt/差分/中值滤波在这个长度以下都没有意义；
% dt 用 median(diff(time)) 而不是 time(2)-time(1)，对轻微不均匀网格更稳。
if nargin < 4
    opts = struct();
end
opts = apply_defaults(opts);
% 兼容性检查：旧的 activation_consistent / auto 配置不再支持。
% target 必须来自 raw excitation；如果调用者还传旧字段，给出明确错误。
if isfield(opts, 'targetSource') && ~strcmpi(opts.targetSource, 'excitation')
    error('clean_excitation_signals:badTargetSource', ...
        ['Only targetSource=''excitation'' is supported. ' ...
         'activation_consistent / activation / auto have been removed: ' ...
         'formal MocoInverse results must provide excitation.']);
end

time = time(:);
nT = numel(time);
if nT < 4
    error('clean_excitation_signals:tooShort', 'Need at least 4 samples.');
end
dt = median(diff(time));
if ~(isfinite(dt) && dt > 0)
    error('clean_excitation_signals:badTime', 'time must be strictly increasing.');
end
fs = 1 / dt;

% xRaw 是必需的：正式 MocoInverse 一定有 excitation control。
% 这里直接报错，不再尝试从 activation 反演。
if isempty(xRaw)
    error('clean_excitation_signals:noExcitation', ...
        ['MocoInverse excitation is required, but xRaw is empty. ' ...
         'The formal MocoInverse solution must contain ' ...
         '"/forceset/<muscle>" excitation columns; this function no longer ' ...
         'reconstructs excitation from activation.']);
end

hasX = true;
hasA = ~isempty(aRaw);

% -------------------------------------------------------------------------
% 1. 输入类型/尺寸检查 + 非有限值处理
% -------------------------------------------------------------------------
% 对 xRaw：要求每一列至少有一个有限值；整列缺失说明该肌肉没有 excitation，
% 直接报错（而不是让后续低通/裁剪产生 NaN）。
if size(xRaw, 1) ~= nT
    error('clean_excitation_signals:sizeMismatch', ...
        'xRaw must have %d rows.', nT);
end
xRaw = double(xRaw);
missingX = all(~isfinite(xRaw), 1);
if any(missingX)
    error('clean_excitation_signals:missingExcitationChannel', ...
        ['Excitation is missing for %d muscle column(s) (all-NaN). ' ...
         'Formal MocoInverse results must provide excitation for every muscle; ' ...
         'refusing to reconstruct it from activation.'], nnz(missingX));
end
xRaw = fill_nans(xRaw);

if hasA
    aRaw = double(aRaw);
    if size(aRaw, 1) ~= nT
        error('clean_excitation_signals:sizeMismatch', ...
            'aRaw must have %d rows.', nT);
    end
    aRaw = fill_nans(aRaw);
end

% 列数必须一致；nMus 由 excitation 决定（target 是 excitation）。
nMus = size(xRaw, 2);
if hasA && size(aRaw, 2) ~= nMus
    error('clean_excitation_signals:muscleMismatch', ...
        'xRaw and aRaw must have the same number of columns.');
end

% -------------------------------------------------------------------------
% 2. 如果只有 excitation，先用激活动力学积分出 activation
% -------------------------------------------------------------------------
% 这一步只是为了让后续特征计算（lCEN/fATN，需要 activation）有输入；
% 它不改变 xRaw，也不代表我们“相信”这个积分结果优于 Moco 的 activation。
% 一阶动力学：da/dt = (x-a)/tau，上升/下降用不同时间常数。
if ~hasA
    aRaw = integrate_activation(xRaw, time, opts.tauAct, opts.tauDeact);
end

% -------------------------------------------------------------------------
% 3. activation 的 despike + 低通
% -------------------------------------------------------------------------
% 先 despike 再低通：如果先低通，spike 会被抹成较宽的“鼓包”，反而更难
% 用 Hampel 检测。despike 只替换短 run，不会削掉持续的真实峰值。
% 最后 clip 到 [0,1]：低通在 0/1 附近可能产生轻微过冲，物理上 activation
% 不能为负或大于 1。
if opts.doDespike
    [aProc, flagA] = despike_columns(aRaw, opts.medianWindow, ...
        opts.hampelSigma, opts.despikeMinSigma, opts.maxBurstSamples);
else
    aProc = aRaw;
    flagA = false(size(aRaw));
end
aClean = zero_phase_lowpass(aProc, fs, opts.lowpassHz, opts.filterOrder);
aClean = min(max(aClean, 0), 1);

% -------------------------------------------------------------------------
% 4. 构造 LASSO 回归 target
% -------------------------------------------------------------------------
% flagX 只用于 QC：记录 raw x 里被判为短 spike 的样本，方便在 QC 图/CSV
% 里报告“原始 excitation 有多脏”。
flagX = false(nT, nMus);
xProc = [];
if hasX
    if opts.doDespike
        [xProc, flagX] = despike_columns(xRaw, opts.medianWindow, ...
            opts.hampelSigma, opts.despikeMinSigma, opts.maxBurstSamples);
    else
        xProc = xRaw;
    end
end

% target 只有一个来源：清洗后的 raw excitation。
% 不再有 activation / activation_consistent 分支。
xClean = zero_phase_lowpass(xProc, fs, opts.lowpassHz, opts.filterOrder);
targetUsed = 'excitation';

if any(~isfinite(xClean(:)))
    error('clean_excitation_signals:missingExcitationChannel', ...
        ['Excitation is required for every muscle; at least one channel ' ...
         'is still non-finite after interpolation.']);
end

% 统一裁剪。floor 必须与 fit_adaptive_lasso_phase_model 的 options.c 一致；
% 否则“y_i <= c 视为 censored”的假设会和 target 的实际下限错位。
xClean = min(max(xClean, opts.floor), opts.ceiling);

% ---- quality control -----------------------------------------------------
aPred = integrate_activation(xClean, time, opts.tauAct, opts.tauDeact);
reconRmse = sqrt(mean((aPred - aClean).^2, 1)).';

nRawX   = nan(1, nMus);
meanRawX = nan(1, nMus);
maxRawX  = nan(1, nMus);
if hasX
    nRawX    = sum(isfinite(xRaw), 1);
    meanRawX = mean(xRaw, 1, 'omitnan');
    maxRawX  = max(xRaw, [], 1, 'omitnan');
end

qc = struct();
qc.targetSource          = targetUsed;
qc.hasRawExcitation      = hasX;
qc.hasRawActivation      = hasA;
qc.fs                    = fs;
qc.dt                    = dt;
qc.activationReconRMSE   = reconRmse;
qc.flaggedExcitation     = flagX;
qc.flaggedActivation     = flagA;
qc.nFlaggedExcitation    = sum(flagX, 1);
qc.pctFlaggedExcitation  = 100 * mean(flagX, 1);
qc.cleanExcitationMin    = min(xClean, [], 1);
qc.cleanExcitationMax    = max(xClean, [], 1);
qc.cleanExcitationMean   = mean(xClean, 1);
qc.cleanActivationMin    = min(aClean, [], 1);
qc.cleanActivationMax    = max(aClean, [], 1);
qc.rawExcitationMean     = meanRawX;
qc.rawExcitationMax      = maxRawX;
qc.nRawExcitationSamples = nRawX;

% qc.stats 是一个 nMus×6 的紧凑数值矩阵，列含义见 statNames。
% 单独给矩阵是为了让接口直接 array2table/writetable 生成 QC CSV；
% 用 cell/struct 存每块肌肉的统计会让 CSV 写出更啰嗦。
% 列顺序：
%   mean_raw_x          原始 excitation 均值（缺失通道为 NaN）
%   max_raw_x           原始 excitation 峰值
%   mean_clean_x        清洗后 target 均值
%   max_clean_x         清洗后 target 峰值
%   pct_flagged         被判为短 spike 的样本百分比
%   activation_recon_rmse  由 xClean 积分回 activation 与 aClean 的 RMSE
qc.statNames = {'mean_raw_x', 'max_raw_x', 'mean_clean_x', ...
                'max_clean_x', 'pct_flagged', 'activation_recon_rmse'};
qc.stats = [meanRawX(:), maxRawX(:), mean(xClean, 1).', ...
            max(xClean, [], 1).', 100 * mean(flagX, 1).', reconRmse];
end

% ========================================================================
% 局部函数 1：填默认参数
% ========================================================================
% 采用“用户字段优先、缺省/空值才填默认”的策略：
% * 空字符串/空矩阵也视为未设置，避免调用者传 struct() 时留下空字段；
% * 不删除用户多给的字段，方便未来扩展；
% * 所有默认值都在这里集中，文档和代码不会分叉。
function opts = apply_defaults(opts)
defaults = struct( ...
    'lowpassHz', 6, ...
    'filterOrder', 4, ...
    'doDespike', true, ...
    'medianWindow', 5, ...
    'hampelSigma', 4, ...
    'despikeMinSigma', 0.005, ...
    'maxBurstSamples', 5, ...
    'floor', 0.01, ...
    'ceiling', 1.0, ...
    'tauAct', 0.01, ...
    'tauDeact', 0.04);
fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}) || isempty(opts.(fn{i}))
        opts.(fn{i}) = defaults.(fn{i});
    end
end
end

% ========================================================================
% 局部函数 2：按列线性插值非有限值
% ========================================================================
% 只处理部分缺失的列；整列全 NaN 原样保留，这是刻意的：
%   * 对 activation，整列缺失后面会由接口报 MissingActivation；
%   * 对可选 excitation，整列缺失是合法的（例如只导出了 5 块肌肉）；
% 如果在这里无条件报错，接口就无法处理部分列文件。
function v = fill_nans(v)
bad = ~isfinite(v);
if ~any(bad(:))
    return;
end
nT = size(v, 1);
for c = 1:size(v, 2)
    col = v(:, c);
    b = bad(:, c);
    if all(b)
        % 整列全缺失：保持 NaN。
        % 对 xRaw，调用前已经显式检查并报错，不会走到这里；
        % 对 aRaw，整列缺失会让后续特征计算报 MissingActivation。
        continue;
    end
    good = find(~b);
    col(b) = interp1(good, col(good), find(b), 'linear', 'extrap');
    v(:, c) = col;
end
end

% ========================================================================
% 局部函数 3：Hampel/中值 despike（只替换短 run）
% ========================================================================
% 原理：
%   1) 用滑动中值 med 作为局部“正常水平”；
%   2) 残差 res = x - med；
%   3) 用 1.4826*MAD 估计残差尺度（对离群点稳健）；
%   4) |res| > nSigma*尺度 判为 spike；
%   5) 只替换连续长度 <= maxRun 的 spike run，持续高值保留。
%
% 为什么需要 minSigma：
%   activation 大量停在 0.01 下限，导致 res 的一半以上是 0，MAD=0。
%   若不加地板，任何微小波动都会满足 |res| > 4*0，flag 率可达 20–30%。
%   minSigma=0.005 表示“小于 0.5% activation 的偏差不当作 spike”。
%
% 输入：
%   x         nT×nMus
%   win       中值窗口（样本数，内部强制为奇数）
%   nSigma    Hampel 阈值倍数
%   minSigma  鲁棒尺度的绝对下限
%   maxRun    只替换长度不超过该值的连续 run
%
% 输出：
%   y         替换后的信号
%   flag      nT×nMus logical，标记“原始信号中被判为 spike”的样本
%             （即使该 run 太长、最终没有替换，flag 仍然为 true，便于 QC）
function [y, flag] = despike_columns(x, win, nSigma, minSigma, maxRun)
%DESPIKE_COLUMNS  Hampel filter that only replaces short spike runs.
y = x;
flag = false(size(x));
nT = size(x, 1);
% 窗口太小/太大都没有意义，直接不处理。
% movmedian 要求窗口为奇数才能关于当前样本对称；偶数时 +1。
if win < 3
    return;
end
if mod(win, 2) == 0
    win = win + 1;
end
if win > nT
    win = nT - 1 + mod(nT, 2);   % largest odd <= nT
end
if win < 3
    return;
end

for c = 1:size(x, 2)
    % 每列独立处理：不同肌肉的幅值/噪声水平不同，不能共享一个尺度。
    col = x(:, c);
    med = movmedian(col, win, 'omitnan');
    res = col - med;
    s = 1.4826 * median(abs(res - median(res)));
    if ~isfinite(s) || s <= 0
        s = max(eps, 1e-9);
    end
    % Many activation samples sit on the 0.01 floor, making the robust scale
    % collapse to zero.  A small absolute floor prevents flagging every normal
    % oscillation as a spike.
    s = max(s, minSigma);
    f = abs(res) > nSigma * s;
    if ~any(f)
        continue;
    end
    runs = logical_runs(f);
    replace = false(nT, 1);
    for r = 1:size(runs, 1)
        if (runs(r, 2) - runs(r, 1) + 1) <= maxRun
            replace(runs(r, 1):runs(r, 2)) = true;
        end
    end
    y(replace, c) = med(replace);
    flag(:, c) = f;
end
end

% ========================================================================
% 局部函数 4：把 logical mask 转成 [start, end] run 列表
% ========================================================================
% 两端各补一个 false，这样第一段/最后一段也能被 diff 捕获：
%   d = 1  -> run 开始
%   d = -1 -> run 结束（end = 位置-1）
% 输出 runs 是 nRun×2，列为 1-based 样本索引。
function runs = logical_runs(mask)
mask = logical(mask(:));
d = diff([false; mask; false]);
starts = find(d == 1);
ends   = find(d == -1) - 1;
runs = [starts(:), ends(:)];
end

% ========================================================================
% 局部函数 5：由 excitation 正向积分 activation
% ========================================================================
% 离散格式（显式 Euler，x 取新时刻）：
%   a_k = a_{k-1} + dt/tau * (x_k - a_{k-1})
%   tau = tauAct   if x_k >= a_{k-1}
%         tauDeact if x_k <  a_{k-1}
% 这是 Moco/Millard 一阶激活动力学的常用离散形式，和 Moco 的
% tau_act=10 ms / tau_deact=40 ms 默认值一致。
%
% 边界：
%   * a(1)=x(1)，与 MocoInverse 的 a(0)=x(0) 约束一致；
%   * 每步 clip 到 [0,1]，防止显式积分在陡峭输入下过冲；
%   * 只用它来“补 activation”，不用于生成 Moco 的 x。
function a = integrate_activation(x, t, tauAct, tauDeact)
% Forward integration of da/dt = (x-a)/tau with x at the new time level.
nT = size(x, 1);
nMus = size(x, 2);
a = zeros(nT, nMus);
a(1, :) = min(max(x(1, :), 0), 1);
for k = 2:nT
    dt = t(k) - t(k - 1);
    tau = tauAct * ones(1, nMus);
    tau(x(k, :) < a(k - 1, :)) = tauDeact;
    a(k, :) = a(k - 1, :) + (dt ./ tau) .* (x(k, :) - a(k - 1, :));
    a(k, :) = min(max(a(k, :), 0), 1);
end
end


