function F = build_lasso_features_from_kinematics(modelFile, time, q, qdot, ...
    activation, muscleNames, qNames, opts)
%BUILD_LASSO_FEATURES_FROM_KINEMATICS 由 q/qdot/activation 重建 26 个 LASSO 特征。
%
% 这个函数在整个链路中的位置
% --------------------------
% 预测仿真控制器 extract_lasso_reflex_features.m 在运行时看到的是：
%
%   1–9   归一化纤维长度 lCEN（同侧 9 块肌肉）
%   10–18 归一化 MTU 力 fATN = F_fiber_along_tendon / Fmax
%         注意：fATN 对“所有肌肉统一缩放 Fmax”是不变的（分子也按同一比例
%         缩放）；只有 lopt/tendon slack/几何或非均匀 Fmax 改变才会改变它。
%   19–26 pelvis_tilt / hip / knee / ankle 的角度与角速度
%
% MocoInverse 的结果里有 q、qdot、activation，但没有 lCEN/fATN。因此需要
% 在“给定运动 + 给定 activation”的条件下，重建肌肉的纤维长度状态并计算
% 这两个量。这个函数就是做这件事。
%
% 为什么不能简单用“每帧静力平衡”
% ------------------------------
% 一个直觉做法是：每帧把 q 设进去，调用 equilibrateMuscles，然后读 lCEN。
% 但 forward simulation 里的 fiber_length 是一个带记忆的状态，不是每帧
% 重新平衡的。我们在 results/sim_result_*.sto 上实测过：
%   * 准静态平衡近似对 rect_fem / vasti 的 lCEN 误差可达 0.13；
%   * 用 <=2 ms 子步长显式积分 fiber_length 微分方程，误差 <1e-3。
% 所以这里采用子步积分，而不是 equilibrate。
%
% 算法
% ----
%   1. 载入预测仿真模型（默认 human0918.osim）；
%   2. 建立模型坐标/肌肉 handle，并把 qNames/q 映射到模型坐标；
%   3. 第一帧：设 activation、fiber_length=lopt，realizeDynamics，
%      equilibrateMuscles 求初始纤维长度；
%   4. 对每一对相邻输出时刻 [k-1, k]：
%        - 把 dt 切成若干 <= substepMax 的子步；
%        - 子步内线性插值 q/qdot/activation；
%        - 设坐标、activation、当前 fiber_length；
%        - realizeDynamics + computeStateVariableDerivatives；
%        - 只对 fiber_length 做显式 Euler：fl += dtSub * dfl/dt；
%        - clip 到 [fiberLengthMin, fiberLengthMax] * lopt，防止数值爆炸；
%   5. 在真实输出时刻计算 lCEN/fATN。
%
% 输入
% ----
%   modelFile    : 预测仿真模型 .osim（默认 human0918；特征必须和运行时一致）
%   time         : nT×1 [s]
%   q, qdot      : nT×nQ 坐标值/速度（rad 或 m；与 qNames 对应）
%   activation   : nT×nMus 肌肉激活（0..1）
%   muscleNames  : 1×nMus OpenSim 全名，例如 {'hamstrings_r',...}
%   qNames       : 1×nQ 与 q/qdot 列对应的坐标名
%   opts         :
%       .substepMax        默认 0.002 s（100 Hz 数据 -> 5 个子步）
%       .equilibrateInitial 默认 true；失败时退化为 lopt 初值
%       .fiberLengthMin    默认 0.20（相对 lopt）
%       .fiberLengthMax    默认 2.50
%       .verbose           默认 true
%
% 输出 F
% ------
%   .time, .muscleNames
%   .lCEN  nT×nMus  归一化纤维长度
%   .fATN  nT×nMus  归一化 MTU 力（可能为小负值，来自 passive force）
%   .modelFile
%   .qNames
%   .diagnostics   积分子步数、lCEN/fATN 范围、equilibrate 失败次数等
%
% OpenSim API 注意点
% ------------------
% * MATLAB 绑定里坐标/肌肉是 Java 对象，索引 0-based；这里统一用 i-1。
% * 坐标值用 Coordinate.setValue(state, val, false) 设置；false 表示不把它
%   设为 prescribed（不锁定）。
% * 肌肉状态用 model.setStateVariableValue(state, path, value) 按名字设置，
%   避免直接操作 state.getY() 时 internal order 与 label order 不一致。
% * 导数用 getStateVariableDerivativeValue，同样按名字取，避免 permutation。
% * 所有 set 完后必须 realizeDynamics，否则 getFiberForceAlongTendon 可能
%   读不到最新状态。

% -------------------------------------------------------------------------
% 0. 参数默认值与输入尺寸检查
% -------------------------------------------------------------------------
% 这里把所有“列数必须对齐”的检查集中做掉。下游的 Java 调用不会给出
% 友好的尺寸错误，提前检查可以避免 MATLAB 在 setValue 处抛 Java 异常。
if nargin < 8 || isempty(opts)
    opts = struct();
end
opts = apply_defaults(opts);

if ~isfile(modelFile)
    error('build_lasso_features_from_kinematics:noModel', ...
        'Model file not found: %s', modelFile);
end
time = time(:);
nT = numel(time);
if nT < 2
    error('build_lasso_features_from_kinematics:tooShort', 'Need >= 2 samples.');
end
if size(q, 1) ~= nT || size(qdot, 1) ~= nT || size(activation, 1) ~= nT
    error('build_lasso_features_from_kinematics:sizeMismatch', ...
        'time, q, qdot and activation must have the same number of rows.');
end
if size(activation, 2) ~= numel(muscleNames)
    error('build_lasso_features_from_kinematics:muscleMismatch', ...
        'activation has %d columns but muscleNames has %d entries.', ...
        size(activation, 2), numel(muscleNames));
end
if size(q, 2) ~= numel(qNames)
    error('build_lasso_features_from_kinematics:coordMismatch', ...
        'q has %d columns but qNames has %d entries.', size(q, 2), numel(qNames));
end

% -------------------------------------------------------------------------
% 1. 载入模型并建立坐标映射
% -------------------------------------------------------------------------
% 注意：这里只把模型自己拥有的坐标设进去；IK 里可能还有 pelvis_list、
% hip_adduction 等 3D 坐标，human0918 没有，直接忽略。这样同一段代码可以
% 同时服务 2D human0918 和 3D SQR_simbody。
nMus = numel(muscleNames);
model = org.opensim.modeling.Model(char(modelFile));
state = model.initSystem();

% coordHandles : 模型坐标 Java 对象
% coordNamesModel : 模型坐标名
% coordCol : 该模型坐标对应 q/qdot 的哪一列；0 表示 IK 里没有这个坐标
%            （对这种坐标，set_coordinates 会设 0）。
cs = model.getCoordinateSet();
nCoordModel = cs.getSize();
coordHandles = cell(1, nCoordModel);
coordNamesModel = cell(1, nCoordModel);
coordCol = zeros(1, nCoordModel);
for i = 1:nCoordModel
    c = cs.get(i - 1);
    cname = char(c.getName());
    coordNamesModel{i} = cname;
    coordHandles{i} = c;
    j = find(strcmp(qNames, cname), 1);
    if ~isempty(j)
        coordCol(i) = j;
    end
end

% -------------------------------------------------------------------------
% 2. 建立肌肉 handle 与状态变量路径
% -------------------------------------------------------------------------
% 用“模型里的肌肉名”去匹配调用者给的 muscleNames，而不是假设顺序一致。
% 这样即使以后模型 ForceSet 顺序改变，特征列仍按调用者指定顺序输出。
%
% 状态变量路径：
%   /forceset/<muscle>/activation
%   /forceset/<muscle>/fiber_length
% setStateVariableValue / getStateVariableDerivativeValue 都支持这种全路径，
% 不需要自己维护 state.getY() 的 internal/label permutation。
ms = model.getMuscles();
modelMuscleNames = cell(1, ms.getSize());
modelMuscleHandles = cell(1, ms.getSize());
for i = 1:ms.getSize()
    m = ms.get(i - 1);
    modelMuscleNames{i} = char(m.getName());
    modelMuscleHandles{i} = m;
end

muscleHandles = cell(1, nMus);
muscleLongNames = cell(1, nMus);   % full path names used by the model
lopt = zeros(1, nMus);
for i = 1:nMus
    j = find(strcmp(modelMuscleNames, muscleNames{i}), 1);
    if isempty(j)
        error('build_lasso_features_from_kinematics:unknownMuscle', ...
            'Muscle "%s" is not present in %s.', muscleNames{i}, modelFile);
    end
    muscleHandles{i} = modelMuscleHandles{j};
    muscleLongNames{i} = modelMuscleNames{j};
    lopt(i) = muscleHandles{i}.getOptimalFiberLength();
end

actPaths = cell(1, nMus);
flPaths  = cell(1, nMus);
for i = 1:nMus
    actPaths{i} = ['/forceset/' muscleLongNames{i} '/activation'];
    flPaths{i}  = ['/forceset/' muscleLongNames{i} '/fiber_length'];
end

% -------------------------------------------------------------------------
% 3. 第一帧：求初始纤维长度
% -------------------------------------------------------------------------
% 先设 q/qdot -> realizePosition；再把 activation 和 fiber_length=lopt 写进
% 状态；realizeDynamics 后调用 equilibrateMuscles。
%
% 为什么初值用 lopt：
%   * 我们不知道 Moco 解里第一帧的 fiber_length；
%   * equilibration 从 lopt 出发通常能收敛到合理平衡点；
%   * 即使 equilibrate 失败，lopt 也是一个不会让后续积分发散的初值。
% 第一帧附近的 transient 会被接口的 edgeMargin（默认 0.05 s）排除在 LASSO
% 拟合之外。
q0 = q(1, :);
qd0 = qdot(1, :);
set_coordinates(1, q0, qd0);
model.realizePosition(state);

fl = lopt;  % start from optimal fiber length, then equilibrate
for i = 1:nMus
    model.setStateVariableValue(state, actPaths{i}, activation(1, i));
    model.setStateVariableValue(state, flPaths{i},  fl(i));
end
model.realizeDynamics(state);

% equilibrateMuscles 对某些肌肉/构型可能抛 MuscleCannotEquilibrate；
% 这时退化为逐肌肉 computeEquilibrium，仍失败就保留 lopt。不要因为一块
% 肌肉失败就让整个接口失败。
nEquilibrateFailures = 0;
if opts.equilibrateInitial
    try
        model.equilibrateMuscles(state);
        model.realizeDynamics(state);
    catch ME
        nEquilibrateFailures = nEquilibrateFailures + 1;
        warning('build_lasso_features_from_kinematics:equilibrateFailed', ...
            'equilibrateMuscles failed (%s); trying per-muscle equilibrium.', ...
            ME.message);
        for i = 1:nMus
            try
                muscleHandles{i}.computeEquilibrium(state);
            catch
                % keep lopt for this muscle
            end
        end
        try
            model.realizeDynamics(state);
        catch
        end
    end
end
for i = 1:nMus
    v = model.getStateVariableValue(state, flPaths{i});
    if isfinite(v) && v > 0
        fl(i) = v;
    end
end
fl = min(max(fl, opts.fiberLengthMin * lopt), opts.fiberLengthMax * lopt);

lCEN = nan(nT, nMus);
fATN = nan(nT, nMus);
nSubTotal = 0;

% -------------------------------------------------------------------------
% 4. 沿给定运动积分纤维长度
% -------------------------------------------------------------------------
% 关键点：
%   * 每个输出区间内部再分子步，是为了让显式 Euler 稳定；
%     100 Hz 数据 + substepMax=2 ms -> 5 个子步。
%   * 子步内 q/qdot/activation 线性插值；fiber_length 是唯一“积分”的状态，
%     q/qdot/activation 都是 prescribed（每步覆盖）。
%   * dfl/dt 只从模型读；activation 的导数不积分，因为 activation 由 Moco
%     给定，不是我们要预测的量。
%   * 每个子步后 clip fiber_length，避免数值不稳定导致负长度/爆炸。
for k = 1:nT
    if k > 1
        dt = time(k) - time(k - 1);
        nSub = max(1, ceil(dt / opts.substepMax));
        dtSub = dt / nSub;
        nSubTotal = nSubTotal + nSub;

        for s = 1:nSub
            alpha = s / nSub;
            qk  = (1 - alpha) * q(k - 1, :)  + alpha * q(k, :);
            qdk = (1 - alpha) * qdot(k - 1, :) + alpha * qdot(k, :);
            ak  = (1 - alpha) * activation(k - 1, :) + alpha * activation(k, :);

            set_coordinates(k, qk, qdk);
            model.realizePosition(state);
            for mm = 1:nMus
                model.setStateVariableValue(state, actPaths{mm}, ak(mm));
                model.setStateVariableValue(state, flPaths{mm},  fl(mm));
            end
            model.realizeDynamics(state);

            % computeStateVariableDerivatives 只填充 state.getYDot()；
            % 这里用 getStateVariableDerivativeValue 按名字取 fiber_length
            % 的导数，避开 internal/label order 的坑。
            model.computeStateVariableDerivatives(state);
            for mm = 1:nMus
                dfl = model.getStateVariableDerivativeValue(state, flPaths{mm});
                if isfinite(dfl)
                    fl(mm) = fl(mm) + dtSub * dfl;
                end
            end
            fl = min(max(fl, opts.fiberLengthMin * lopt), opts.fiberLengthMax * lopt);
        end
    end

    % 记录真实输出时刻（不是子步时刻）：
    % 子步只用于推进 fiber_length；最终特征必须在 tOut 上取值，才能和
    % target/相位/权重逐样本对齐。
    set_coordinates(k, q(k, :), qdot(k, :));
    model.realizePosition(state);
    for mm = 1:nMus
        model.setStateVariableValue(state, actPaths{mm}, activation(k, mm));
        model.setStateVariableValue(state, flPaths{mm},  fl(mm));
    end
    model.realizeDynamics(state);

    for mm = 1:nMus
        lCEN(k, mm) = muscleHandles{mm}.getNormalizedFiberLength(state);
        fATN(k, mm) = muscleHandles{mm}.getFiberForceAlongTendon(state) ...
                      / muscleHandles{mm}.getMaxIsometricForce();
    end
end

F = struct();
F.time        = time;
F.muscleNames = muscleNames;
F.lCEN        = lCEN;
F.fATN        = fATN;
F.modelFile   = char(modelFile);
F.qNames      = qNames(:).';
F.diagnostics = struct( ...
    'nEquilibrateFailures', nEquilibrateFailures, ...
    'meanSubstepsPerFrame', nSubTotal / max(1, nT - 1), ...
    'lCENmin', min(lCEN(:)), 'lCENmax', max(lCEN(:)), ...
    'fATNmin', min(fATN(:)), 'fATNmax', max(fATN(:)), ...
    'fiberLengthMinFactor', opts.fiberLengthMin, ...
    'fiberLengthMaxFactor', opts.fiberLengthMax);

if opts.verbose
    fprintf(['[features] %s: lCEN [%.3f, %.3f], fATN [%.3f, %.3f], ' ...
             'mean substeps/frame = %.1f\n'], ...
        char(modelFile), F.diagnostics.lCENmin, F.diagnostics.lCENmax, ...
        F.diagnostics.fATNmin, F.diagnostics.fATNmax, ...
        F.diagnostics.meanSubstepsPerFrame);
end

    % ---- nested helper: prescribe q / qdot for one frame ----------------
    % 嵌套函数共享父函数的 state/coordHandles/coordCol，避免每帧重复查找
    % Java 对象。k 只用于可读性/调试，不参与索引（因此标 INUSD）。
    % 模型里有、但 qNames 里没有的坐标统一设 0：对 human0918 只有 9 个
    % 坐标，全部都能匹配；对 3D 模型则不会漏设坐标。
    function set_coordinates(k, qRow, qdRow) %#ok<INUSD>
        for cc = 1:nCoordModel
            h = coordHandles{cc};
            if coordCol(cc) > 0
                h.setValue(state, qRow(coordCol(cc)), false);
                h.setSpeedValue(state, qdRow(coordCol(cc)));
            else
                h.setValue(state, 0, false);
                h.setSpeedValue(state, 0);
            end
        end
    end
end

% ========================================================================
% 局部函数：默认参数
% ========================================================================
% substepMax 是最关键参数：
%   * 0.002 s 在 100 Hz 数据上给 5 个子步，验证误差 <1e-3；
%   * 调大（如 0.005）会更快，但 rect_fem/vasti 的 fiber_length 误差会增大；
%   * 调小（如 0.001）更准，但 MATLAB Java 调用次数线性增加。
% fiberLengthMin/Max 是数值安全边界，不是生理模型参数；正常步态 lCEN
% 大约 0.5–1.4，2.5 已经非常宽松。
function opts = apply_defaults(opts)
defaults = struct( ...
    'substepMax', 0.002, ...
    'equilibrateInitial', true, ...
    'fiberLengthMin', 0.20, ...
    'fiberLengthMax', 2.50, ...
    'verbose', true);
fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}) || isempty(opts.(fn{i}))
        opts.(fn{i}) = defaults.(fn{i});
    end
end
end
