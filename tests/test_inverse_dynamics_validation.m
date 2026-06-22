function test_inverse_dynamics_validation()
% test_inverse_dynamics_validation.m
%
% Pipeline:
%   1. 通过 OpenSim API 读取 UN.sto
%   2. 提取逆动力学所需的数据 (coordinate values + 地反力)
%   3. 执行逆动力学，得到 hip_flexion_r / knee_angle_r / ankle_angle_r 关节力矩
%   4. 将 ID 计算结果与 UN.sto 中已有的 .moment 数据画图对比
%
% 依赖: OpenSim 4.x MATLAB API + human0918.osim 模型

import org.opensim.modeling.*

%% ====== 0. 路径与常量 ======
[scriptDir, ~, ~] = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..');

modelPath   = fullfile(projectRoot, 'model', 'human0918.osim');
stoPath     = fullfile(projectRoot, 'Sparse Group LASSO Validation', 'UN.sto');
resultsDir  = fullfile(projectRoot, 'results');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

assert(isfile(modelPath), '模型文件不存在: %s', modelPath);
assert(isfile(stoPath),   'STO  文件不存在: %s', stoPath);

fprintf('========== 逆动力学验证测试 ==========\n');
fprintf('模型 : %s\n', modelPath);
fprintf('运动 : %s\n\n', stoPath);

%% ====== 1. 用 TimeSeriesTable 读取 STO ======
fprintf('[1/4] 读取 STO 文件 ...\n');
tbl = TimeSeriesTable(stoPath);
nFrames = tbl.getNumRows();
nCols   = tbl.getNumColumns();

% 列标签
colLabels = tbl.getColumnLabels();
allLabels = cell(1, nCols);
for i = 0 : nCols - 1
    allLabels{i + 1} = char(colLabels.get(i));
end

% 时间列
time = zeros(nFrames, 1);
indCol = tbl.getIndependentColumn();
for r = 0 : nFrames - 1
    time(r + 1) = indCol.get(r);
end
dt = mean(diff(time));
fprintf('  帧数 = %d, dt = %.5f s, 时长 = %.2f s\n', nFrames, dt, time(end));

% ---- 辅助函数: 按标签读列 ----
    function vec = readCol(label)
        idx = find(strcmp(allLabels, label), 1);
        if isempty(idx)
            error('STO 中找不到列: %s', label);
        end
        dep = tbl.getDependentColumnAtIndex(idx - 1);
        vec = zeros(nFrames, 1);
        for r = 0 : nFrames - 1, vec(r + 1) = dep.get(r); end
    end

    function [vec, found] = tryReadCol(label)
        idx = find(strcmp(allLabels, label), 1);
        if isempty(idx)
            vec = nan(nFrames, 1); found = false; return;
        end
        dep = tbl.getDependentColumnAtIndex(idx - 1);
        vec = zeros(nFrames, 1);
        for r = 0 : nFrames - 1, vec(r + 1) = dep.get(r); end
        found = true;
    end

%% ====== 2. 提取 ID 所需数据 ======
fprintf('[2/4] 提取坐标和地反力数据 ...\n');

% --- 2a. 地反力 (ground reaction forces) ---
% SCONE 命名: leg1_r.grf_x/y/z, leg1_r.grm_x/y/z, leg1_r.cop_x/y/z
%             leg0_l.grf_x/y/z, leg0_l.grm_x/y/z, leg0_l.cop_x/y/z
% grf/cop 在 ground 坐标系中; grm 是 COP 处自由力矩 (不双重计入 force×COP)
legNames = {'leg1_r', 'leg0_l'};
legBodies = {'', ''};          % 模型 body 名称, Step 3 中回填
grfForceVals  = cell(2, 1);   % {Nframes x 3}
grfMomentVals = cell(2, 1);   % {Nframes x 3}
grfCopVals    = cell(2, 1);   % {Nframes x 3}
legHasData    = false(2, 1);

for s = 1 : 2
    ln = legNames{s};
    [~, found] = tryReadCol([ln, '.grf_x']);  %#ok<ASGLU>
    if ~found, continue; end
    legHasData(s) = true;

    grfForceVals{s}  = [readCol([ln, '.grf_x']), readCol([ln, '.grf_y']), readCol([ln, '.grf_z'])];
    grfCopVals{s}    = [readCol([ln, '.cop_x']), readCol([ln, '.cop_y']), readCol([ln, '.cop_z'])];

    [~, hasMx] = tryReadCol([ln, '.grm_x']);
    [~, hasMy] = tryReadCol([ln, '.grm_y']);
    [~, hasMz] = tryReadCol([ln, '.grm_z']);
    if hasMx && hasMy && hasMz
        grfMomentVals{s} = [readCol([ln, '.grm_x']), readCol([ln, '.grm_y']), readCol([ln, '.grm_z'])];
        fprintf('  已提取 %s: GRF + GRM + COP\n', ln);
    else
        grfMomentVals{s} = zeros(nFrames, 3);
        fprintf('  已提取 %s: GRF + COP (无 GRM)\n', ln);
    end
end

% --- 2b. STO 参考力矩 ---
stoMomentLabels = {'hip_flexion_r.moment', 'knee_angle_r.moment', 'ankle_angle_r.moment'};
stoMoments = zeros(nFrames, 3);
for i = 1 : 3
    stoMoments(:, i) = readCol(stoMomentLabels{i});
end
fprintf('  参考力矩: %s, %s, %s\n', stoMomentLabels{:});

%% ====== 3. 执行逆动力学 ======
fprintf('[3/4] 执行逆动力学 ...\n');

% --- 加载模型, 移除肌肉 ---
model = Model(modelPath);

forceSet = model.updForceSet();
removeList = {};
for i = forceSet.getSize() - 1 : -1 : 0
    f = forceSet.get(i);
    clsName = char(f.getConcreteClassName());
    if contains(clsName, 'Muscle')
        removeList{end + 1} = char(f.getName()); %#ok<AGROW>
        forceSet.remove(i);
    end
end
fprintf('  已移除 %d 块肌肉: %s...\n', length(removeList), strjoin(removeList(1:min(3,end)), ', '));

% --- 获取模型坐标 / body 信息 ---
coordSet = model.getCoordinateSet();
nCoord = coordSet.getSize();
coordNameList = cell(1, nCoord);
for i = 0 : nCoord - 1
    coordNameList{i + 1} = char(coordSet.get(i).getName());
end
coordMap = containers.Map(coordNameList, num2cell(0 : nCoord - 1));

bodySet = model.getBodySet();
nBodies = bodySet.getSize();
calcnRFull = ''; calcnLFull = '';
for i = 0 : nBodies - 1
    bn = char(bodySet.get(i).getName());
    if contains(bn, 'calcn') && contains(bn, '_r'), calcnRFull = bn; end
    if contains(bn, 'calcn') && contains(bn, '_l'), calcnLFull = bn; end
end
fprintf('  calcn_r body: "%s", calcn_l body: "%s"\n', calcnRFull, calcnLFull);

legBodies{1} = calcnRFull;   % leg1_r → 右脚
legBodies{2} = calcnLFull;   % leg0_l → 左脚

% --- STO→模型 坐标名映射 (别名处理) ---
sto2modelName = containers.Map();
knownAliases = { ...
    'knee_angle_r',      'knee_flexion_r'; ...
    'knee_angle_l',      'knee_flexion_l'; ...
    'ankle_angle_r',     'ankle_dorsiflexion_r'; ...
    'ankle_angle_l',     'ankle_dorsiflexion_l'};

    function modelName = resolveModelCoord(stoName)
        if sto2modelName.isKey(stoName)
            modelName = sto2modelName(stoName); return;
        end
        if coordMap.isKey(stoName)
            modelName = stoName;
            sto2modelName(stoName) = stoName; return;
        end
        for a = 1 : size(knownAliases, 1)
            if strcmp(stoName, knownAliases{a, 1}) && coordMap.isKey(knownAliases{a, 2})
                modelName = knownAliases{a, 2};
                sto2modelName(stoName) = modelName;
                fprintf('  坐标名映射: STO "%s" → 模型 "%s"\n', stoName, modelName);
                return;
            end
        end
        for k = 0 : nCoord - 1
            mn = coordNameList{k + 1};
            if contains(mn, stoName) || contains(stoName, mn)
                modelName = mn;
                sto2modelName(stoName) = modelName;
                fprintf('  坐标名映射(模糊): STO "%s" → 模型 "%s"\n', stoName, modelName);
                return;
            end
        end
        modelName = '';
    end

% --- 为模型每一个坐标尝试匹配 STO 列 ---
% coordMatch{}.stoName — STO 列名 (q 列)
% coordMatch{}.modelIdx — 模型坐标 0-based index
% coordMatch{}.stoQ / stoU / stoA — Nframes x 1 数据
fprintf('  匹配 STO 列 → 模型坐标 ...\n');
coordMatch = {};
targetJointModelIdx = zeros(1, 3);  % hip, knee, ankle 在模型中的 0-based 索引
targetJointStoName = {'hip_flexion_r', 'knee_angle_r', 'ankle_angle_r'};

for ci = 0 : nCoord - 1
    modelCN = coordNameList{ci + 1};
    
    % 从模型坐标名反向查找 STO 列:
    % 1) 模型名本身在 STO 中?
    [~, foundQ] = tryReadCol(modelCN);
    % 2) 别名反向查找: STO 名 → 模型名, 已知 knee_flexion → knee_angle
    stoCN = '';
    if foundQ
        stoCN = modelCN;
    else
        for a = 1 : size(knownAliases, 1)
            if strcmp(modelCN, knownAliases{a, 2})
                [~, foundAlt] = tryReadCol(knownAliases{a, 1});
                if foundAlt
                    stoCN = knownAliases{a, 1};
                    foundQ = true;
                end
                break;
            end
        end
    end
    if ~foundQ, continue; end
    
    % 读取 q 和 u (速率先从 _u 列读, 缺失则差分 q)
    qVec = readCol(stoCN);
    [uVec, foundU] = tryReadCol([stoCN, '_u']);
    if ~foundU
        uVec = finDiff(qVec, dt);
    end
    aVec = finDiff(uVec, dt);
    
    entry.stoName = stoCN;
    entry.modelIdx = ci;
    entry.stoQ = qVec;
    entry.stoU = uVec;
    entry.stoA = aVec;
    coordMatch{end + 1} = entry;  %#ok<AGROW>
    
    % 记录三个目标关节的索引
    for j = 1 : 3
        mn = resolveModelCoord(targetJointStoName{j});
        if strcmp(modelCN, mn)
            targetJointModelIdx(j) = ci;
        end
    end
end
nMatched = length(coordMatch);
fprintf('  成功匹配 %d / %d 个模型坐标\n', nMatched, nCoord);
for j = 1 : 3
    if targetJointModelIdx(j) < 0
        warning('目标关节 "%s" 未匹配到模型坐标', targetJointStoName{j});
    end
end

% --- 用 PrescribedForce 注入 GRF ---
fprintf('  构建 PrescribedForce 对象 ...\n');
nLegsAdded = 0;
for s = 1 : 2
    if ~legHasData(s), continue; end
    bodyName = legBodies{s};
    if isempty(bodyName), continue; end
    
    Fx = grfForceVals{s}(:, 1);  Fy = grfForceVals{s}(:, 2);  Fz = grfForceVals{s}(:, 3);
    Px = grfCopVals{s}(:, 1);    Py = grfCopVals{s}(:, 2);    Pz = grfCopVals{s}(:, 3);
    % Mx = grfMomentVals{s}(:, 1); My = grfMomentVals{s}(:, 2); Mz = grfMomentVals{s}(:, 3);
    
    pfName = sprintf('grf_%s', legNames{s});
    pf = PrescribedForce(pfName, bodySet.get(bodyName));
    pf.setForceFunctions( buildPLF(Fx, time), buildPLF(Fy, time), buildPLF(Fz, time));
    pf.setPointFunctions( buildPLF(Px, time), buildPLF(Py, time), buildPLF(Pz, time));
    pf.setForceIsInGlobalFrame(true);
    pf.setPointIsInGlobalFrame(true);
    model.addForce(pf);
    nLegsAdded = nLegsAdded + 1;
end
fprintf('  已添加 %d 个 PrescribedForce (GRF+GRM+COP)\n', nLegsAdded);
model.finalizeConnections();

st = model.initSystem();
idSolver = InverseDynamicsSolver(model);

% 预分配 udot Vector (复用, 每帧只改值)
udot = Vector(nCoord, 0.0);

% --- 逐帧逆动力学 ---
fprintf('  逐帧计算 (%d 帧)...\n', nFrames);
idMoments = zeros(nFrames, 3);
hw = waitbar(0, '逆动力学计算中...');

% 找出 pelvis_tx 和 calcn_r body 在匹配列表中的索引, 用于诊断打印
diagPelvisTxIdx = -1;
diagCalcnRBodyName = legBodies{1};  % leg1_r → calcn_r
for m = 1 : nMatched
    if strcmp(coordNameList{coordMatch{m}.modelIdx + 1}, 'pelvis_tx')
        diagPelvisTxIdx = m;
    end
end

% 诊断帧: 第1帧, 中间帧, 末尾帧
diagFrames = [1, round(nFrames / 2), nFrames];

% 诊断: 打印 COP 首尾值, 确认其在全局坐标系
if legHasData(1)
    cop1_start = grfCopVals{1}(1, 1);
    cop1_end   = grfCopVals{1}(end, 1);
    fprintf('  [诊断] leg1_r cop_x 首帧=%.4f, 末帧=%.4f (%s全局系)\n', ...
        cop1_start, cop1_end, ...
        ternary(abs(cop1_end - cop1_start) > 1, '✅ 随前进增大, 确认', '⚠ 未增大, 可能非全局系'));
end

for frame = 1 : nFrames
    if mod(frame, 200) == 0
        waitbar(frame / nFrames, hw, sprintf('帧 %d / %d', frame, nFrames));
    end

    st.setTime(time(frame));

    % 设置所有匹配坐标的 value 和 speed
    % 先全部用 false (不逐次装配), 最后对关键坐标用 true 触发一次性约束投影
    for m = 1 : nMatched
        ci = coordMatch{m}.modelIdx;
        coordSet.get(ci).setValue(st, coordMatch{m}.stoQ(frame), false);
        coordSet.get(ci).setSpeedValue(st, coordMatch{m}.stoU(frame));
    end

    % 去掉 assemble: 对给定运动学的 ID, 不应运行装配优化
    % 直接 realize 到 Velocity (内部先 Position→Velocity, 满足约束但不改动独立坐标)
    model.realizeVelocity(st);

    % ---- 诊断打印 (frame=1/mid/end) ----
    if ismember(frame, diagFrames)
        pelvisSetVal = NaN;
        pelvisStVal  = NaN;
        footGx = NaN;
        copX   = NaN;
        if diagPelvisTxIdx > 0
            pelvisSetVal = coordMatch{diagPelvisTxIdx}.stoQ(frame);
            pelvisStVal  = coordSet.get(coordMatch{diagPelvisTxIdx}.modelIdx).getValue(st);
        end
        if ~isempty(diagCalcnRBodyName)
            try
                footGx = bodySet.get(diagCalcnRBodyName).getPositionInGround(st).get(0);
            catch
            end
        end
        if legHasData(1)
            copX = grfCopVals{1}(frame, 1);
        end
        fprintf('  [诊断] frame=%4d | pelvis_tx 设定=%.4f 状态=%.4f | foot_gx=%.4f cop_x=%.4f\n', ...
            frame, pelvisSetVal, pelvisStVal, footGx, copX);
    end

    % 填充所有匹配坐标的 udot
    for m = 1 : nMatched
        ci = coordMatch{m}.modelIdx;
        udot.set(ci, coordMatch{m}.stoA(frame));
    end

    try
        tau = idSolver.solve(st, udot);
        for j = 1 : 3
            cj = targetJointModelIdx(j);
            if cj >= 0
                idMoments(frame, j) = tau.get(cj);
            end
        end
    catch ME_id
        try
            model.realizeDynamics(st);
            for j = 1 : 3
                cj = targetJointModelIdx(j);
                if cj >= 0
                    idMoments(frame, j) = coordSet.get(cj).getMoment(st);
                end
            end
        catch ME_fb
            if frame == 1
                warning('ID 失败 (帧 %d): %s', frame, ME_fb.message);
                warning('原始错误: %s', ME_id.message);
            end
        end
    end
end

close(hw);
fprintf('  逆动力学计算完成!\n');

%% ====== 4. 画图对比 ======
fprintf('[4/4] 画图对比 ...\n');

jointTitles = {'Hip Flexion R', 'Knee Angle R', 'Ankle Angle R'};
jointTags   = {'hip\_flexion\_r', 'knee\_angle\_r', 'ankle\_angle\_r'};

figure('Color', 'w', 'Position', [50, 50, 1300, 850], 'Name', 'ID Validation');

for j = 1 : 3
    subplot(3, 1, j);

    ref  = stoMoments(:, j);
    calc = idMoments(:, j);

    plot(time, ref,  'b-',  'LineWidth', 1.5); hold on;
    plot(time, calc, 'r--', 'LineWidth', 1.5);

    xlabel('Time (s)', 'FontSize', 11);
    ylabel('Moment (N\cdotm)', 'FontSize', 11);
    title(sprintf('%s — Joint Moment Comparison', jointTitles{j}), 'FontSize', 12);
    legend({'STO Reference', 'ID Computed'}, 'Location', 'best', 'FontSize', 10);
    grid on; box on;

    % RMSE & r
    mask = ~isnan(calc) & ~isnan(ref);
    if sum(mask) > 5
        e = ref(mask) - calc(mask);
        rmse = sqrt(mean(e.^2));
        rho  = corr(ref(mask), calc(mask));
        xl = xlim; yl = ylim;
        text(xl(1) + 0.02 * range(xl), yl(2) - 0.12 * range(yl), ...
            sprintf('RMSE = %.4f N·m    r = %.4f', rmse, rho), ...
            'FontSize', 10, 'BackgroundColor', [1 1 0.85]);
    end
end

sgtitle('Inverse Dynamics Validation: STO Reference vs. ID Computation', ...
    'FontSize', 14, 'FontWeight', 'bold');

%% ====== 保存结果 ======
savePath = fullfile(resultsDir, 'id_validation_results.mat');
save(savePath, 'time', 'stoMoments', 'idMoments', 'jointTags', 'jointTitles');
fprintf('\n结果已保存至: %s\n', savePath);
fprintf('========== 测试完成 ==========\n');

end

%% ====== 辅助函数 ======
function d = finDiff(x, dt)
% 数值微分: 中心差分 + 边界单侧差分
n = length(x);
d = zeros(n, 1);
if n == 1, return; end
if n == 2
    d(:) = (x(2) - x(1)) / dt; return;
end
d(1)     = (x(2) - x(1)) / dt;
d(2:n-1) = (x(3:n) - x(1:n-2)) / (2 * dt);
d(n)     = (x(n) - x(n-1)) / dt;
end

function plf = buildPLF(dataVec, time)
% 构建 PiecewiseLinearFunction; 零/空数据用 Constant(0) 提高效率
    import org.opensim.modeling.*
    if isempty(dataVec) || max(abs(dataVec)) < 1e-12
        plf = Constant(0);
    else
        plf = PiecewiseLinearFunction();
        for r = 1 : length(dataVec)
            plf.addPoint(time(r), dataVec(r));
        end
    end
end

function s = ternary(cond, t, f)
% 三元运算符: cond 为 true 返回 t, 否则 f
    if cond, s = t; else, s = f; end
end
