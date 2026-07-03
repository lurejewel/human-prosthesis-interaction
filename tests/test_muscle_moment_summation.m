function test_muscle_moment_summation()
% test_muscle_moment_summation.m
%
% Pipeline:
%   1. 通过 OpenSim API 读取 UN.sto
%   2. 对右侧腿, 读取所有肌肉绕各个 coordinate 的 moment
%      (格式: <muscle>_r.<coord>.moment, 例如 rect_fem_r.knee_extension_r.moment)
%   3. 按 coordinate 加总肌肉力矩, 膝关节减去被动约束力 (knee_r.torque)
%   4. 与 UN.sto 中的净关节力矩 (hip_flexion_r.moment 等) 画图对比
%
% 依赖: OpenSim 4.x MATLAB API

import org.opensim.modeling.*

%% ====== 0. 路径 ======
[scriptDir, ~, ~] = fileparts(mfilename('fullpath'));
projectRoot = fullfile(scriptDir, '..');

stoPath     = fullfile(projectRoot, 'Sparse Group LASSO Validation', 'UN.sto');
resultsDir  = fullfile(projectRoot, 'results');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

assert(isfile(stoPath), 'STO 文件不存在: %s', stoPath);

fprintf('========== 肌肉力矩求和验证 ==========\n');
fprintf('运动 : %s\n\n', stoPath);

%% ====== 1. 用 TimeSeriesTable 读取 STO ======
fprintf('[1/4] 读取 STO 文件 ...\n');
tbl = TimeSeriesTable(stoPath);
nFrames = tbl.getNumRows();
nCols   = tbl.getNumColumns();

colLabels = tbl.getColumnLabels();
allLabels = cell(1, nCols);
for i = 0 : nCols - 1
    allLabels{i + 1} = char(colLabels.get(i));
end

time = zeros(nFrames, 1);
indCol = tbl.getIndependentColumn();
for r = 0 : nFrames - 1
    time(r + 1) = indCol.get(r);
end
fprintf('  帧数 = %d, 时长 = %.2f s\n', nFrames, time(end));

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

%% ====== 2. 搜索并读取所有右侧肌肉绕各关节的 moment ======
fprintf('[2/4] 搜索右侧肌肉 moment 列 ...\n');

% 右侧肌肉列表 (human0918 模型中存在的肌肉, _r 后缀)
muscleNames = {'hamstrings_r', 'bifemsh_r', 'glut_max_r', 'iliopsoas_r', ...
               'rect_fem_r', 'vasti_r', 'gastroc_r', 'soleus_r', 'tib_ant_r'};

% 目标关节坐标
targetCoords = {'hip_flexion_r', 'knee_extension_r', 'ankle_dorsiflexion_r'};

% 搜索所有 <muscle>_r.<coord>.moment 列
% momentMap{coordIdx} = Nframes x nM 矩阵, 每列是一个肌肉的 moment
% momentLabels{coordIdx} = cell array of column labels
sumMoments = zeros(nFrames, 3);  % 按坐标累加
momentDetail = cell(3, 1);       % 各坐标下每个肌肉的贡献明细
momentLabels = cell(3, 1);

for j = 1 : 3
    coord = targetCoords{j};
    foundCols = {};
    foundMusc = {};
    for m = 1 : length(muscleNames)
        label = [muscleNames{m}, '.', coord, '.moment'];
        if any(strcmp(allLabels, label))
            foundCols{end + 1} = label;           %#ok<AGROW>
            foundMusc{end + 1} = muscleNames{m};  %#ok<AGROW>
        end
    end
    
    if isempty(foundCols)
        warning('  未找到任何肌肉对 %s 的 moment 列', coord);
        continue;
    end
    
    nFound = length(foundCols);
    detail = zeros(nFrames, nFound);
    for k = 1 : nFound
        detail(:, k) = readCol(foundCols{k});
    end
    sumMoments(:, j) = sum(detail, 2);
    momentDetail{j} = detail;
    momentLabels{j} = foundCols;
    
    fprintf('  %s: %d 块肌肉贡献', coord, nFound);
    for k = 1 : nFound
        fprintf(', %s', foundMusc{k});
    end
    fprintf('\n');
end

%% ====== 3. 读取 STO 中的净关节力矩 ======
fprintf('[3/4] 读取 STO 净关节力矩 ...\n');

stoMomentLabels = {'hip_flexion_r.moment', 'knee_extension_r.moment', 'ankle_dorsiflexion_r.moment'};
stoMoments = zeros(nFrames, 3);
for j = 1 : 3
    stoMoments(:, j) = readCol(stoMomentLabels{j});
end
fprintf('  已读取: %s\n', strjoin(stoMomentLabels, ', '));

% 读取膝关节被动约束力
kneePassiveTorque = readCol('knee_r.torque');
fprintf('  已读取: knee_r.torque (膝关节被动约束力)\n');

% 修正后的肌肉力矩总和: 膝关节减去被动约束力
sumMomentsAdjusted = sumMoments;
sumMomentsAdjusted(:, 2) = sumMoments(:, 2) - kneePassiveTorque;

%% ====== 4. 画图对比 ======
fprintf('[4/4] 画图对比 ...\n');

jointTitles = {'Hip Flexion R', 'Knee Angle R', 'Ankle Angle R'};

% ---- 图1: 净力矩 vs 肌肉力矩总和 (膝关节已减去被动约束力 knee_r.torque) ----
figure('Color', 'w', 'Position', [50, 50, 1300, 850], 'Name', 'Muscle Moment Summation');

for j = 1 : 3
    subplot(3, 1, j);
    
    plot(time, stoMoments(:, j), 'b-', 'LineWidth', 1.5); hold on;
    plot(time, sumMomentsAdjusted(:, j), 'r--', 'LineWidth', 1.5);
    
    xlabel('Time (s)', 'FontSize', 11);
    ylabel('Moment (N\cdotm)', 'FontSize', 11);
    if j == 2
        title(sprintf('%s — Net Joint Moment vs Sum (Muscle - Passive Torque)', jointTitles{j}), ...
            'FontSize', 12);
    else
        title(sprintf('%s — Net Joint Moment vs Sum of Muscle Moments', jointTitles{j}), ...
            'FontSize', 12);
    end
    legend({'STO Net Moment', 'Sum (Adjusted)'}, 'Location', 'best', 'FontSize', 10);
    grid on; box on;
    
    mask = ~isnan(sumMomentsAdjusted(:, j)) & ~isnan(stoMoments(:, j));
    if sum(mask) > 5
        e = stoMoments(mask, j) - sumMomentsAdjusted(mask, j);
        rmse = sqrt(mean(e.^2));
        rho  = corr(stoMoments(mask, j), sumMomentsAdjusted(mask, j));
        xl = xlim; yl = ylim;
        text(xl(1) + 0.02 * range(xl), yl(2) - 0.12 * range(yl), ...
            sprintf('RMSE = %.4f N·m    r = %.4f', rmse, rho), ...
            'FontSize', 10, 'BackgroundColor', [1 1 0.85]);
    end
end

sgtitle('Muscle Moment Summation: Net Joint Moment vs Sum of Muscle Moments (Knee - Passive Torque)', ...
    'FontSize', 14, 'FontWeight', 'bold');

% ---- 图2: 各肌肉 moment 堆叠图 (膝关节含被动约束力) ----
figure('Color', 'w', 'Position', [100, 100, 1300, 850], ...
    'Name', 'Individual Muscle Moments');

for j = 1 : 3
    subplot(3, 1, j);
    
    detail = momentDetail{j};
    if isempty(detail)
        title(sprintf('%s — 无肌肉 moment 数据', jointTitles{j}));
        continue;
    end
    
    % 对于膝关节 (j=2)，在堆叠中加入被动约束力 (-knee_r.torque)
    if j == 2
        detailWithPassive = [detail, -kneePassiveTorque];
        h = area(time, detailWithPassive);
        cmap = lines(size(detailWithPassive, 2));
        for k = 1 : size(detailWithPassive, 2)
            h(k).FaceColor = cmap(k, :);
            h(k).EdgeColor = 'none';
            h(k).FaceAlpha = 0.7;
        end
    else
        h = area(time, detail);
        cmap = lines(size(detail, 2));
        for k = 1 : size(detail, 2)
            h(k).FaceColor = cmap(k, :);
            h(k).EdgeColor = 'none';
            h(k).FaceAlpha = 0.7;
        end
    end
    
    hold on;
    plot(time, stoMoments(:, j), 'k-', 'LineWidth', 2);
    plot(time, sumMomentsAdjusted(:, j), 'k--', 'LineWidth', 1.5);
    
    xlabel('Time (s)', 'FontSize', 11);
    ylabel('Moment (N\cdotm)', 'FontSize', 11);
    title(sprintf('%s — Individual Muscle Contributions', jointTitles{j}), ...
        'FontSize', 12);
    
    % 图例
    legLabels = momentLabels{j};
    shortNames = cell(size(legLabels));
    for k = 1 : length(legLabels)
        tokens = regexp(legLabels{k}, '^(.+?)\.', 'tokens', 'once');
        if ~isempty(tokens), shortNames{k} = tokens{1}; else shortNames{k} = legLabels{k}; end
    end
    if j == 2
        shortNames{end + 1} = '-PassiveTorque';
    end
    legend([shortNames, {'STO Net', 'Sum (Adj)'}], 'Location', 'bestoutside', 'FontSize', 8, 'Interpreter','none');
    grid on; box on;
end

sgtitle('Individual Muscle Moment Contributions (Stacked, Knee includes Passive Torque)', ...
    'FontSize', 14, 'FontWeight', 'bold');

%% ====== 保存 ======
savePath = fullfile(resultsDir, 'muscle_moment_summation_results.mat');
save(savePath, 'time', 'stoMoments', 'sumMoments', 'sumMomentsAdjusted', ...
    'kneePassiveTorque', 'momentDetail', ...
    'momentLabels', 'targetCoords', 'jointTitles');
fprintf('\n结果已保存至: %s\n', savePath);
fprintf('========== 测试完成 ==========\n');

end
