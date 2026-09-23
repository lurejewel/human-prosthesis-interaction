function S = read_moco_results(opt)
%READ_MOCO_RESULTS  读取并可视化 MocoInverse 的输出文件（纯 MATLAB，无需 OpenSim 接口）。
%
%   用法
%       S = read_moco_results();                              % 读取默认 outputs/ 并绘图
%       S = read_moco_results(struct('Plot', false));         % 只读数据不绘图
%       S = read_moco_results(struct('Tag', 'millard'));      % 读取带 _millard 后缀的一组结果
%       S = read_moco_results(struct('SaveFigures', true, ...
%                                    'FigDir', 'outputs/figures'));   % 保存 PNG
%
%   读取的文件（默认目录 = 本文件同级的 outputs/，即 experiments/moco/outputs/）
%       moco_inverse_solution[_tag].sto      完整解：状态 + 控制（一行 94 列）
%           - 状态(56)：18 个 /forceset/<肌肉>/activation
%                     + 19 个 /jointset/.../value  + 19 个 /speed
%           - 控制(37)：18 个 /forceset/<肌肉>        = 肌肉兴奋 excitation
%                     +  6 个 /forceset/residual_*  = 骨盆残差作动器
%                     + 13 个 /forceset/reserve_*   = 储备作动器
%       moco_inverse_activations[_tag].sto    仅肌肉激活（state，方便与 EMG 对比）
%       moco_inverse_excitations[_tag].sto    仅肌肉兴奋（control）
%       moco_inverse_residuals[_tag].csv      残差/储备作动器峰值报告
%       moco_inverse_activation_report[_tag].csv  每块肌肉峰值/平均激活
%       moco_inverse_summary[_tag].txt        求解摘要
%
%   返回结构体 S 的主要字段
%       S.time             nT×1   时间 (s)
%       S.muscles          1×18   cellstr 肌肉名
%       S.activation       nT×18  肌肉激活 a（状态；da/dt = (x-a)/tau_a）
%       S.excitation       nT×18  肌肉兴奋 x（控制，与 S.muscles 顺序一致）
%       S.coordNames       1×19   cellstr 坐标名
%       S.coordValue       nT×19  坐标值（转动为 rad，平动为 m）
%       S.coordSpeed       nT×19  坐标速度（rad/s 或 m/s）
%       S.actuatorNames    1×19   cellstr 残差/储备作动器名
%       S.actuatorControl  nT×19  残差/储备控制轨迹
%       S.residuals        table  残差/储备峰值报告（来自 CSV）
%       S.activationReport table  肌肉激活报告（来自 CSV）
%       S.summary          string 求解摘要文本
%       S.solution         table  完整解表（列名已清洗为合法变量名）
%       S.labels           cellstr 与 S.solution 列一一对应的原始列名
%       S.ikValue          nT×19  喂给 Moco 的 IK（滤波+弧度化，插值到 S.time）
%       S.ikSpeed          nT×19  其数值微分速度
%       S.ikRawValue       nT×19  原始未滤波 IK（度->弧度，插值到 S.time）
%       S.ikRawSpeed       nT×19  其数值微分速度
%
%   图 2/图 3 会把 Moco 的坐标轨迹与 IK 画在一起：实线 = Moco 解，虚线 = 喂给
%   Moco 的滤波后 IK，点线 = 原始未滤波 IK；子图标题给出 Moco 与滤波 IK 的 RMS
%   差异。MocoInverse 的坐标是被样条约束的，所以与滤波 IK 的差异通常极小；与
%   原始 IK 的可见差异主要来自 6 Hz 低通滤波。
%
%   关于「兴奋 vs 激活」：MocoInverse 同时解两者。兴奋 x 是控制变量，激活 a 是
%   状态变量，由一阶激活动力学联系：da/dt = (x - a)/tau_a（并约束 a(0)=x(0)）。
%   完整 solution 文件中两列都在；单独导出的 *_excitations.sto 只含兴奋列。

    arguments
        opt.OutputDir (1,1) string = ""
        opt.Tag (1,1) string = ""
        opt.Plot (1,1) logical = true
        opt.SaveFigures (1,1) logical = false
        opt.FigDir (1,1) string = ""
    end

    % ---------------------------------------------------------------------
    % 路径
    % ---------------------------------------------------------------------
    thisDir = fileparts(mfilename('fullpath'));
    if opt.OutputDir == ""
        outDir = fullfile(thisDir, "outputs");
    else
        outDir = char(opt.OutputDir);
    end
    if opt.Tag == ""
        suf = "";
    else
        suf = "_" + opt.Tag;
    end

    S.files.solution  = fullfile(outDir, sprintf('moco_inverse_solution%s.sto',  suf));
    S.files.activations  = fullfile(outDir, sprintf('moco_inverse_activations%s.sto', suf));
    S.files.excitations  = fullfile(outDir, sprintf('moco_inverse_excitations%s.sto', suf));
    S.files.residuals    = fullfile(outDir, sprintf('moco_inverse_residuals%s.csv',  suf));
    S.files.activationReport = fullfile(outDir, sprintf('moco_inverse_activation_report%s.csv', suf));
    S.files.summary      = fullfile(outDir, sprintf('moco_inverse_summary%s.txt', suf));

    if ~exist(S.files.solution, 'file')
        error('read_moco_results:noSolution', ...
              ['未找到解文件: %s\n请先在 experiments/moco/ 下运行 ' ...
               'moco_inverse.py 生成输出。'], S.files.solution);
    end

    % ---------------------------------------------------------------------
    % 读取 .sto 文件
    % ---------------------------------------------------------------------
    [S.solution, S.labels] = readStoTable(S.files.solution);
    S.time = S.solution.time;

    if exist(S.files.activations, 'file')
        S.activationsTable = readStoTable(S.files.activations);
    end
    if exist(S.files.excitations, 'file')
        S.excitationsTable = readStoTable(S.files.excitations);
    end

    % ---------------------------------------------------------------------
    % 按列名分类，拆出矩阵字段
    % ---------------------------------------------------------------------
    data  = table2array(S.solution);          % 第 1 列是 time
    labels = S.labels;

    actIdx = [];  actNames = {};
    excIdx = [];  excNames = {};
    valIdx = [];  spdIdx  = [];  coordNames = {};
    actuIdx = []; actuNames = {};

    for k = 2:numel(labels)
        parts = strsplit(labels{k}, '/');
        switch numel(parts)
            case 4   % /forceset/<肌肉>/activation
                if strcmp(parts{4}, 'activation')
                    actIdx(end+1) = k; %#ok<AGROW>
                    actNames{end+1} = parts{3}; %#ok<AGROW>
                end
            case 5   % /jointset/<joint>/<coord>/<value|speed>
                if strcmp(parts{5}, 'value')
                    valIdx(end+1) = k; %#ok<AGROW>
                    coordNames{end+1} = parts{4}; %#ok<AGROW>
                elseif strcmp(parts{5}, 'speed')
                    spdIdx(end+1) = k; %#ok<AGROW>
                end
            case 3   % /forceset/<name>：兴奋或残差/储备作动器
                if startsWith(parts{3}, 'residual_') || startsWith(parts{3}, 'reserve_')
                    actuIdx(end+1) = k; %#ok<AGROW>
                    actuNames{end+1} = parts{3}; %#ok<AGROW>
                else
                    excIdx(end+1) = k; %#ok<AGROW>
                    excNames{end+1} = parts{3}; %#ok<AGROW>
                end
        end
    end

    S.muscles   = actNames;
    S.activation = data(:, actIdx);
    S.coordNames = coordNames;
    S.coordValue = data(:, valIdx);
    S.coordSpeed = data(:, spdIdx);
    S.actuatorNames = actuNames;
    S.actuatorControl = data(:, actuIdx);

    % 兴奋列按 S.muscles 的顺序重排（保证与 S.activation 一一对应）
    S.excitation = nan(size(S.activation));
    for m = 1:numel(S.muscles)
        j = find(strcmp(excNames, S.muscles{m}), 1);
        if ~isempty(j)
            S.excitation(:, m) = data(:, excIdx(j));
        end
    end
    if any(isnan(S.excitation(:)))
        warning('read_moco_results:missingExcitation', ...
                '部分肌肉的兴奋列缺失（旧版输出？），对应位置为 NaN。');
    end

    % ---------------------------------------------------------------------
    % IK 参考运动学（图 2/3 对比用）：优先用喂给 Moco 的滤波+弧度化文件
    % ---------------------------------------------------------------------
    ikFiltered = fullfile(outDir, 'ik_radians.mot');
    ikRawPath  = fullfile(thisDir, '..', 'data', 'SQR_walking', 'level_walking_ik.mot');
    S.ikValue = [];  S.ikSpeed = [];  S.ikTime = [];
    S.ikRawValue = [];  S.ikRawSpeed = [];
    [S.ikValue, S.ikSpeed, S.ikTime] = loadIKCurves(ikFiltered, S);
    [S.ikRawValue, S.ikRawSpeed]     = loadIKCurves(ikRawPath, S);

    % ---------------------------------------------------------------------
    % 读取 CSV 与摘要
    % ---------------------------------------------------------------------
    if exist(S.files.residuals, 'file')
        S.residuals = readtable(S.files.residuals);
    else
        S.residuals = [];
    end
    if exist(S.files.activationReport, 'file')
        S.activationReport = readtable(S.files.activationReport);
    else
        S.activationReport = [];
    end
    if exist(S.files.summary, 'file')
        S.summary = readlines(S.files.summary);
        S.summary = S.summary(strlength(S.summary) > 0);
    else
        S.summary = string.empty;
    end

    % ---------------------------------------------------------------------
    % 控制台报告
    % ---------------------------------------------------------------------
    fprintf('\n==== MocoInverse 结果报告 ====\n');
    fprintf('输出目录 : %s\n', outDir);
    fprintf('时间范围 : %.3f - %.3f s（%d 个时间点）\n', S.time(1), S.time(end), numel(S.time));
    if ~isempty(S.summary)
        fprintf('%s\n', S.summary);
    end
    fprintf('肌肉数   : %d\n', numel(S.muscles));
    if ~isempty(S.activationReport)
        nSat = sum(S.activationReport.peak_activation > 0.98);
        fprintf('激活饱和的肌肉 : %d / %d（峰值 > 0.98）\n', nSat, height(S.activationReport));
        if nSat > 0
            sat = S.activationReport(S.activationReport.peak_activation > 0.98, :);
            for r = 1:height(sat)
                fprintf('    %-12s  peak = %.3f\n', sat.muscle{r}, sat.peak_activation(r));
            end
        end
    end
    if ~isempty(S.residuals)
        isRes = startsWith(S.residuals.actuator, 'residual_');
        nNear = sum(abs(S.residuals.peak_control(isRes)) > 0.9);
        fprintf('残差作动器接近饱和(>0.9) : %d / %d\n', nNear, sum(isRes));
        [~, io] = maxk(S.residuals.peak_force, 3);
        fprintf('出力最大的 3 个作动器:\n');
        for r = 1:numel(io)
            fprintf('    %-45s  peak force = %8.2f\n', S.residuals.actuator{io(r)}, ...
                    S.residuals.peak_force(io(r)));
        end
    end

    % Moco 坐标轨迹与 IK 的 RMS 差异
    if ~isempty(S.ikValue) && ~isempty(S.coordValue)
        isT = endsWith(S.coordNames, ["_tx", "_ty", "_tz"]);
        fprintf('\nMoco 运动学 vs IK（滤波后，即喂给 Moco 的数据）RMS 差异:\n');
        for c = 1:numel(S.coordNames)
            d = sqrt(mean((S.coordValue(:, c) - S.ikValue(:, c)).^2));
            if isT(c)
                fprintf('    %-24s %9.4g m\n', S.coordNames{c}, d);
            else
                fprintf('    %-24s %9.4g deg\n', S.coordNames{c}, rad2deg(d));
            end
        end
        if ~isempty(S.ikRawValue)
            dRaw = sqrt(mean((S.coordValue - S.ikRawValue).^2, 1));
            fprintf('vs 原始(未滤波) IK 的最大 RMS 差异: %.4g deg（主要是滤波所致）\n', ...
                    max(dRaw(~isT)) * 180 / pi);
        end
    end

    % ---------------------------------------------------------------------
    % 绘图
    % ---------------------------------------------------------------------
    if opt.Plot
        f1 = plotMuscles(S);
        f23 = plotKinematics(S);
        f4 = plotActuators(S);
        if opt.SaveFigures
            if opt.FigDir == ""
                figDir = fullfile(outDir, "figures");
            else
                figDir = char(opt.FigDir);
            end
            if ~exist(figDir, 'dir'), mkdir(figDir); end
            figs = [f1; f23(:); f4];
            names = ["01_activation_vs_excitation"; "02_kinematics_value"; ...
                     "03_kinematics_speed"; "04_residual_reserve"];
            for f = 1:numel(figs)
                if ~isempty(figs(f)) && isgraphics(figs(f))
                    exportgraphics(figs(f), fullfile(figDir, names(f) + ".png"), ...
                                   'Resolution', 150);
                end
            end
            fprintf('图片已保存到: %s\n', figDir);
        end
    end
end

% =========================================================================
function f = plotMuscles(S)
% 图 1：18 块肌肉的激活（实线，状态）与兴奋（虚线，控制）
    f = figure('Name', '肌肉激活 vs 兴奋 (activation vs excitation)', 'NumberTitle', 'off');
    tl = tiledlayout(f, 6, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:numel(S.muscles)
        nexttile; hold on; grid on;
        plot(S.time, S.activation(:, k), '-',  'LineWidth', 1.4, 'Color', [0.00 0.45 0.74]);
        plot(S.time, S.excitation(:, k), '--', 'LineWidth', 1.0, 'Color', [0.85 0.33 0.10]);
        yline(1.0, 'k:');
        yline(0.98, 'r:', 'Alpha', 0.6);
        ylim([0, 1.1]);
        title(S.muscles{k}, 'Interpreter', 'none');
        if k == 1
            legend({'activation (state)', 'excitation (control)'}, ...
                   'Location', 'southeast', 'FontSize', 7);
        end
        if mod(k - 1, 3) == 0, ylabel('(-)'); end
        if k > 15, xlabel('t (s)'); end
    end
    linkaxes(tl.Children, 'x');
end

% =========================================================================
function f = plotKinematics(S)
% 图 2/3：关节角度（deg）与角速度（deg/s），Moco 与 IK 对比；
%        实线 = Moco，虚线 = 滤波后 IK（喂给 Moco 的），点线 = 原始未滤波 IK；
%        平动自由度保持 m 与 m/s
    isTrans = endsWith(S.coordNames, ["_tx", "_ty", "_tz"]);
    haveIK  = ~isempty(S.ikValue);
    haveRaw = ~isempty(S.ikRawValue);

    f2 = figure('Name', '坐标值 Moco vs IK (value)', 'NumberTitle', 'off');
    t2 = tiledlayout(f2, 5, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:numel(S.coordNames)
        nexttile; hold on; grid on;
        if isTrans(k)
            plot(S.time, S.coordValue(:, k), '-', 'LineWidth', 1.2);
            if haveIK,  plot(S.time, S.ikValue(:, k), '--', 'LineWidth', 1.2); end
            if haveRaw, plot(S.time, S.ikRawValue(:, k), ':', ...
                             'Color', [0.5 0.5 0.5], 'LineWidth', 1.2); end
            ylabel('m');  unit = 'm';
        else
            plot(S.time, rad2deg(S.coordValue(:, k)), '-', 'LineWidth', 1.2);
            if haveIK,  plot(S.time, rad2deg(S.ikValue(:, k)), '--' ...
                            , 'LineWidth', 1.2); end
            if haveRaw, plot(S.time, rad2deg(S.ikRawValue(:, k)), ':', ...
                             'Color', [0.5 0.5 0.5], 'LineWidth', 1.2); end
            ylabel('deg');  unit = 'deg';
        end
        if k == 1
            leg = {'Moco'};
            if haveIK,  leg{end+1} = 'IK (filtered)'; end
            if haveRaw, leg{end+1} = 'IK (raw)'; end
            legend(leg, 'Location', 'best', 'FontSize', 7);
        end
        if haveIK
            d = sqrt(mean((S.coordValue(:, k) - S.ikValue(:, k)).^2));
            if ~isTrans(k), d = rad2deg(d); end
            title(sprintf('%s  (RMS %.3g %s)', S.coordNames{k}, d, unit), ...
                  'Interpreter', 'none');
        else
            title(S.coordNames{k}, 'Interpreter', 'none');
        end
        if k > 15, xlabel('t (s)'); end
    end
    linkaxes(t2.Children, 'x');

    f3 = figure('Name', '坐标速度 Moco vs IK (speed)', 'NumberTitle', 'off');
    t3 = tiledlayout(f3, 5, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:numel(S.coordNames)
        nexttile; hold on; grid on;
        if isTrans(k)
            plot(S.time, S.coordSpeed(:, k), '-', 'LineWidth', 1.2);
            if haveIK,  plot(S.time, S.ikSpeed(:, k), '--', 'LineWidth', 1.2); end
            if haveRaw, plot(S.time, S.ikRawSpeed(:, k), ':', ...
                             'Color', [0.5 0.5 0.5], 'LineWidth', 1.2); end
            ylabel('m/s');
        else
            plot(S.time, rad2deg(S.coordSpeed(:, k)), '-', 'LineWidth', 1.2);
            if haveIK,  plot(S.time, rad2deg(S.ikSpeed(:, k)), '--', ...
                    'LineWidth', 1.2); end
            if haveRaw, plot(S.time, rad2deg(S.ikRawSpeed(:, k)), ':', ...
                             'Color', [0.5 0.5 0.5], 'LineWidth', 1.2); end
            ylabel('deg/s');
        end
        title(S.coordNames{k}, 'Interpreter', 'none');
        if k > 15, xlabel('t (s)'); end
    end
    linkaxes(t3.Children, 'x');
    f = [f2; f3];
end

% =========================================================================
function f = plotActuators(S)
% 图 4：残差/储备作动器使用情况
    if isempty(S.residuals), f = gobjects(0); return; end
    res = S.residuals;
    isRes = startsWith(res.actuator, 'residual_');
    short = regexprep(res.actuator, '^(residual|reserve)_jointset_', '');

    f = figure('Name', '残差/储备作动器', 'NumberTitle', 'off');
    t4 = tiledlayout(f, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile; hold on; grid on;
    bar(categorical(short(isRes)), res.peak_control(isRes), 'FaceColor', [0.30 0.75 0.93]);
    yline(1.0, 'r--', '控制边界 \pm1', 'LabelVerticalAlignment', 'bottom');
    yline(-1.0, 'r--');
    title('残差作动器：峰值控制（边界 ±1，接近 1 表示容量不足/模型不匹配）');
    ylabel('peak |control|');
    xtickangle(35);

    nexttile; hold on; grid on;
    bar(categorical(short(~isRes)), res.peak_force(~isRes), 'FaceColor', [0.93 0.69 0.13]);
    title('储备作动器：峰值出力（N·m 或 N；数值大说明肌肉力量不足）');
    ylabel('peak force');
    xtickangle(35);
end

% =========================================================================
function [T, labels] = readStoTable(path)
%READSTOTABLE  读取 OpenSim 的 .sto/.mot 文本表（不依赖 OpenSim 接口）。
%   返回 table T（变量名已清洗为合法标识符）与原始列名 cellstr labels。
    fid = fopen(path, 'r');
    if fid < 0
        error('read_moco_results:cannotOpen', '无法打开文件: %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));

    lines = {};
    while ~feof(fid)
        l = fgetl(fid);
        if ~ischar(l), break; end
        lines{end+1} = l; %#ok<AGROW>
    end

    iEnd = find(strcmpi(strtrim(lines), 'endheader'), 1);
    if isempty(iEnd)
        error('read_moco_results:badFile', '文件缺少 endheader: %s', path);
    end
    iLab = [];
    for k = iEnd + 1 : numel(lines)
        if ~isempty(strtrim(lines{k}))
            iLab = k;
            break;
        end
    end
    if isempty(iLab)
        error('read_moco_results:badFile', '文件缺少列名行: %s', path);
    end
    labels = strsplit(strtrim(lines{iLab}));

    nCol = numel(labels);
    rows = zeros(numel(lines) - iLab, nCol);
    nRow = 0;
    for k = iLab + 1 : numel(lines)
        if isempty(strtrim(lines{k})), continue; end
        vals = str2double(strsplit(strtrim(lines{k})));
        if numel(vals) ~= nCol
            warning('read_moco_results:badRow', ...
                    '第 %d 行有 %d 列，期望 %d 列，已跳过。', k, numel(vals), nCol);
            continue;
        end
        nRow = nRow + 1;
        rows(nRow, :) = vals;
    end
    rows = rows(1:nRow, :);

    vn = cellfun(@(s) sanitizeVarName(s), labels, 'UniformOutput', false);
    T = array2table(rows, 'VariableNames', vn);
end

% =========================================================================
function vn = sanitizeVarName(s)
% 把 '/forceset/soleus_r/activation' 这类路径清洗成合法变量名
    vn = regexprep(s, '[^A-Za-z0-9_]', '_');
    vn = regexprep(vn, '^_+', '');
    if isempty(vn) || ~isletter(vn(1))
        vn = ['x_' vn]; %#ok<AGROW>
    end
end

% =========================================================================
function [vals, spds, tRef] = loadIKCurves(path, S)
%LOADIKCURVES  读取 IK 坐标文件并按 S.coordNames 插值到 S.time。
%   自动识别 inDegrees 头（度->弧度只作用于转动坐标，平动 _tx/_ty/_tz 不动）。
    vals = []; spds = []; tRef = [];
    if ~exist(path, 'file'), return; end

    [ikT, ikLabels] = readStoTable(path);
    tRef = ikT.time;
    degFlag = isInDegrees(path);
    isT = endsWith(S.coordNames, ["_tx", "_ty", "_tz"]);

    vals = nan(numel(S.time), numel(S.coordNames));
    for c = 1:numel(S.coordNames)
        j = find(strcmp(ikLabels, S.coordNames{c}), 1);
        if isempty(j), continue; end
        v = ikT{:, j};
        if degFlag && ~isT(c)
            v = deg2rad(v);
        end
        if ~issorted(tRef)
            [tRef, ord] = sort(tRef);
            v = v(ord);
        end
        vals(:, c) = interp1(tRef, v, S.time, 'linear', 'extrap');
    end

    spds = nan(size(vals));
    dt = diff(S.time);
    for c = 1:size(vals, 2)
        spds(2:end, c) = diff(vals(:, c)) ./ dt;
    end
    spds(1, :) = spds(2, :);
end

% =========================================================================
function tf = isInDegrees(path)
%ISINDEGREES  检查 .sto/.mot 头部是否有 inDegrees=yes（不依赖 OpenSim 接口）
    tf = false;
    fid = fopen(path, 'r');
    if fid < 0, return; end
    cleaner = onCleanup(@() fclose(fid));
    n = 0;
    while ~feof(fid) && n < 30
        l = fgetl(fid);
        if ~ischar(l), break; end
        n = n + 1;
        tok = regexpi(strtrim(l), '^inDegrees\s*=\s*(yes|no)$', 'tokens', 'once');
        if ~isempty(tok)
            tf = strcmpi(tok{1}, 'yes');
            return;
        end
        if strcmpi(strtrim(l), 'endheader'), return; end
    end
end
