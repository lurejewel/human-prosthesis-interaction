function S = read_moco_results(varargin)
%READ_MOCO_RESULTS  读取并可视化 MocoInverse 的输出文件。
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
%       S.coordSpeedSource char   速度来源："moco solution /speed" 或
%                                 "ik_radians.mot (central difference)"
%       S.actuatorNames    1×19   cellstr 残差/储备作动器名
%       S.actuatorControl  nT×19  残差/储备控制轨迹
%       S.residuals        table  残差/储备峰值报告（来自 CSV）
%       S.activationReport table  肌肉激活报告（来自 CSV）
%       S.semg             nT×18  处理后的 sEMG 包络（0–1；无实测或落在
%                                 sEMG 处理窗口 170–240 s 之外的为 NaN）
%       S.semgChannels     1×18  对应 sEMG 通道名（无实测的肌肉为空）
%       S.gain             18×1   A-form free gain（来自 moco_inverse_gains_<tag>.csv；
%                                 无 sEMG 的肌肉为 NaN，no-gain run 为 1）
%       S.gainSource       string "free gain"/"free gain (mean over segments)"/"fixed gain=1"
%       S.s                18×1   B-form scale s（来自 moco_inverse_s_used_<tag>.csv）
%       S.scaleKind        "A_gain" / "B_s"
%       S.semgScaled       nT×18  A: S.semg .* S.gain；B: S.semg ./ S.s
%                                 （图 1 的绿色曲线；B 的等值条件是 u=e/s）
%       S.grf             struct  GRF: t, right.vx/vy/vz, left.vx/vy/vz
%       S.summary          string 求解摘要文本
%       S.solution         table  完整解表（列名已清洗为合法变量名）
%       S.labels           cellstr 与 S.solution 列一一对应的原始列名
%       S.ikValue          nT×19  喂给 Moco 的 IK（滤波+弧度化，插值到 S.time）
%       S.ikSpeed          nT×19  其数值微分速度
%       S.ikRawValue       nT×19  原始未滤波 IK（度->弧度，插值到 S.time）
%       S.ikRawSpeed       nT×19  其数值微分速度
%
%   图 2 会把 Moco 的坐标轨迹与 IK 画在一起：实线 = Moco 解，虚线 = 喂给
%   Moco 的滤波后 IK，点线 = 原始未滤波 IK；子图标题给出 Moco 与滤波 IK 的 RMS
%   差异。MocoInverse 的坐标是被样条约束的，所以与滤波 IK 的差异通常极小；与
%   原始 IK 的可见差异主要来自 6 Hz 低通滤波。
%
%   关于「兴奋 vs 激活」：MocoInverse 同时解两者。兴奋 x 是控制变量，激活 a 是
%   状态变量，由一阶激活动力学联系：da/dt = (x - a)/tau_a（并约束 a(0)=x(0)）。
%   完整 solution 文件中两列都在；单独导出的 *_excitations.sto 只含兴奋列。

    % Accept either a struct of options (the documented form) or
    % name-value pairs.  The previous `arguments` block only accepted
    % name-value pairs, so read_moco_results(struct('Tag','x')) failed
    % with "invalid argument at position 1. Name required."
    if nargin == 0
        opt = struct();
    elseif nargin == 1 && isstruct(varargin{1})
        opt = varargin{1};
    else
        if mod(numel(varargin), 2) ~= 0
            error('read_moco_results:args', ...
                  'Options must be a struct or name-value pairs.');
        end
        opt = struct();
        for k = 1:2:numel(varargin)
            name = varargin{k};
            if ~(ischar(name) || isstring(name))
                error('read_moco_results:args', ...
                      'Option names must be character vectors or strings.');
            end
            opt.(char(name)) = varargin{k+1};
        end
    end

    % Defaults + type normalisation.
    if ~isfield(opt, 'OutputDir') || isempty(opt.OutputDir)
        opt.OutputDir = "";
    else
        opt.OutputDir = string(opt.OutputDir);
    end
    if ~isfield(opt, 'Tag') || isempty(opt.Tag)
        opt.Tag = "";
    else
        opt.Tag = string(opt.Tag);
    end
    if ~isfield(opt, 'Plot') || isempty(opt.Plot)
        opt.Plot = true;
    else
        opt.Plot = logical(opt.Plot);
    end
    if ~isfield(opt, 'SaveFigures') || isempty(opt.SaveFigures)
        opt.SaveFigures = false;
    else
        opt.SaveFigures = logical(opt.SaveFigures);
    end
    if ~isfield(opt, 'FigDir') || isempty(opt.FigDir)
        opt.FigDir = "";
    else
        opt.FigDir = string(opt.FigDir);
    end
    if ~isfield(opt, 'RightOnly') || isempty(opt.RightOnly)
        opt.RightOnly = true;
    else
        opt.RightOnly = logical(opt.RightOnly);
    end
    S.rightOnly = opt.RightOnly;

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
    S.files.gains        = fullfile(outDir, sprintf('moco_inverse_gains%s.csv', suf));
    S.files.sUsed        = fullfile(outDir, sprintf('moco_inverse_s_used%s.csv', suf));
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

    % 右腿-only：所有 *_r 肌肉（包括无 sEMG reference 的右腿肌肉）。
    S.musclePlotIdx = find(endsWith(S.muscles, '_r'));
    if isempty(S.musclePlotIdx), S.musclePlotIdx = 1:numel(S.muscles); end
    S.actuatorPlotIdx = find(cellfun(@(s) isempty(regexp(s, '_l($|_)', 'once')), S.actuatorNames));
    if isempty(S.actuatorPlotIdx), S.actuatorPlotIdx = 1:numel(S.actuatorNames); end

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
    [S.ikValue, S.ikSpeed, S.ikTime, ikNames] = loadIKCurves(ikFiltered, S);
    [S.ikRawValue, S.ikRawSpeed, ~, ~] = loadIKCurves(ikRawPath, S);

    % --sagittal-cost 的解只有 18 个 activation state，没有 /jointset/.../value
    % 和 /speed 列。这种解用 ik_radians.mot 的坐标值 + 中心差分速度填充
    % S.coordValue / S.coordSpeed，这样 S.coordSpeed 和 RMS 报告仍然可用。
    if isempty(S.coordNames) && ~isempty(S.ikValue)
        S.coordNames       = ikNames;
        S.coordValue       = S.ikValue;
        S.coordSpeed       = S.ikSpeed;
        S.coordSpeedSource = "ik_radians.mot (central difference)";
        S.coordValueSource = "ik fallback (solution has no Moco coords)";
    else
        S.coordValueSource = "moco solution /value";
        S.coordSpeedSource = "moco solution /speed";
    end
    % IK fallback 后 coordNames 才被填充，必须在这里重算索引。
    S.coordPlotIdx = find(cellfun(@(s) isempty(regexp(s, '_l($|_)', 'once')), S.coordNames));
    if isempty(S.coordPlotIdx), S.coordPlotIdx = 1:numel(S.coordNames); end

    % 处理后的 sEMG 包络（只有 10 块肌肉有实测；其余为 NaN）
    [S.semg, S.semgChannels] = loadSemg(fullfile(thisDir, 'outputs', 'semg'), S);

    % B 形式：moco_inverse_s_used_<tag>.csv（segment,muscle,s）。
    % (e - s u)^2 的等值条件是 u = e/s，所以 control 空间里画 sEMG / s。
    % A 形式：moco_inverse_gains_<tag>.csv，画 sEMG .* g。
    measured = any(~isnan(S.semg), 1)';
    if isfile(S.files.sUsed)
        S.s = loadSValues(S.files.sUsed, S.muscles);
        S.s(~measured) = NaN;
        S.semgScaled = S.semg ./ S.s';
        S.scaleKind = "B_s";
        S.semgPlotLabel = "sEMG / s";
    else
        [S.gain, S.gainSource] = loadGains(S.files.gains, S.muscles);
        if contains(S.gainSource, "fixed")
            S.gain(~measured) = NaN;       % 无 sEMG 的肌肉不标 gain
        end
        S.gain(measured & isnan(S.gain)) = 1.0;
        S.semgScaled = S.semg .* S.gain';
        S.scaleKind = "A_gain";
        if contains(S.gainSource, "free")
            S.semgPlotLabel = "sEMG × free gain";
        else
            S.semgPlotLabel = "sEMG (gain=1)";
        end
    end

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

    % 地面反力（从 summary 的 grf=... 找这次运行实际用的 GRF 文件）
    S.grf = loadGRF(outDir, S);

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

    % Moco 坐标轨迹与 IK 的 RMS 差异。sagittal-cost 解的 coordValue 本身就是
    % 从 IK 来的，这时比较没有意义，跳过。
    if ~isempty(S.ikValue) && ~isempty(S.coordValue) && ...
            S.coordSpeedSource == "moco solution /speed"
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
        f2 = plotKinematics(S);
        fGRF = plotGRF(S);
        f4 = plotActuators(S);
        if opt.SaveFigures
            if opt.FigDir == ""
                figDir = fullfile(outDir, "figures");
            else
                figDir = char(opt.FigDir);
            end
            if ~exist(figDir, 'dir'), mkdir(figDir); end
            figs = {f1, f2, fGRF, f4};
            names = ["01_activation_vs_excitation"; "02_kinematics_value"; ...
                     "03_grf"; "04_residual_reserve"];
            for f = 1:numel(figs)
                if ~isempty(figs{f}) && isgraphics(figs{f})
                    exportgraphics(figs{f}, fullfile(figDir, names(f) + ".png"), ...
                                   'Resolution', 150);
                end
            end
            fprintf('图片已保存到: %s\n', figDir);
        end
    end
end

% =========================================================================
function f = plotMuscles(S)
% 图 1：默认只画右腿 5 块肌肉（可用 RightOnly=false 画全部）
    f = figure('Name', '肌肉激活 vs 兴奋 (activation vs excitation)', 'NumberTitle', 'off');
    idx = S.musclePlotIdx;
    if isfield(S,'rightOnly') && S.rightOnly, idx = S.musclePlotIdx; else, idx = 1:numel(S.muscles); end
    n = numel(idx);
    tl = tiledlayout(f, ceil(n/3), min(3,n), 'TileSpacing', 'compact', 'Padding', 'compact');
    for kk = 1:n
        k = idx(kk);
        nexttile; hold on; grid on;
        plot(S.time, S.activation(:, k), '-',  'LineWidth', 1.4, 'Color', [0.00 0.45 0.74]);
        plot(S.time, S.excitation(:, k), '--', 'LineWidth', 1.0, 'Color', [0.85 0.33 0.10]);
        if isfield(S, 'semgScaled') && ~isempty(S.semgScaled)
            semgPlot = S.semgScaled(:, k);
        elseif isfield(S, 'semg') && ~isempty(S.semg)
            semgPlot = S.semg(:, k);
        else
            semgPlot = [];
        end
        hasEmg = ~isempty(semgPlot) && any(~isnan(semgPlot));
        if hasEmg
            plot(S.time, semgPlot, '-', 'LineWidth', 1.0, ...
                 'Color', [0.15 0.60 0.20]);
        end
        yline(1.0, 'k:');
        yline(0.98, 'r:', 'Alpha', 0.6);
        ylim([0, 1.1]);

        % 每块肌肉标题里标注 scale：B 用 s，A 用 g
        if isfield(S, 'scaleKind') && S.scaleKind == "B_s"
            sv = NaN;
            if isfield(S, 's') && numel(S.s) >= k
                sv = S.s(k);
            end
            if isnan(sv)
                scaleTxt = 's=n/a';
            else
                scaleTxt = sprintf('s=%.2f', sv);
            end
        else
            g = NaN;
            if isfield(S, 'gain') && numel(S.gain) >= k
                g = S.gain(k);
            end
            if isnan(g)
                scaleTxt = 'g=n/a';
            elseif isfield(S, 'gainSource') && contains(S.gainSource, "mean")
                scaleTxt = sprintf('g=%.2f (mean)', g);
            elseif isfield(S, 'gainSource') && contains(S.gainSource, "fixed")
                scaleTxt = sprintf('g=%.2f (fixed)', g);
            else
                scaleTxt = sprintf('g=%.2f', g);
            end
        end
        title(sprintf('%s  (%s)', S.muscles{k}, scaleTxt), 'Interpreter', 'none');
        if kk == 1
            leg = {'activation (state)', 'excitation (control)'};
            if hasEmg
                if isfield(S, 'semgPlotLabel')
                    leg{end+1} = char(S.semgPlotLabel);
                else
                    leg{end+1} = 'sEMG (processed)';
                end
            end
            legend(leg, 'Location', 'southeast', 'FontSize', 7);
        end
        if mod(kk - 1, 3) == 0, ylabel('(-)'); end
        if kk > n - 3, xlabel('t (s)'); end
    end
    linkaxes(tl.Children, 'x');
end

% =========================================================================
function f = plotKinematics(S)
% 图 2：关节角度（deg），Moco 与 IK 对比；
%        实线 = Moco，虚线 = 滤波后 IK（喂给 Moco 的），点线 = 原始未滤波 IK；
%        平动自由度保持 m
    if isempty(S.coordNames)
        f = gobjects(0);
        return;
    end
    isTrans = endsWith(S.coordNames, ["_tx", "_ty", "_tz"]);
    haveIK  = ~isempty(S.ikValue);
    haveRaw = ~isempty(S.ikRawValue);
    fallback = isfield(S,'coordValueSource') && contains(S.coordValueSource, "fallback");
    idx = 1:numel(S.coordNames);
    if isfield(S,'rightOnly') && S.rightOnly, idx = S.coordPlotIdx; end
    n = numel(idx);
    if n == 0
        f = gobjects(0);
        return;
    end
    f2 = figure('Name', '坐标值 Moco vs IK (value)', 'NumberTitle', 'off');
    t2 = tiledlayout(f2, ceil(n/4), min(4,n), 'TileSpacing', 'compact', 'Padding', 'compact');
    for kk = 1:n
        k = idx(kk);
        nexttile; hold on; grid on;
        if isTrans(k)
            if ~fallback, plot(S.time, S.coordValue(:, k), '-', 'LineWidth', 1.2); end
            if haveIK,  plot(S.time, S.ikValue(:, k), '--', 'LineWidth', 1.2); end
            if haveRaw, plot(S.time, S.ikRawValue(:, k), ':', ...
                             'Color', [0.5 0.5 0.5], 'LineWidth', 1.2); end
            ylabel('m');  unit = 'm';
        else
            if ~fallback, plot(S.time, rad2deg(S.coordValue(:, k)), '-', 'LineWidth', 1.2); end
            if haveIK,  plot(S.time, rad2deg(S.ikValue(:, k)), '--' ...
                            , 'LineWidth', 1.2); end
            if haveRaw, plot(S.time, rad2deg(S.ikRawValue(:, k)), ':', ...
                             'Color', [0.5 0.5 0.5], 'LineWidth', 1.2); end
            ylabel('deg');  unit = 'deg';
        end
        if kk == 1
            if fallback, leg = {'IK (filtered)'}; else, leg = {'Moco'}; end
            if haveIK && ~fallback, leg{end+1} = 'IK (filtered)'; end
            if haveRaw, leg{end+1} = 'IK (raw)'; end
            legend(leg, 'Location', 'best', 'FontSize', 7);
        end
        if fallback
            title(S.coordNames{k}, 'Interpreter', 'none');
        elseif haveIK
            d = sqrt(mean((S.coordValue(:, k) - S.ikValue(:, k)).^2));
            if ~isTrans(k), d = rad2deg(d); end
            title(sprintf('%s  (RMS %.3g %s)', S.coordNames{k}, d, unit), ...
                  'Interpreter', 'none');
        else
            title(S.coordNames{k}, 'Interpreter', 'none');
        end
        if kk > n - 4, xlabel('t (s)'); end
    end
    if fallback
        sgtitle('IK reference only: this solution has no Moco coordinate states', ...
                'Interpreter', 'none', 'FontSize', 9);
    end
    linkaxes(t2.Children, 'x');

    f = f2;
end

% =========================================================================
function f = plotActuators(S)
% 图 4：残差/储备作动器使用情况
    if isempty(S.residuals), f = gobjects(0); return; end
    res = S.residuals;
    if isfield(S,'rightOnly') && S.rightOnly
        res = res(S.actuatorPlotIdx, :);
    end
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
function f = plotGRF(S)
% 图 3：地面反力 vx / vy / vz，right 实线，left 虚线
    f = figure('Name', '地面反力 GRF (vx, vy, vz)', 'NumberTitle', 'off');
    tl = tiledlayout(f, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    if ~isfield(S, 'grf') || isempty(S.grf) || isempty(S.grf.t) || ...
            ~isfield(S.grf, 'right') || isempty(S.grf.right.vx)
        nexttile; axis off;
        text(0.5, 0.5, '未找到 GRF 文件', 'HorizontalAlignment', 'center');
        return;
    end
    t = S.grf.t;
    mask = t >= S.time(1) & t <= S.time(end);
    comps = {'vx', 'vy', 'vz'};
    ylabs = {'v_x (N)', 'v_y (N)', 'v_z (N)'};
    for c = 1:3
        nexttile; hold on; grid on;
        plot(t(mask), S.grf.right.(comps{c})(mask), '-', 'LineWidth', 1.3);
        ylabel(ylabs{c});
        if c == 1
            legend({'right'}, 'Location', 'best');
            title('地面反力（GRF，右腿）');
        end
        if c == 3, xlabel('t (s)'); end
    end
    linkaxes(tl.Children, 'x');
end

% =========================================================================
function grf = loadGRF(outDir, S)
%LOADGRF  读取本次运行实际使用的 GRF（summary 里的 grf=...），
%   解析 ExternalLoads XML 的 right/left force_identifier，返回 vx/vy/vz。
    grf = struct('t', [], 'xml', "", 'mot', "");
    grf.right = struct('vx', [], 'vy', [], 'vz', []);
    grf.left  = struct('vx', [], 'vy', [], 'vz', []);

    xml = "";
    if ~isempty(S.summary)
        line = S.summary(startsWith(S.summary, "grf="));
        if ~isempty(line)
            xml = extractAfter(line(1), "grf=");
        end
    end
    % Always prefer the heel-strike-corrected GRF, if it exists, so the figure
    % shows the same (fixed) shear forces that were fed to MocoInverse.
    fixedXml = fullfile(fileparts(mfilename('fullpath')), '..', 'data', ...
                        'SQR_walking', 'level_walking_grf_fixed.xml');
    if isfile(fixedXml), xml = fixedXml; end
    if xml == "" || ~isfile(xml)
        d = dir(fullfile(outDir, "grf_window_*.xml"));
        if isempty(d), return; end
        [~, i] = max([d.datenum]);
        xml = fullfile(outDir, d(i).name);
    end
    if ~isfile(xml), return; end
    grf.xml = string(xml);

    txt = fileread(xml);
    m = regexp(txt, '<datafile>\s*(.*?)\s*</datafile>', 'tokens', 'once');
    if isempty(m), return; end
    mot = m{1};
    if ~isfile(mot)
        mot = fullfile(fileparts(xml), mot);
    end
    if ~isfile(mot), return; end
    grf.mot = string(mot);

    [T, labels] = readStoTable(mot);
    grf.t = T.time;

    blocks = regexp(txt, '<ExternalForce name="([^"]+)">(.*?)</ExternalForce>', 'tokens');
    for b = 1:numel(blocks)
        nm = lower(blocks{b}{1});
        body = blocks{b}{2};
        fi = regexp(body, '<force_identifier>\s*(.*?)\s*</force_identifier>', 'tokens', 'once');
        if isempty(fi), continue; end
        pref = strtrim(fi{1});       % e.g. 'ground_force2_v'
        for c = {'x', 'y', 'z'}
            j = find(strcmp(labels, [pref c{1}]), 1);
            if isempty(j), continue; end
            if contains(nm, 'right')
                grf.right.(['v' c{1}]) = T{:, j};
            elseif contains(nm, 'left')
                grf.left.(['v' c{1}]) = T{:, j};
            end
        end
    end
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
function [vals, channels] = loadSemg(semgDir, S)
%LOADSEMG  读取处理后的 sEMG 包络，并插值到 S.time。
%   文件: <semgDir>/level_walking_semg_processed.mat，变量 t / emgEnvNorm。
%   注意: 该文件只覆盖 sEMG 处理窗口(默认 170~240 s)，窗口外不给外推值(返回 NaN)。
%   通道顺序（1..10）:
%       rect_fem_r, vasti_r, hamstrings_r, tib_ant_r, gastroc_r,
%       rect_fem_l, vasti_l, hamstrings_l, tib_ant_l, gastroc_l
%   没有对应 sEMG 的肌肉保持 NaN。
    vals = nan(numel(S.time), numel(S.muscles));
    channels = cell(1, numel(S.muscles));
    f = fullfile(semgDir, 'level_walking_semg_processed.mat');
    if ~exist(f, 'file'), return; end
    d = load(f, 't', 'emgEnvNorm');
    if ~isfield(d, 't') || ~isfield(d, 'emgEnvNorm'), return; end
    t = d.t(:);
    E = d.emgEnvNorm;
    if size(E, 2) < 10, return; end
    chMap = {'rect_fem_r','vasti_r','hamstrings_r','tib_ant_r','gastroc_r', ...
             'rect_fem_l','vasti_l','hamstrings_l','tib_ant_l','gastroc_l'};
    for k = 1:numel(S.muscles)
        j = find(strcmpi(chMap, S.muscles{k}), 1);
        if isempty(j), continue; end
        vals(:, k) = interp1(t, E(:, j), S.time, 'linear', NaN);   % 窗口外保持 NaN
        channels{k} = chMap{j};
    end
end

% =========================================================================
function s = loadSValues(sPath, muscles)
%LOADSVALUES  读取 moco_inverse_s_used_<tag>.csv 的 per-muscle s。
%   文件列: segment,muscle,s。多段时取平均；缺失肌肉保持 NaN。
    s = nan(numel(muscles), 1);
    if ~isfile(sPath), return; end
    try
        T = readtable(sPath);
    catch
        return;
    end
    vars = lower(string(T.Properties.VariableNames));
    mi = find(vars == "muscle", 1);
    si = find(vars == "s" | vars == "scale", 1);
    if isempty(mi) || isempty(si), return; end
    srcMus = string(T{:, mi}); srcMus = srcMus(:);
    srcS = T{:, si};
    if iscell(srcS), srcS = cell2mat(srcS); end
    for k = 1:numel(muscles)
        rows = srcMus == string(muscles{k});
        if any(rows)
            s(k) = mean(srcS(rows), 'omitnan');
        end
    end
end

% =========================================================================
function [gains, source] = loadGains(gainPath, muscles)
%LOADGAINS  读取 moco_inverse_gains_<tag>.csv 的 per-muscle gain。
%   文件列: segment,muscle,gain。若有多段，取各段平均（图 1 是整段图，
%   分段 gain 的逐点切换不值得引入；平均会在标题里注明）。
%   没有 gains 文件（no-gain run）时返回全 1，source="fixed gain=1"。
    gains = nan(numel(muscles), 1);
    source = "fixed gain=1";
    if ~isfile(gainPath), gains(:) = 1.0; return; end
    try
        T = readtable(gainPath);
    catch
        gains(:) = 1.0; return;
    end
    vars = lower(string(T.Properties.VariableNames));
    mi = find(vars == "muscle", 1);
    gi = find(vars == "gain", 1);
    if isempty(mi) || isempty(gi)
        gains(:) = 1.0; return;
    end
    srcMus  = string(T{:, mi});
    srcMus  = srcMus(:);
    srcGain = T{:, gi};
    if iscell(srcGain), srcGain = cell2mat(srcGain); end
    si = find(vars == "segment", 1);
    nSeg = 1;
    if ~isempty(si), nSeg = numel(unique(T{:, si})); end
    for k = 1:numel(muscles)
        rows = srcMus == string(muscles{k});
        if any(rows)
            gains(k) = mean(srcGain(rows), 'omitnan');
        end
    end
    if nSeg > 1
        source = "free gain (mean over segments)";
    else
        source = "free gain";
    end
end

% =========================================================================
function [vals, spds, tRef, names] = loadIKCurves(path, S)
%LOADIKCURVES  读取 IK 坐标文件并按 S.coordNames 插值到 S.time。
%   自动识别 inDegrees 头（度->弧度只作用于转动坐标，平动 _tx/_ty/_tz 不动）。
%   如果 S.coordNames 为空（--sagittal-cost 的解只有激活状态），则使用 IK
%   文件里的全部坐标列，并把 names 返回给调用者。
%   速度用中心差分（gradient）计算，比原来的前向 diff 更准。
    vals = []; spds = []; tRef = []; names = {};
    if ~exist(path, 'file'), return; end

    [ikT, ikLabels] = readStoTable(path);
    tRef = ikT.time;
    degFlag = isInDegrees(path);
    if isempty(S.coordNames)
        names = ikLabels(~strcmpi(ikLabels, 'time'));
    else
        names = S.coordNames;
    end
    isT = endsWith(names, ["_tx", "_ty", "_tz"]);

    vals = nan(numel(S.time), numel(names));
    for c = 1:numel(names)
        j = find(strcmp(ikLabels, names{c}), 1);
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
    dt = median(diff(S.time));
    for c = 1:size(vals, 2)
        spds(:, c) = gradient(vals(:, c), dt);
    end
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
