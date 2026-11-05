function T = read_sto_mot_table(filePath)
%READ_STO_MOT_TABLE 读取 OpenSim / SCONE 的 .sto/.mot 文本表（不依赖 OpenSim API）。
%
% 这个函数存在的理由
% ------------------
% 本仓库的 MocoInverse 输出（.sto）、IK 输出（.mot）、GRF（.mot）都只是
% “带一段文本头的、制表符分隔的数值表”。如果为了读一个文件就启动
% org.opensim.modeling.TimeSeriesTable，会有两个问题：
%   1) 某些机器上 OpenSim-MATLAB 接口未配置好，连 QC 都做不了；
%   2) 我们只需要 time + 列名 + 数值矩阵，没必要创建 Java 对象。
% 因此这里实现一个纯 MATLAB 读取器，供接口、QC、相位检测共用。
%
% 支持的两种头格式
% ----------------
% (a) 标准 OpenSim .sto/.mot：
%         header 行（key=value，可能有多行）
%         endheader
%         time<TAB>col1<TAB>col2...
%         数值行...
%
% (b) SCONE 的 .sto：
%         头行每一行后面都被制表符 padding 到整行宽度；
%         例如 "nRows=2001\t\t\t...".
%         这些 padding 对 strtrim() 无害，因此本函数可以直接读 UN.sto。
%
% 重要约定
% --------
% * 第一列始终被当作时间列，列名通常为 "time"；
% * 本函数不做 deg/rad 转换，也不做单位换算。调用者若需要角度转换，
%   应检查 T.inDegrees 并只对转动坐标做 deg2rad（见
%   moco_inverse_to_lasso_dataset.m 的 IK 分支）；
% * 本函数不依赖 Signal Processing / Statistics 工具箱。
%
% 输入
% ----
%   filePath : char 或 string，.sto/.mot 文件路径。
%
% 输出结构 T
% ----------
%   T.time        nT×1  double，第一列时间；
%   T.labels      1×nC  cellstr，包含 "time" 在内的全部列名；
%   T.data        nT×nC double，包含时间列在内的全部数值；
%   T.nRows       实际读到的数据行数；
%   T.nColumns    列数（包含时间列）；
%   T.inDegrees   逻辑值，头部是否写了 inDegrees=yes；
%   T.header      cellstr，endheader 之前的所有头行（去掉尾部空白）；
%   T.file        char，绝对路径。
%
% 失败模式
% --------
% * 文件不存在 / 打不开；
% * 找不到 endheader；
% * endheader 之后找不到列名行；
% * 找不到任何数值行。
% 这些都会给出带文件名和原因的错误，而不是静默返回空表。

% -------------------------------------------------------------------------
% 0. 输入校验
% -------------------------------------------------------------------------
if nargin < 1 || isempty(filePath)
    error('read_sto_mot_table:badInput', ...
        'filePath must be a non-empty char/string.');
end
filePath = char(filePath);
if ~isfile(filePath)
    error('read_sto_mot_table:notFound', 'File not found: %s', filePath);
end
T.file = filePath;

fid = fopen(filePath, 'r');
if fid < 0
    error('read_sto_mot_table:cannotOpen', 'Cannot open file: %s', filePath);
end
% onCleanup 保证即使后面 error() 提前退出，文件句柄也会被关闭。
% 这对批处理/循环读很多文件很重要，否则 Windows 上会留下锁。
cleaner = onCleanup(@() fclose(fid));

% -------------------------------------------------------------------------
% 1. 扫描头部，直到 endheader，并解析 inDegrees
% -------------------------------------------------------------------------
% 注意：
% * 不能用 fscanf/textscan 直接跳过头，因为头部行数不固定；
% * SCONE 的头行可能带大量尾部 \t，strtrim 后再比较最稳；
% * inDegrees 只记录，不在此处转换。
header = {};
inDegrees = false;
endHeaderLine = '';
while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line)
        % fgetl 到 EOF 返回 -1（double），不是 char。
        break;
    end
    s = strtrim(line);
    if strcmpi(s, 'endheader')
        endHeaderLine = line;
        break;
    end
    header{end+1} = s; %#ok<AGROW>
    tok = regexpi(s, '^inDegrees\s*=\s*(yes|no)$', 'tokens', 'once');
    if ~isempty(tok)
        inDegrees = strcmpi(tok{1}, 'yes');
    end
end
if isempty(endHeaderLine)
    error('read_sto_mot_table:noEndHeader', ...
        'No "endheader" line found in %s.', filePath);
end

% -------------------------------------------------------------------------
% 2. 读取列名行（endheader 之后的第一个非空行）
% -------------------------------------------------------------------------
% 有些写出的文件在 endheader 和列名之间会多一个空行，所以这里循环跳过空行。
% 列名用空白（空格或制表符）分隔；OpenSim 列名本身不允许空格。
labelLine = '';
while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end
    if ~isempty(strtrim(line))
        labelLine = line;
        break;
    end
end
if isempty(labelLine)
    error('read_sto_mot_table:noLabels', ...
        'No column-label line found after endheader in %s.', filePath);
end

labels = regexp(strtrim(labelLine), '\s+', 'split');
labels = labels(~cellfun(@isempty, labels));
nC = numel(labels);
if nC < 2
    error('read_sto_mot_table:badLabels', ...
        'Expected at least 2 columns in %s, found %d.', filePath, nC);
end

% -------------------------------------------------------------------------
% 3. 读取数值矩阵
% -------------------------------------------------------------------------
% textscan 比逐行 str2double 快很多，尤其对 70k 行的 GRF/IK 文件。
% 'CollectOutput', true 把 nC 个 %f 字段收集成一个 nT×nC 矩阵。
% 'CommentStyle','#' 是防御性的：OpenSim 某些文件可能在数据区加注释。
%
% 注意：如果某行数值个数不足 nC，textscan 会提前停止，data 行数会少于
% 文件末尾的实际行数。Moco/IK/GRF 的正常输出不会出现这种情况；若出现，
% 调用者可通过 T.nRows 与头部 nRows 对比来发现。
fmt = repmat('%f', 1, nC);
C = textscan(fid, fmt, 'CollectOutput', true, 'CommentStyle', '#');
data = C{1};
if isempty(data)
    error('read_sto_mot_table:noData', 'No numeric rows found in %s.', filePath);
end

T.header     = header;
T.inDegrees  = inDegrees;
T.labels     = labels;
T.data       = data;
T.time       = data(:, 1);
T.nRows      = size(data, 1);
T.nColumns   = nC;

% 防御性检查：如果 textscan 实际读到的列数和列名数不一致，给出警告。
% 正常文件不会触发；手工生成的 .sto 如果少写/多写一列，这里会提醒。
if size(data, 2) ~= nC
    warning('read_sto_mot_table:columnMismatch', ...
        'Read %d columns but header has %d labels in %s.', ...
        size(data, 2), nC, filePath);
end
end
