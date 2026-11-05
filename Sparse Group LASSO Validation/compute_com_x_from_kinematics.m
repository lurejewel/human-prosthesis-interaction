function C = compute_com_x_from_kinematics(modelFile, time, q, qdot, qNames)
%COMPUTE_COM_X_FROM_KINEMATICS 计算 pelvis / calcn 质心的全局 X 坐标。
%
% 为什么需要这个函数
% ------------------
% 旧 LASSO 输入里的 `leg1_r.state` 来自 SCONE 的步态状态机。MocoInverse
% 没有这个列，所以接口必须自己重建 0..4 相位。重建规则必须和预测仿真
% 里的 `functions/core/cal_gait_phase.m` 完全一致，否则 LASSO 的相位标签
% 和运行时控制器看到的相位不同。
%
% `cal_gait_phase.m` 的规则不仅用 GRF，还用：
%   * pelvis 质心的全局 X 坐标；
%   * calcn_r / calcn_l 质心的全局 X 坐标。
% 本函数就是计算这三个 X 坐标。
%
% 计算方式与 cal_gait_phase.m 相同：
%
%   pelvisCOM  = pelvis.findStationLocationInGround(state, pelvis.getMassCenter())
%   calcnRCOM  = calcn_r.findStationLocationInGround(state, calcn_r.getMassCenter())
%   calcnLCOM  = calcn_l.findStationLocationInGround(state, calcn_l.getMassCenter())
%
% 只取返回 Vec3 的第 0 个分量（X），因为 cal_gait_phase.m 只用 X 判断
% “脚在骨盆前方还是后方”。
%
% 输入
% ----
%   modelFile : 预测仿真模型 .osim（默认 human0918）
%   time      : nT×1 [s]
%   q, qdot   : nT×nQ 坐标值/速度
%   qNames    : 1×nQ 与 q/qdot 列对应的坐标名
%
% 输出 C
% ------
%   C.time      nT×1
%   C.pelvis    nT×1  pelvis 质心全局 X [m]
%   C.calcnR    nT×1  calcn_r 质心全局 X [m]
%   C.calcnL    nT×1  calcn_l 质心全局 X [m]
%   C.modelFile char
%
% 注意
% ----
% * 这里只需要位置，不需要动力学；realizePosition 足够。
% * 只设置模型拥有的坐标；IK 里多出来的 3D 坐标（pelvis_list 等）对
%   human0918 不存在，直接忽略。
% * q/qdot 必须有限；上层接口已经检查过。

if ~isfile(modelFile)
    error('compute_com_x_from_kinematics:noModel', ...
        'Model file not found: %s', modelFile);
end
time = time(:);
nT = numel(time);
if size(q, 1) ~= nT || size(qdot, 1) ~= nT
    error('compute_com_x_from_kinematics:sizeMismatch', ...
        'q/qdot must have the same number of rows as time.');
end
if size(q, 2) ~= numel(qNames)
    error('compute_com_x_from_kinematics:coordMismatch', ...
        'q has %d columns but qNames has %d entries.', ...
        size(q, 2), numel(qNames));
end

model = org.opensim.modeling.Model(char(modelFile));
state = model.initSystem();

% ---- 坐标 handle 与 q 列映射 -------------------------------------------
cs = model.getCoordinateSet();
nCoord = cs.getSize();
coordHandles = cell(1, nCoord);
coordCol = zeros(1, nCoord);
for i = 1:nCoord
    c = cs.get(i - 1);
    coordHandles{i} = c;
    j = find(strcmp(qNames, char(c.getName())), 1);
    if ~isempty(j)
        coordCol(i) = j;
    end
end

% ---- 三个 body 的质心局部坐标 -----------------------------------------
bodyPelvis = model.getBodySet().get('pelvis');
bodyCalcnR = model.getBodySet().get('calcn_r');
bodyCalcnL = model.getBodySet().get('calcn_l');
pelvisCOM_local = bodyPelvis.getMassCenter();
calcnRCOM_local = bodyCalcnR.getMassCenter();
calcnLCOM_local = bodyCalcnL.getMassCenter();

C = struct();
C.time   = time;
C.pelvis = zeros(nT, 1);
C.calcnR = zeros(nT, 1);
C.calcnL = zeros(nT, 1);
C.modelFile = char(modelFile);

for k = 1:nT
    % 设置所有模型坐标；qNames 里没有的坐标设为 0。
    for i = 1:nCoord
        if coordCol(i) > 0
            coordHandles{i}.setValue(state, q(k, coordCol(i)), false);
            coordHandles{i}.setSpeedValue(state, qdot(k, coordCol(i)));
        else
            coordHandles{i}.setValue(state, 0, false);
            coordHandles{i}.setSpeedValue(state, 0);
        end
    end
    model.realizePosition(state);

    % 与 cal_gait_phase.m 完全相同的调用；get(0) 是 Vec3 的 X 分量。
    C.pelvis(k) = bodyPelvis.findStationLocationInGround( ...
        state, pelvisCOM_local).get(0);
    C.calcnR(k) = bodyCalcnR.findStationLocationInGround( ...
        state, calcnRCOM_local).get(0);
    C.calcnL(k) = bodyCalcnL.findStationLocationInGround( ...
        state, calcnLCOM_local).get(0);
end
end
