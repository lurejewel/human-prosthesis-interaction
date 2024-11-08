function [phaseR, phaseL] = cal_gait_phase(model, state, GRF)
%% materials for gait phase detection
% init gait phases
phaseR = -1;
phaseL = -1;

% body weight, stance threshold
BW = -model.getTotalMass(state) * model.getGravity.get(1); % body weight (N)
stanceTh = 0.23137978; % for heel strike detection

% Normal GRF (暂时先通过地反力文件读取，后期改为接触力反算的形式)
grfR = GRF(1)/BW; % normal component of GRF, right-foot side
grfL = GRF(2)/BW; % left-foot side

% X-axis of pelvis COM
pelvis = model.getBodySet().get('pelvis');
pelvisCOMLocal = pelvis.getMassCenter();
pelvisCOMGlobal = pelvis.findStationLocationInGround(state, pelvisCOMLocal);
pelvisCOM = pelvisCOMGlobal.get(0);

% X axis of calcn COM
calcnR = model.getBodySet().get('calcn_r');
calcnRCOMLocal = calcnR.getMassCenter();
calcnRCOMGlobal = calcnR.findStationLocationInGround(state, calcnRCOMLocal);
calcnRCOM = calcnRCOMGlobal.get(0);
calcnL = model.getBodySet().get('calcn_l');
calcnLCOMLocal = calcnL.getMassCenter();
calcnLCOMGlobal = calcnL.findStationLocationInGround(state, calcnLCOMLocal);
calcnLCOM = calcnLCOMGlobal.get(0);

%% gait detection
% [Early Stance]: GRF >= stanceTh && calcnCOM >= pelvisCOM
% [Late Stance]:  GRF >= stanceTh && calcnCOM < pelvisCOM && opposite GRF <
% stanceTh
% [Liftoff]:      GRF >= stanceTh && calcnCOM < pelvisCOM && (opposite GRF 
% > stanceTh || calcnCOM+1 < pelvisCOM)
% [Swing]:        GRF < stanceTh && calcnCOM < pelvisCOM
% [Landing]:      GRF < stanceTh && calcnCOM >= pelvisCOM

% right
if grfR >= stanceTh
    if calcnRCOM >= pelvisCOM
        phaseR = 0; % Early Stance
    elseif grfL > stanceTh || calcnRCOM+1 < pelvisCOM
        phaseR = 2; % Liftoff
    elseif grfL < stanceTh
        phaseR = 1; % Late Stance
    end
elseif calcnRCOM < pelvisCOM
    phaseR = 3; % Swing
elseif calcnRCOM >= pelvisCOM
    phaseR = 4; % Landing
end

% left
if grfL >= stanceTh
    if calcnLCOM >= pelvisCOM
        phaseL = 0; % Early Stance
    elseif grfR > stanceTh || calcnLCOM+1 < pelvisCOM
        phaseL = 2; % LiftOff
    elseif grfR < stanceTh
        phaseL = 1; % Late Stance
    end
elseif calcnLCOM < pelvisCOM
    phaseL = 3; % Swing
elseif calcnLCOM >= pelvisCOM
    phaseL = 4; % Landing
end

% check
if phaseR == -1 || phaseL == -1
    error('invalid gait phase detected.')
end

end