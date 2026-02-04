function modelInfo = cal_gait_phase(model, modelInfo, state, frameIndex)
%% materials for gait phase detection
% state, grf
% state = modelInfo.state;
GRF = [modelInfo.dy.grf.fyr(frameIndex), modelInfo.dy.grf.fyl(frameIndex)];

% init gait phases
phaseR = -1;
phaseL = -1;

% body weight, stance threshold
BW = -model.getTotalMass(state) * model.getGravity.get(1); % body weight (N)
stanceTh = 0.23137978; % for heel strike detection

% Normal GRF
grfR = GRF(1)/BW; % normal component of GRF, right-foot side
grfL = GRF(2)/BW; % left-foot side

% X-axis of pelvis COM
pelvis = model.getBodySet().get('pelvis');
pelvisCOMLocal = pelvis.getMassCenter();
pelvisCOMGlobal = pelvis.findStationLocationInGround(state, pelvisCOMLocal);
pelvisCOM = pelvisCOMGlobal.get(0);

% X-axis of calcn COM
calcnR = model.getBodySet().get('calcn_r');
calcnRCOMLocal = calcnR.getMassCenter();
calcnRCOMGlobal = calcnR.findStationLocationInGround(state, calcnRCOMLocal);
calcnRCOM = calcnRCOMGlobal.get(0);
calcnL = model.getBodySet().get('calcn_l');
calcnLCOMLocal = calcnL.getMassCenter();
calcnLCOMGlobal = calcnL.findStationLocationInGround(state, calcnLCOMLocal);
calcnLCOM = calcnLCOMGlobal.get(0);

%% gait detection

if frameIndex == 1
% first frame:
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

    % after first frame, using iterative method:
    % [Early Stance->Late Stance]: calcnCOM < pelvisCOM
    % [Late Stance->Liftoff]:      opposite GRF >= stanceTh || calcnCOM+1 < pelvisCOM
    % [Liftoff->Swing]:            GRF < stanceTh
    % [Swing->Landing]:            calcnCOM >= pelvisCOM
    % [Landing->Early Stance]:     GRF >= stanceTh
elseif frameIndex > 1

    phaseR = modelInfo.dy.phase.r(frameIndex-1);
    phaseL = modelInfo.dy.phase.l(frameIndex-1);

    % right
    switch phaseR
        case 0 % Early Stance -> Late Stance
            if calcnRCOM < pelvisCOM
                phaseR = 1;
            end
        case 1 % Late Stance -> Liftoff
            if grfL >= stanceTh || calcnRCOM+1 < pelvisCOM
                phaseR = 2;
            end
        case 2 % Liftoff -> Swing
            if grfR < stanceTh
                phaseR = 3;
            end
        case 3 % Swing -> Landing
            if calcnRCOM > pelvisCOM
                phaseR = 4;
            end
        case 4 % Landing -> Early Stance
            if grfR >= stanceTh
                phaseR = 0;
            end
        otherwise
            error('invalid gait phase.')
    end
    
    % left
    switch phaseL
        case 0 % Early Stance -> Late Stance
            if calcnLCOM < pelvisCOM
                phaseL = 1;
            end
        case 1 % Late Stance -> Liftoff
            if grfR >= stanceTh || calcnLCOM+1 < pelvisCOM
                phaseL = 2;
            end
        case 2 % Liftoff -> Swing
            if grfL < stanceTh
                phaseL = 3;
            end
        case 3 % Swing -> Landing
            if calcnLCOM > pelvisCOM
                phaseL = 4;
            end
        case 4 % Landing -> Early Stance
            if grfL >= stanceTh
                phaseL = 0;
            end
        otherwise
            error('invalid gait phase.')
    end

else
    error('invalid frame index.')
end

% check
if phaseR == -1 || phaseL == -1
    error('invalid gait phase detected.')
end

modelInfo.dy.phase.r(frameIndex) = phaseR;
modelInfo.dy.phase.l(frameIndex) = phaseL;

end