function modelInfo = cal_grf(simCache, modelInfo, state, frameIndex)
% Name: cal_grf
% Description: read the foot-ground contact forces of the model and reduce
%   them to ONE resultant ground reaction force (GRF) plus ONE centre of
%   pressure (CoP) per foot, both expressed in the ground frame.
%
%   In human0918.osim each foot's contact is split into three independent
%   HuntCrossleyForce components (heel / toe / tip, see foot_contact_defs).
%   OpenSim's getRecordValues() returns the force the contact sphere exerts
%   ON the ground, so the reaction force on the body is the negated vector
%   and the resultant is a plain sum:
%       F_grf = - sum_k F_ground,k
%
%   The CoP is the normal-force-weighted application point of the resultant:
%       x_cop = sum_k (Fy_k * x_k) / sum_k Fy_k ,   y_cop = 0
%   where x_k is the horizontal ground position of contact point k and the
%   CoP is taken to lie ON the ground plane (y = 0), which is where the
%   contact points actually are.  The CoP is undefined (NaN) while the foot
%   is unloaded, i.e. during swing.

contactDef = simCache.contactDef;
nSides  = numel(contactDef.sides);
nPts    = contactDef.nPoints;
minLoad = 1e-3;   % [N] below this the foot counts as unloaded -> CoP = NaN

fx   = zeros(1, nSides);
fy   = zeros(1, nSides);
copx = nan(1, nSides);

for i = 1 : nSides
    pts = simCache.footContact{i};
    sumFy   = 0;   % sum of the normal components (CoP weighting)
    sumMx   = 0;   % sum of Fy_k * x_k
    for j = 1 : nPts
        f   = pts(j).force.getRecordValues(state);
        fxk = -f.get(0);
        fyk = -f.get(1);

        fx(i) = fx(i) + fxk;
        fy(i) = fy(i) + fyk;

        if fyk > 0   % only the (compressive) normal part defines the CoP
            pg     = pts(j).frame.findStationLocationInGround(state, pts(j).loc);
            sumMx  = sumMx + fyk * pg.get(0);
            sumFy  = sumFy + fyk;
        end
    end

    if sumFy > minLoad
        copx(i) = sumMx / sumFy;
    end
end

% side 1 = right, side 2 = left (see foot_contact_defs)
modelInfo.dy.grf.fxr(frameIndex)   = fx(1);
modelInfo.dy.grf.fyr(frameIndex)   = fy(1);
modelInfo.dy.grf.copxr(frameIndex) = copx(1);
modelInfo.dy.grf.fxl(frameIndex)   = fx(2);
modelInfo.dy.grf.fyl(frameIndex)   = fy(2);
modelInfo.dy.grf.copxl(frameIndex) = copx(2);

end

