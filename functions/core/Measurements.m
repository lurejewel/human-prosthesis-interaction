% -------------------------------------------------------------------------
% Name: Measurements.m
% Description: definition of of measurements Class.
% Member variables:
%
% Member functions:
% - Measurements(): construction function.
%
% - gait_completeness_measure(): fitness evaluation at stage 1. The longer
% distance the model achieves, better fitness it gets.
%
% - gait_velocity_measure(weight): term that evaluates whether the model
% walks/runs near the desired velocity.
%
% - knee_limit_measure(weight): term that evaluates the coordinate limit forces
% of the knee joint, which occurs at hyperextension / hyperflexion.
%
% - ankle_limit_measure(weight): term that evaluates the coordinate limit forces
% of the ankle joint, which occurs at hyperextension / hyperflexion.
%
% - limit_force_calculation(ang, vel, angUpTh, angLowTh): calculate
% coordinate limit force given the joint angle and joint velocity in a
% piece of time.
% -------------------------------------------------------------------------

classdef Measurements < handle

    properties
        modelInfo   % handle reference to ModelInfo (no deep-copy)
        mass
        gravity

        simConfig
        lastTime
        npts
        heelstrikeEventR
        heelstrikeEventL
        muscleMass
        muscleLambda
        distance

    end

    methods
        function obj = Measurements(modelInfo)
            obj.modelInfo = modelInfo;  % store handle, no copy
            obj.mass    = modelInfo.st.model.totalMass;
            obj.gravity = modelInfo.st.model.gravity;
            obj.simConfig = modelInfo.st.simInfo;
            obj.lastTime  = modelInfo.dy.lastTime;
            obj.npts = find(all(~isnan(modelInfo.dy.labelHistory)), 1, 'last');

            % obj.muscleMass   = modelInfo.st.muscle.mass;
            % obj.muscleLambda = modelInfo.st.muscle.lambda;

            % pre-allocate and collect heelstrike events
            obj.heelstrikeEventR = zeros(1, obj.npts);
            obj.heelstrikeEventL = zeros(1, obj.npts);
            nHS_R = 0; nHS_L = 0;
            for i = 2 : obj.npts
                if modelInfo.dy.phase.r(i-1) == 4 && modelInfo.dy.phase.r(i) == 0
                    nHS_R = nHS_R + 1;
                    obj.heelstrikeEventR(nHS_R) = i;
                end
                if modelInfo.dy.phase.l(i-1) == 4 && modelInfo.dy.phase.l(i) == 0
                    nHS_L = nHS_L + 1;
                    obj.heelstrikeEventL(nHS_L) = i;
                end
            end
            obj.heelstrikeEventR = obj.heelstrikeEventR(1:nHS_R);
            obj.heelstrikeEventL = obj.heelstrikeEventL(1:nHS_L);

            obj.distance = modelInfo.dy.labelHistory( ...
                modelInfo.st.model.map('pelvis_tx/value'), obj.npts);

        end

        function fit = gait_completeness_measure(obj, weight)
            fit = weight * (obj.simConfig.endTime - obj.lastTime)^2;
        end

        function fit = gait_velocity_measure(obj, weight)
            th = 0.05;
            nR = length(obj.heelstrikeEventR);
            nL = length(obj.heelstrikeEventL);
            velBuffer = nan(nR + nL - 2, 1);
            if isempty(velBuffer)
                fit = weight * 1;
                return;
            end

            sh  = obj.modelInfo.dy.labelHistory;
            sm  = obj.modelInfo.st.model.map;
            idx = sm('pelvis_tx/value');
            dt  = obj.simConfig.stepTime;

            for i = 1 : nR - 1
                velBuffer(i) = (sh(idx, obj.heelstrikeEventR(i+1)) - sh(idx, obj.heelstrikeEventR(i))) ...
                             / ((obj.heelstrikeEventR(i+1) - obj.heelstrikeEventR(i)) * dt);
            end
            offset = max(0, nR - 1);
            for i = 1 : nL - 1
                velBuffer(i + offset) = (sh(idx, obj.heelstrikeEventL(i+1)) - sh(idx, obj.heelstrikeEventL(i))) ...
                                      / ((obj.heelstrikeEventL(i+1) - obj.heelstrikeEventL(i)) * dt);
            end
            velErrAvgN = mean(abs(velBuffer - obj.simConfig.speed)) / obj.simConfig.speed;
            fit = weight * max(0, velErrAvgN - th);
        end

        function fit = knee_limit_measure(obj, weight)
            th = 5; angUpTh = -5; angLowTh = -120;

            sh  = obj.modelInfo.dy.labelHistory;
            sm  = obj.modelInfo.st.model.map;
            n   = obj.npts;

            kneeAngLeft  = rad2deg(sh(sm('knee_extension_l/value'), 1:n));
            kneeAngRight = rad2deg(sh(sm('knee_extension_r/value'), 1:n));
            kneeVelLeft  = rad2deg(sh(sm('knee_extension_l/speed'), 1:n));
            kneeVelRight = rad2deg(sh(sm('knee_extension_r/speed'), 1:n));

            fAvg = mean(abs(Measurements.limit_force_calculation(kneeAngLeft, kneeVelLeft, angUpTh, angLowTh)));
            fitL = weight * max(0, fAvg - th);
            fAvg = mean(abs(Measurements.limit_force_calculation(kneeAngRight, kneeVelRight, angUpTh, angLowTh)));
            fitR = weight * max(0, fAvg - th);
            fit = fitL + fitR;
        end

        function fit = ankle_limit_measure(obj, weight)
            th = 5; angUpTh = 20; angLowTh = -45;

            sh  = obj.modelInfo.dy.labelHistory;
            sm  = obj.modelInfo.st.model.map;
            n   = obj.npts;

            ankleAngLeft  = rad2deg(sh(sm('ankle_dorsiflexion_l/value'), 1:n));
            ankleAngRight = rad2deg(sh(sm('ankle_dorsiflexion_r/value'), 1:n));
            ankleVelLeft  = rad2deg(sh(sm('ankle_dorsiflexion_l/speed'), 1:n));
            ankleVelRight = rad2deg(sh(sm('ankle_dorsiflexion_r/speed'), 1:n));

            fAvg = mean(abs(Measurements.limit_force_calculation(ankleAngLeft, ankleVelLeft, angUpTh, angLowTh)));
            fitL = weight * max(0, fAvg - th);
            fAvg = mean(abs(Measurements.limit_force_calculation(ankleAngRight, ankleVelRight, angUpTh, angLowTh)));
            fitR = weight * max(0, fAvg - th);
            fit = fitL + fitR;
        end

        function fit = grf_limit_measure(obj, weight)
            th = 1.4;
            n  = obj.npts;
            m  = obj.mass;
            g  = obj.gravity;
            grfRN = obj.modelInfo.dy.grf.fyr(1:n) / m / g;
            grfLN = obj.modelInfo.dy.grf.fyl(1:n) / m / g;
            fit = weight * (mean((grfLN - th) .* (grfLN > th)) ...
                          + mean((grfRN - th) .* (grfRN > th)));
        end

        function fit = effort_measure(obj, weight)

            sh      = obj.modelInfo.dy.labelHistory;
            sm      = obj.modelInfo.st.model.map;
            pelvisX = sh(sm('pelvis_tx/value'), obj.npts);
            m       = obj.mass;

            fit = weight * obj.modelInfo.dy.totalMetabolic / m / pelvisX;

            % nMuscles = numel(obj.muscleMass);
            % Effort   = zeros(1, obj.npts);
            % mi       = obj.modelInfo.dy.muscle;
            % massVec  = obj.muscleMass;
            % lambdaVec = obj.muscleLambda;
            % 
            % for i = 2 : obj.npts
            %     dE = zeros(1, nMuscles);
            %     for j = 1 : nMuscles
            %         m  = massVec(j);
            %         la = lambdaVec(j);
            %         u  = mi.exc(j, i);
            %         a  = mi.act(j, i);
            %         lCEN = mi.lCEN(j, i);
            %         vCE  = max(-mi.vCE(j, i), 0);
            %         fMTU = mi.fMTU(j, i);
            %         fCE  = mi.fCE(j, i);
            % 
            %         dA = m * Measurements.fAFcn(la, u);
            %         dM = m * Measurements.gFcn(lCEN) * Measurements.fMFcn(la, a);
            %         dS = 0.25 * fMTU * vCE;
            %         dW = fCE * vCE;
            %         dE(j) = dA + dM + dS + dW;
            %     end
            %     Effort(i) = sum(dE);
            % end
            % 
            % fit = weight * (mean(Effort) / obj.mass + 1.51);
        end
    end
    methods (Static)
        function force = limit_force_calculation(ang, vel, angUpTh, angLowTh)
            kUp = 10;
            kLow = 10;
            damp = 0.2;

            kUp_ = kUp * ( atan(10*(ang-angUpTh))/pi + 0.5 );
            kLow_ = kLow * ( -atan(10*(ang-angLowTh))/pi + 0.5 );
            fUpLim = kUp_ .* ( angUpTh - ang );
            fLowLim = kLow_ .* ( angLowTh - ang );
            fDamp = -damp * (kUp_/kUp + kLow_/kLow) .* vel;
            force = fUpLim + fLowLim + fDamp;

        end

        function out = fAFcn(lambda, u)
            out = 40 * lambda * sin(pi*u/2) + 133 * (1-lambda) * (1-cos(pi*u/2));

        end

        function out = gFcn(lCE)
            if lCE > 0 && lCE <= 0.5
                out = 0.5;
            elseif lCE > 0.5 && lCE <= 1.0
                out = lCE;
            elseif lCE > 1.0 && lCE <= 1.5
                out = -2*lCE + 3;
            elseif lCE > 1.5
                out = 0;
            else
                error('muscle fiber length cannot be non-positive.');

            end

        end

        function out = fMFcn(lambda, a)
            out = 74 * lambda * sin(pi*a/2) + 111 * (1-lambda) * (1-cos(pi*a/2));

        end

    end
end