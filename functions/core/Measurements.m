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

classdef Measurements

    properties
        mass
        gravity

        simConfig
        lastTime % duration of the simulation
        npts
        state
        heelstrikeEventR
        heelstrikeEventL
        grf = struct('normalR', [], 'normalL', [])
        stateMap
        muscleinfo
        distance

    end

    methods
        function obj = Measurements(modelInfo)
            obj.mass = modelInfo.st.model.totalMass;
            obj.gravity = 9.80665;
            obj.simConfig = modelInfo.st.simInfo;
            obj.lastTime = modelInfo.dy.lastTime;
            obj.npts = find(all(~isnan(modelInfo.dy.stateHistory)),1,'last');

            obj.stateMap = modelInfo.st.model.map;
            obj.state = modelInfo.dy.stateHistory(:,1:obj.npts);
            
            obj.grf.normalR = modelInfo.dy.grf.fyr(1:obj.npts);
            obj.grf.normalL = modelInfo.dy.grf.fyl(1:obj.npts);

            obj.muscleinfo.mass = modelInfo.st.muscle.mass;
            obj.muscleinfo.lambda = modelInfo.st.muscle.lambda;
            obj.muscleinfo.u = modelInfo.dy.muscle.exc(:,1:obj.npts);
            obj.muscleinfo.a = modelInfo.dy.muscle.act(:,1:obj.npts);
            obj.muscleinfo.F_MTU = modelInfo.dy.muscle.fMTU(:,1:obj.npts);
            obj.muscleinfo.F_CE = modelInfo.dy.muscle.fCE(:,1:obj.npts);
            obj.muscleinfo.l_CEN = modelInfo.dy.muscle.lCEN(:,1:obj.npts);
            obj.muscleinfo.v_CE = modelInfo.dy.muscle.vCE(:,1:obj.npts);

            obj.distance = obj.state(obj.stateMap('pelvis_tx/value'),end);

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

        end

        function fit = gait_completeness_measure(obj, weight)
            % fit = obj.simConfig.endTime*obj.simConfig.speed - obj.distance; % desired distance - actual distance (in meters)
            fit = weight * (obj.simConfig.endTime - obj.lastTime)^2; %

        end

        function fit = gait_velocity_measure(obj, weight)
            th = 0.05; % threshold of velocity range; if the error between the desired velocity and the actual walking velocity is less than the threshold, then no penalty is imposed
            velBuffer = nan(length(obj.heelstrikeEventR)+length(obj.heelstrikeEventL)-2, 1);
            if isempty(velBuffer)
                fit = weight * 1;
            else
                for i = 1 : length(obj.heelstrikeEventR)-1
                    velBuffer(i) = ( obj.state(obj.stateMap('pelvis_tx/value'),obj.heelstrikeEventR(i+1)) - obj.state(obj.stateMap('pelvis_tx/value'),obj.heelstrikeEventR(i)) ) / ( (obj.heelstrikeEventR(i+1)-obj.heelstrikeEventR(i)) * obj.simConfig.stepTime );
                end
                if isempty(i)
                    i = 0;
                end
                i_ = i;
                for i = 1 : length(obj.heelstrikeEventL)-1
                    velBuffer(i+i_) = ( obj.state(obj.stateMap('pelvis_tx/value'),obj.heelstrikeEventL(i+1)) - obj.state(obj.stateMap('pelvis_tx/value'),obj.heelstrikeEventL(i)) ) / ( (obj.heelstrikeEventL(i+1)-obj.heelstrikeEventL(i)) * obj.simConfig.stepTime );
                end
                velErrAvgN = mean(abs(velBuffer-obj.simConfig.speed)) / obj.simConfig.speed;
                fit = weight * velErrAvgN * (velErrAvgN > th);
            end
            % velErrSum = 0; % initialize sum of error of walking velocity
            % for i = 1 : obj.npts
            %     velErr = obj.state(obj.stateMap('pelvis_tx/speed'),i) - obj.simConfig.speed;
            %     velErrSum = velErrSum + velErr * (abs(velErr)>th);
            % end
            % velErrAvg = velErrSum / obj.npts; % average error of walking velocity during the simulation
            %
            % fit = weight * (velErrAvg + obj.npts<obj.lastTime/obj.simConfig.stepTime-1);

        end

        function fit = knee_limit_measure(obj, weight)
            th = 5; % in Newton; if the average knee limit force is less than th, then no penalty is imposed
            angUpTh = -3; % in deg; upper bound of the knee angle
            angLowTh = -120; % in deg; lower bound of the knee angle
            % kUp = 2;
            % kLow = 2;
            % damp = 0.2;
            kneeAngLeft = rad2deg(obj.state(obj.stateMap('knee_flexion_l/value'),:));
            kneeAngRight = rad2deg(obj.state(obj.stateMap('knee_flexion_r/value'),:));
            kneeVelLeft = rad2deg(obj.state(obj.stateMap('knee_flexion_l/speed'),:));
            kneeVelRight = rad2deg(obj.state(obj.stateMap('knee_flexion_r/speed'),:));

            % kUp_ = kUp * ( atan(10*(kneeAngLeft-angUpTh))/pi + 0.5 );
            % kLow_ = kLow * ( atan(10*(kneeAngLeft-angLowTh))/pi + 0.5 );
            % fUpLim = kUp_ .* ( angUpTh - kneeAngLeft );
            % fLowLim = kLow_ .* ( angLowTh - kneeAngLeft );
            % fDamp = -damp * (kUp_/kUp + kLow_/kLow_) .* kneeVelLeft;
            % f = fUpLim + fLowLim + fDamp;
            %
            % fAvg = sum(abs(f))/length(kneeAngLeft); % same as mean(abs(f))?
            fAvg = mean(abs(Measurements.limit_force_calculation(kneeAngLeft, kneeVelLeft, angUpTh, angLowTh)));
            fitL = weight * fAvg * (fAvg>th);

            % kUp_ = kUp * ( atan(10*(kneeAngRight-angUpTh))/pi + 0.5 );
            % kLow_ = kLow * ( atan(10*(kneeAngRight-angLowTh))/pi + 0.5 );
            % fUpLim = kUp_ .* ( angUpTh - kneeAngRight );
            % fLowLim = kLow_ .* ( angLowTh - kneeAngRight );
            % fDamp = -damp * (kUp_/kUp + kLow_/kLow_) .* kneeVelRight;
            % f = fUpLim + fLowLim + fDamp;
            %
            % fAvg = sum(abs(f))/length(kneeAngRight); % same as mean(abs(f))?
            fAvg = mean(abs(Measurements.limit_force_calculation(kneeAngRight, kneeVelRight, angUpTh, angLowTh)));
            fitR = weight * fAvg * (fAvg>th);

            fit = fitL + fitR;

        end

        function fit = ankle_limit_measure(obj, weight)
            th = 5;
            angUpTh = 20;
            angLowTh = -45;

            ankleAngLeft = rad2deg(obj.state(obj.stateMap('ankle_dorsiflexion_l/value'),:));
            ankleAngRight = rad2deg(obj.state(obj.stateMap('ankle_dorsiflexion_r/value'),:));
            ankleVelLeft = rad2deg(obj.state(obj.stateMap('ankle_dorsiflexion_l/speed'),:));
            ankleVelRight = rad2deg(obj.state(obj.stateMap('ankle_dorsiflexion_r/speed'),:));

            fAvg = mean(abs(Measurements.limit_force_calculation(ankleAngLeft, ankleVelLeft, angUpTh, angLowTh)));
            fitL = weight * fAvg * (fAvg>th);

            fAvg = mean(abs(Measurements.limit_force_calculation(ankleAngRight, ankleVelRight, angUpTh, angLowTh)));
            fitR = weight * fAvg * (fAvg>th);

            fit = fitL + fitR;

        end

        function fit = grf_limit_measure(obj, weight)
            th = 1.4; % penalty imposed when vertical force exceeds thx body weights
            grfRN = obj.grf.normalR / obj.mass / obj.gravity;
            grfLN = obj.grf.normalL / obj.mass / obj.gravity;

            fitL = weight * mean((grfLN-th).*(grfLN>th));
            fitR = weight * mean((grfRN-th).*(grfRN>th));

            fit = fitL + fitR;

        end

        function fit = effort_measure(obj, weight)

            Effort = zeros(1, obj.npts);
            nMuscles = numel(obj.muscleinfo.mass);
            massVec  = obj.muscleinfo.mass;
            lambdaVec = obj.muscleinfo.lambda;

            for i = 2 : obj.npts
                dE = zeros(1, nMuscles);
                for j = 1 : nMuscles
                    m = massVec(j);
                    lambda = lambdaVec(j);
                    u = obj.muscleinfo.u(j,i);
                    a = obj.muscleinfo.a(j,i);
                    lCEN = obj.muscleinfo.l_CEN(j,i);
                    vCE = max(-obj.muscleinfo.v_CE(j,i), 0);
                    fMTU = obj.muscleinfo.F_MTU(j,i);
                    fCE = obj.muscleinfo.F_CE(j,i);

                    % dE = dA + dM + dS + dW
                    dA = m * Measurements.fAFcn(lambda, u); % muscle activation heat rate
                    dM = m * Measurements.gFcn(lCEN) * Measurements.fMFcn(lambda, a); % muscle maintenance heat rate
                    dS = 0.25 * fMTU *vCE; % muscle shortening heat rate
                    dW = fCE * vCE; % positive mechanical work rate
                    dE(j) = dA + dM + dS + dW; % muscular metabolic expenditure of the j-th muscle at i-th timestep

                end
                Effort(i) = sum(dE); % overall muscular metabolic expenditure at i-th timestep

            end

            fit = weight * (mean(Effort)/obj.mass + 1.51); % 1.51: basal metabolic energy rate

        end
    end
    methods (Static)
        function force = limit_force_calculation(ang, vel, angUpTh, angLowTh)
            kUp = 2;
            kLow = 2;
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