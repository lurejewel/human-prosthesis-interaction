% -------------------------------------------------------------------------
% Name: ModelInfo.m
% Description: definition of modelInfo Class.
% Member variables:
% 
% Member functions:
% - ModelInfo(): construction function.
% - 
% -------------------------------------------------------------------------

classdef ModelInfo

    properties
        % muscle states
        muscleExcitations
        muscleActivations
        muscleForces
        muscleFiberForcesATN % ATN: along tendon, normalized
        muscleFiberLengthN
        muscleReflexDelay
        muscleReflex
        % gait phases
        phaseR
        phaseL
        % others
        grf
        state
        time
            end

    methods
        function obj = ModelInfo(model, simConfig, muscleReflexParaArray)

            % number of muscles
            muscleNum = model.getMuscles.getSize;

            % number of points
            npts = ceil((simConfig.endTime - simConfig.startTime) / simConfig.stepTime) + 1;

            % muslce states: excitation, activation, muscle fiber forces along tendon(or not?)
            delay = round(0.01/simConfig.stepTime); % 10 ms delay for muscle reflex mechanism to take effect
            obj.muscleReflexDelay = delay;
            obj.muscleExcitations = nan(muscleNum, npts+delay);
            obj.muscleActivations = nan(muscleNum, npts);
            obj.muscleForces = nan(muscleNum, npts);
            obj.muscleFiberForcesATN = nan(muscleNum, npts);
            obj.muscleFiberLengthN = nan(muscleNum, npts);

            % muscle reflex parameters
            obj.muscleReflex.tib.KL                     = muscleReflexParaArray(1);
            obj.muscleReflex.tib.L0                     = muscleReflexParaArray(2);
            obj.muscleReflex.tib_sol.KF                 = muscleReflexParaArray(3);
            obj.muscleReflex.sol.KF                     = muscleReflexParaArray(4);
            obj.muscleReflex.gas.KF                     = muscleReflexParaArray(5);
            obj.muscleReflex.ili_pelvis_tilt.KP         = muscleReflexParaArray(6);
            obj.muscleReflex.ili_pelvis_tilt.KV         = muscleReflexParaArray(7);
            obj.muscleReflex.ili_pelvis_tilt.C0         = muscleReflexParaArray(8);
            obj.muscleReflex.ili.C0                     = muscleReflexParaArray(9);
            obj.muscleReflex.ili.KL                     = muscleReflexParaArray(10);
            obj.muscleReflex.ili.L0                     = muscleReflexParaArray(11);
            obj.muscleReflex.ili_pelvis_tilt.P02        = muscleReflexParaArray(12);
            obj.muscleReflex.ili_pelvis_tilt.KP2        = muscleReflexParaArray(13);
            obj.muscleReflex.ili_pelvis_tilt.KV2        = muscleReflexParaArray(14);
            obj.muscleReflex.ili_ham.KL                 = muscleReflexParaArray(15);
            obj.muscleReflex.ili_ham.L0                 = muscleReflexParaArray(16);
            obj.muscleReflex.ham_pelvis_tilt.KP         = muscleReflexParaArray(17);
            obj.muscleReflex.ham_pelvis_tilt.KV         = muscleReflexParaArray(18);
            obj.muscleReflex.ham_pelvis_tilt.C0         = muscleReflexParaArray(19);
            % obj.muscleReflex.ham.KF                     = muscleReflexParaArray(20);
            obj.muscleReflex.ham_glu.KF                 = muscleReflexParaArray(20);
            obj.muscleReflex.glu_pelvis_tilt.KP         = muscleReflexParaArray(21);
            obj.muscleReflex.glu_pelvis_tilt.KV         = muscleReflexParaArray(22);
            obj.muscleReflex.glu_pelvis_tilt.C0         = muscleReflexParaArray(23);
            obj.muscleReflex.glu.KF                     = muscleReflexParaArray(24);
            % obj.muscleReflex.glu.C0                     = muscleReflexParaArray(25);
            obj.muscleReflex.vas.KF1                    = muscleReflexParaArray(25); % 26 -> 25
            obj.muscleReflex.vas.KF2                    = muscleReflexParaArray(26);
            obj.muscleReflex.vas.C0                     = muscleReflexParaArray(27);
            obj.muscleReflex.vas_knee.pos_max           = muscleReflexParaArray(28);

            % phases
            obj.phaseR = nan(1, npts);
            obj.phaseL = nan(1, npts);

            % ground reaction forces
            obj.grf.normalR = nan(1, npts);
            obj.grf.normalL = nan(1, npts);
            obj.grf.frictionR = nan(1, npts);
            obj.grf.frictionL = nan(1, npts);

            % time series
            obj.time = simConfig.startTime : simConfig.stepTime : simConfig.endTime;

        end
    end

end
