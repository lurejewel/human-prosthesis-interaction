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
        % model information
        mass
        g
        % muscle-related information
        muscleNames
        muscleExcitations
        muscleActivations
        muscleForces
        muscleFiberForcesATN % ATN: along tendon, normalized
        muscleFiberLengthN
        muscleReflexDelay
        muscleReflex
        OptimalFiberLengths
        OptimalFiberForces
        % gait phases
        phaseR
        phaseL
        % others
        grf
        % state
        stateHistory
        time
    end

    methods
        function obj = ModelInfo(model, simConfig, muscleReflexParaArray)

            % number of muscles
            muscleNum = model.getMuscles.getSize;

            % number of points
            npts = ceil((simConfig.endTime - simConfig.startTime) / simConfig.stepTime)+1;

            % muslce states: excitation, activation, muscle fiber forces along tendon(or not?)
            delay = round(0.01/simConfig.stepTime); % 10 ms delay for muscle reflex mechanism to take effect
            obj.muscleNames = cell(muscleNum, 1);
            obj.muscleReflexDelay = delay;
            obj.muscleExcitations = nan(muscleNum, npts+delay);
            obj.muscleActivations = nan(muscleNum, npts);
            obj.muscleForces = nan(muscleNum, npts);
            obj.muscleFiberForcesATN = nan(muscleNum, npts);
            obj.muscleFiberLengthN = nan(muscleNum, npts);

            % muscle names
            for muscleIndex = 0 : muscleNum-1
                obj.muscleNames{muscleIndex+1} = char(model.getMuscles.get(muscleIndex).getName);
            end

            % muscle reflex parameters
            obj.muscleReflex = obj.read_muscleReflex_array(muscleReflexParaArray);

            % optimal fiber lengths后续可以改成数组形式
            obj.OptimalFiberLengths.hamstrings = model.getMuscles.get('hamstrings_r').getOptimalFiberLength;
            obj.OptimalFiberLengths.glut_max = model.getMuscles.get('glut_max_r').getOptimalFiberLength;
            obj.OptimalFiberLengths.ilipsoas = model.getMuscles.get('ilipsoas_r').getOptimalFiberLength;
            obj.OptimalFiberLengths.vasti = model.getMuscles.get('vasti_r').getOptimalFiberLength;
            obj.OptimalFiberLengths.gastrocnemius = model.getMuscles.get('gastroc_r').getOptimalFiberLength;
            obj.OptimalFiberLengths.soleus = model.getMuscles.get('soleus_r').getOptimalFiberLength;
            obj.OptimalFiberLengths.tibia = model.getMuscles.get('tibia_r').getOptimalFiberLength;

            % optimal fiber forces
            obj.OptimalFiberForces.hamstrings = model.getMuscles.get('hamstrings_r').getMaxIsometricForce;
            obj.OptimalFiberForces.glut_max = model.getMuscles.get('glut_max_r').getMaxIsometricForce;
            obj.OptimalFiberForces.ilipsoas = model.getMuscles.get('ilipsoas_r').getMaxIsometricForce;
            obj.OptimalFiberForces.vasti = model.getMuscles.get('vasti_r').getMaxIsometricForce;
            obj.OptimalFiberForces.gastrocnemius = model.getMuscles.get('gastroc_r').getMaxIsometricForce;
            obj.OptimalFiberForces.soleus = model.getMuscles.get('soleus_r').getMaxIsometricForce;
            obj.OptimalFiberForces.tibia = model.getMuscles.get('tibia_r').getMaxIsometricForce;

            % phases
            obj.phaseR = nan(1, npts);
            obj.phaseL = nan(1, npts);

            % ground reaction forces
            obj.grf.normalR = nan(1, npts);
            obj.grf.normalL = nan(1, npts);
            obj.grf.frictionR = nan(1, npts);
            obj.grf.frictionL = nan(1, npts);

            % state history
            state = model.initSystem(); % this is a must before calling model.getNumStateVariables()
            obj.stateHistory = nan(model.getNumStateVariables, npts);

            % time series
            obj.time = simConfig.startTime : simConfig.stepTime : simConfig.endTime;

            % basic model information
            obj.mass = model.getTotalMass(state);
            obj.g = -model.getGravity.get(1);
            
        end

        function muscleReflex = read_muscleReflex_array(~, muscleReflexParaArray)
            % Name: read_muscleReflex_array
            % Description: convert the muscle-reflex parameters from array 
            % form to struct form.
            muscleReflex.tib.KL                     = muscleReflexParaArray(1);
            muscleReflex.tib.L0                     = muscleReflexParaArray(2);
            muscleReflex.tib_sol.KF                 = muscleReflexParaArray(3);
            muscleReflex.sol.KF                     = muscleReflexParaArray(4);
            muscleReflex.gas.KF                     = muscleReflexParaArray(5);
            muscleReflex.ili_pelvis_tilt.KP         = muscleReflexParaArray(6);
            muscleReflex.ili_pelvis_tilt.KV         = muscleReflexParaArray(7);
            muscleReflex.ili_pelvis_tilt.C0         = muscleReflexParaArray(8);
            muscleReflex.ili.C0                     = muscleReflexParaArray(9);
            muscleReflex.ili.KL                     = muscleReflexParaArray(10);
            muscleReflex.ili.L0                     = muscleReflexParaArray(11);
            muscleReflex.ili_pelvis_tilt.P02        = muscleReflexParaArray(12);
            muscleReflex.ili_pelvis_tilt.KP2        = muscleReflexParaArray(13);
            muscleReflex.ili_pelvis_tilt.KV2        = muscleReflexParaArray(14);
            muscleReflex.ili_ham.KL                 = muscleReflexParaArray(15);
            muscleReflex.ili_ham.L0                 = muscleReflexParaArray(16);
            muscleReflex.ham_pelvis_tilt.KP         = muscleReflexParaArray(17);
            muscleReflex.ham_pelvis_tilt.KV         = muscleReflexParaArray(18);
            muscleReflex.ham_pelvis_tilt.C0         = muscleReflexParaArray(19);
            % muscleReflex.ham.KF                     = muscleReflexParaArray(20);
            muscleReflex.ham_glu.KF                 = muscleReflexParaArray(20);
            muscleReflex.glu_pelvis_tilt.KP         = muscleReflexParaArray(21);
            muscleReflex.glu_pelvis_tilt.KV         = muscleReflexParaArray(22);
            muscleReflex.glu_pelvis_tilt.C0         = muscleReflexParaArray(23);
            muscleReflex.glu.KF                     = muscleReflexParaArray(24);
            % muscleReflex.glu.C0                     = muscleReflexParaArray(25);
            muscleReflex.vas.KF1                    = muscleReflexParaArray(25); % 26 -> 25
            muscleReflex.vas.KF2                    = muscleReflexParaArray(26);
            muscleReflex.vas.C0                     = muscleReflexParaArray(27);
            muscleReflex.vas_knee.pos_max           = muscleReflexParaArray(28);

        end

    end

end
