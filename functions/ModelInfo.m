% -------------------------------------------------------------------------
% Name: ModelInfo.m
% Description: definition of modelInfo Class.
% Member variables:
%
% Member functions:
% - ModelInfo(): construction function.
% - reset_record(): reset the class to its original state (no record except
% the initial states)
%
% NOTE:
% org.opensim.modeling.model and org.opensim.modeling.state should NOT
% appear in this class.
% -------------------------------------------------------------------------

classdef ModelInfo < handle

    properties

        % static properties; initialized once before optimization
        st
        % ├── muscle                    STATIC MUSCLE PROPERTIES
        % │   ├── names
        % │   ├── lambda                    muscle type I fiber composition
        % │   ├── mass
        % │   ├── lopt                      optimal fiber length
        % │   ├── fopt                      optimal fiber force
        % │   └── map                       abbr -> index 目前只支持lambda取值，其他数组需要改写
        % │   └── delay                     time delay of electrical signals from the central nervous system to the muscle fibers
        % ├── model
        % │   ├── totalMass
        % │   ├── initPose
        % │   └── map
        % └── simInfo
        %     ├── stepTime
        %     ├── endTime
        %     ├── speed
        %     ├── timeSeries
        %     ├── [DEPRECATED] npts
        %     └── [DEPRECATED] gmax                      max #generation
        

        % dynamic properties; updated at every forward simulation
        dy
        % ├── muscle                    DYNAMIC MUSCLE PROPERTIES
        % │   ├── para                      structral parameters for the muscle-reflex controller
        % │   ├── arx                       para in array form
        % │   ├── arz
        % │   ├── exc                       muscle excitations (u)
        % │   ├── act                       muscle activations (a)
        % │   ├── fMTU                      muscle-tendon unit forces, needed for calculation muscle effort
        % │   ├── fCE                       muscle fiber forces, needed for calculation muscle effort
        % │   ├── fATN                      muscle fiber forces along tendon, normalized, needed for calculation muscle excitation
        % │   ├── lCEN                      muscle fiber lengths, normalized
        % │   └── vCE                       velocity of lCEN
        % ├── grf                       GROUND REACTION FORCES
        % │   ├── fyr                       normal reaction force for the right leg
        % │   ├── fxr                       friction force for the right leg
        % │   ├── fyl                       normal reaction force for the left leg
        % │   ├── fxl                       friction force for the left leg
        % ├── phase                     GAIT PHASES
        % │   ├── r
        % │   └── l
        % ├── stateHistory
        % └── lastTime                      the instant when model falls down; =endTime if the model can walk stably

        % % model information
        % mass % total mass of the model
        % 
        % % muscle-related information
        % muscleNames
        % muscleMap % map: abbreviation -> index 目前只支持lambda的取值，其他数组需要改写
        % muscleExcitations % u
        % muscleActivations % a
        % muscleForces % 其实是fiber force
        % muscleActiveForces % F^{CE}
        % muscleFiberForcesATN % AT: along tendon; N: normalized。这个才是真正的muscle force(归一化后）
        % muscleFiberLengthN % \tilde l^{CE}, or say lCEN
        % muscleFiberVelocity % v^{CE}
        % muscleReflexDelay % time delay of electrical signals from the central nervous system to the muscle fibers
        % muscleReflexParam % structure of muscle reflex parameters
        % muscleLambda % type I ratio
        % muscleMass
        % arx % array of muscle reflex parameters
        % arz
        % OptimalFiberLengths
        % OptimalFiberForces
        % 
        % % gait phases
        % phaseR
        % phaseL
        % 
        % % others
        % grf
        % 
        % % state
        % stateHistory
        % stateMap
        % timeSeries
        % lastTime % the instant when model falls down; equals the simConfig.endTime if the model can walk stably
        % % stage % =1/2 for the 2-stage optimization
        % g % number of generations
        % gMax % maximum number of generation        
    end

    methods
        function obj = ModelInfo(modelStaticProp)
            % Class Name: ModelInfo
            % Description: init static (obj.st.*) and dynamic (obj.dy.* ) parameters

            obj.st = modelStaticProp;

            nMuscles = numel(obj.st.muscle.names);
            npts = ceil(obj.st.simInfo.endTime / obj.st.simInfo.stepTime)+1;
            nStates = obj.st.model.map.Count;
            obj.dy.muscle.exc = nan(nMuscles, npts+obj.st.muscle.delay);
            obj.dy.muscle.act = nan(nMuscles, npts);
            obj.dy.muscle.fMTU = nan(nMuscles, npts);
            obj.dy.muscle.fCE = nan(nMuscles, npts);
            obj.dy.muscle.fATN = nan(nMuscles, npts);
            obj.dy.muscle.lCEN = nan(nMuscles, npts);
            obj.dy.muscle.vCE = nan(nMuscles, npts);
            obj.dy.stateHistory = nan(nStates, npts);
            obj.dy.grf.fyr = nan(1, npts);
            obj.dy.grf.fxr = nan(1, npts);
            obj.dy.grf.fyl = nan(1, npts);
            obj.dy.grf.fxl = nan(1, npts);
            obj.dy.phase.r = nan(1, npts);
            obj.dy.phase.l = nan(1, npts);
            obj.dy.lastTime = -1;

            % % number of points
            % npts = ceil(simConfig.endTime / simConfig.stepTime)+1;
            % 
            % % muslce states: excitation, activation, muscle fiber forces along tendon(or not?)
            % delay = round(0.01/simConfig.stepTime); % 10 ms delay for muscle reflex mechanism to take effect
            % muscleNum = model.getMuscles.getSize; % number of muscles
            % [obj.muscleLambda, obj.muscleMass, obj.muscleMap] = read_muscle_lambda_and_mass(['muscle_info_' num2str(muscleNum) '.csv']); 
            % obj.muscleNames = cell(muscleNum, 1);
            % obj.muscleReflexDelay = delay;
            % obj.muscleExcitations = nan(muscleNum, npts+delay);
            % obj.muscleActivations = nan(muscleNum, npts);
            % obj.muscleForces = nan(muscleNum, npts);
            % obj.muscleActiveForces = nan(muscleNum, npts);
            % obj.muscleFiberForcesATN = nan(muscleNum, npts);
            % obj.muscleFiberLengthN = nan(muscleNum, npts);
            % obj.muscleFiberVelocity = nan(muscleNum, npts);
            % 
            % % muscle names
            % for muscleIndex = 0 : muscleNum-1
            %     obj.muscleNames{muscleIndex+1} = char(model.getMuscles.get(muscleIndex).getName);
            % end
            % 
            % % muscle reflex parameters
            % obj.read_muscleReflex_array(initPara);
            % % optimal fiber lengths后续可以改成数组+muscleNameMap的形式
            % obj.OptimalFiberLengths.hamstrings = model.getMuscles.get('hamstrings_r').getOptimalFiberLength;
            % obj.OptimalFiberLengths.glut_max = model.getMuscles.get('glut_max_r').getOptimalFiberLength;
            % obj.OptimalFiberLengths.iliopsoas = model.getMuscles.get('iliopsoas_r').getOptimalFiberLength;
            % obj.OptimalFiberLengths.vasti = model.getMuscles.get('vasti_r').getOptimalFiberLength;
            % obj.OptimalFiberLengths.gastrocnemius = model.getMuscles.get('gastroc_r').getOptimalFiberLength;
            % obj.OptimalFiberLengths.soleus = model.getMuscles.get('soleus_r').getOptimalFiberLength;
            % obj.OptimalFiberLengths.tibia = model.getMuscles.get('tibia_r').getOptimalFiberLength;
            % 
            % % optimal fiber forces
            % obj.OptimalFiberForces.hamstrings = model.getMuscles.get('hamstrings_r').getMaxIsometricForce;
            % obj.OptimalFiberForces.glut_max = model.getMuscles.get('glut_max_r').getMaxIsometricForce;
            % obj.OptimalFiberForces.iliopsoas = model.getMuscles.get('iliopsoas_r').getMaxIsometricForce;
            % obj.OptimalFiberForces.vasti = model.getMuscles.get('vasti_r').getMaxIsometricForce;
            % obj.OptimalFiberForces.gastrocnemius = model.getMuscles.get('gastroc_r').getMaxIsometricForce;
            % obj.OptimalFiberForces.soleus = model.getMuscles.get('soleus_r').getMaxIsometricForce;
            % obj.OptimalFiberForces.tibia = model.getMuscles.get('tibia_r').getMaxIsometricForce;
            % 
            % % phases
            % obj.phaseR = nan(1, npts);
            % obj.phaseL = nan(1, npts);
            % 
            % % ground reaction forces
            % obj.grf.normalR = nan(1, npts);
            % obj.grf.normalL = nan(1, npts);
            % obj.grf.frictionR = nan(1, npts);
            % obj.grf.frictionL = nan(1, npts);
            % 
            % % state history
            % state = model.initSystem(); % this is a must before calling model.getNumStateVariables()
            % obj.stateHistory = nan(model.getNumStateVariables, npts);
            % keys = [];
            % for i = 0 : state.getNQ-1 % value of jointset
            %     keys = [keys, string([char(model.getCoordinateSet.get(i)), '/value'])];
            % end
            % for i = 0 : state.getNU-1 % speed of jointset
            %     keys = [keys, string([char(model.getCoordinateSet.get(i)), '/speed'])];
            % end
            % for i = state.getNQ+state.getNU : state.getNQ+state.getNU+state.getNZ-1 % value of forceset
            %     keys = [keys, string(model.getStateVariableNames.get(i))];
            % end
            % obj.stateMap = containers.Map(keys, 1:length(keys));            
            % 
            % % time series, last time, stage
            % obj.timeSeries = 0 : simConfig.stepTime : simConfig.endTime;
            % obj.lastTime = -1;
            % % obj.stage = 1;
            % obj.g = 1;
            % obj.gMax = 1000;
            % 
            % % basic model information
            % obj.mass = model.getTotalMass(state);
            
        end

        function read_muscleReflex_array(obj, arx, arz)
            % Description: convert the muscle-reflex parameters from array 
            % form to struct form.
            % 后面可以考虑用数组+containers.map代替，方便debug和调优（该函数后续删去）
            obj.dy.muscle.arx = arx;
            obj.dy.muscle.arz = arz;

            muscleReflex.tib.KL                     = arx(1);
            muscleReflex.tib.L0                     = arx(2);
            muscleReflex.tib_sol.KF                 = arx(3);
            muscleReflex.sol.KF                     = arx(4);
            muscleReflex.gas.KF                     = arx(5);
            muscleReflex.ili_pelvis_tilt.KP         = arx(6);
            muscleReflex.ili_pelvis_tilt.KV         = arx(7);
            muscleReflex.ili_pelvis_tilt.C0         = arx(8);
            muscleReflex.ili.C0                     = arx(9);
            muscleReflex.ili.KL                     = arx(10);
            muscleReflex.ili.L0                     = arx(11);
            muscleReflex.ili_pelvis_tilt.P02        = arx(12);
            muscleReflex.ili_pelvis_tilt.KP2        = arx(13);
            muscleReflex.ili_pelvis_tilt.KV2        = arx(14);
            muscleReflex.ili_ham.KL                 = arx(15);
            muscleReflex.ili_ham.L0                 = arx(16);
            muscleReflex.ham_pelvis_tilt.KP         = arx(17);
            muscleReflex.ham_pelvis_tilt.KV         = arx(18);
            muscleReflex.ham_pelvis_tilt.C0         = arx(19);
            % muscleReflex.ham.KF                     = muscleReflexParaArray(20);
            muscleReflex.ham_glu.KF                 = arx(20);
            muscleReflex.glu_pelvis_tilt.KP         = arx(21);
            muscleReflex.glu_pelvis_tilt.KV         = arx(22);
            muscleReflex.glu_pelvis_tilt.C0         = arx(23);
            muscleReflex.glu.KF                     = arx(24);
            % muscleReflex.glu.C0                     = muscleReflexParaArray(25);
            muscleReflex.vas.KF1                    = arx(25); % 26 -> 25
            muscleReflex.vas.KF2                    = arx(26);
            muscleReflex.vas.C0                     = arx(27);
            muscleReflex.vas_knee.pos_max           = arx(28);

            % obj.muscleReflexParam = muscleReflex;
            obj.dy.muscle.para = muscleReflex;

        end

        function reset_record(obj)
            
            obj.dy.muscle.exc(:) = nan;
            obj.dy.muscle.act(:) = nan;
            obj.dy.muscle.fMTU(:) = nan;
            obj.dy.muscle.fCE(:) = nan;
            obj.dy.muscle.fATN(:,2:end) = nan;
            obj.dy.muscle.lCEN(:,2:end) = nan;
            obj.dy.muscle.vCE(:) = nan;
            obj.dy.stateHistory(:) = nan;
            obj.dy.grf.fyr(:) = nan;
            obj.dy.grf.fxr(:) = nan;
            obj.dy.grf.fyl(:) = nan;
            obj.dy.grf.fxl(:) = nan;
            obj.dy.phase.r(:) = nan;
            obj.dy.phase.l(:) = nan;
            obj.dy.lastTime = -1;
            % obj.muscleExcitations(:) = nan;
            % obj.muscleActivations(:) = nan;
            % obj.muscleForces(:) = nan;
            % obj.muscleActiveForces(:) = nan;
            % obj.muscleFiberVelocity(:) = nan;
            % obj.stateHistory(:) = nan;
            % obj.grf.normalR(:) = nan;
            % obj.grf.normalL(:) = nan;
            % obj.grf.frictionR(:) = nan;
            % obj.grf.frictionL(:) = nan;
            % obj.phaseR(:) = nan;
            % obj.phaseL(:) = nan;
            % 
            % f = obj.muscleFiberForcesATN(:,1);
            % obj.muscleFiberForcesATN(:) = nan;
            % obj.muscleFiberForcesATN(:,1) = f;
            % 
            % l = obj.muscleFiberLengthN(:,1);
            % obj.muscleFiberLengthN(:) = nan;
            % obj.muscleFiberLengthN(:,1) = l;

        end


    end

end
