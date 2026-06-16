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

        % LASSO controller properties (set once per optimisation run)
        reflexParamMap   % struct from load_lasso_reflex_controller
        reflexTemplate   % struct with masks, phaseIds, metadata
        reflexParams     % struct with unpacked beta/bias cells

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
            obj.dy.hasRun  = false;

        end

        function read_muscleReflex_array(obj, arx)
            % Description: convert the muscle-reflex parameters from array
            %   form to struct form.  Supports both the legacy hand-crafted
            %   controller and the new LASSO linear-phase controller.
            %   - LASSO path: uses obj.reflexParamMap + obj.reflexTemplate
            %   - Legacy path: uses muscle_reflex_param_defs()
            % obj.dy.muscle.arx = arx;
            % obj.dy.muscle.arz = arz;

            if ~isempty(obj.reflexParamMap) && ~isempty(obj.reflexTemplate)
                % ---- LASSO linear-phase controller ----
                obj.reflexParams = unpack_lasso_reflex_params( ...
                    arx, obj.reflexParamMap, obj.reflexTemplate);
            else
                error('Legacy hand-crafted code commented.');
                % % ---- legacy hand-crafted controller ----
                % defs = muscle_reflex_param_defs();
                % muscleReflex = struct();
                % for i = 1 : numel(defs)
                %     d = defs(i);
                %     muscleReflex.(d.path{1}).(d.path{2}) = arx(d.idx);
                % end
                % obj.dy.muscle.para = muscleReflex;
            end
        end

        % end

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
            obj.dy.hasRun  = false;

        end

    end

end
