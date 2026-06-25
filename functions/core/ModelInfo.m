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
        % │   ├── [DEPRECATED] act          muscle activations (a)
        % │   ├── [DEPRECATED] fMTU         muscle-tendon unit forces
        % │   ├── [DEPRECATED] fCE          muscle fiber forces
        % │   ├── fATN                      muscle fiber forces along tendon, normalized, needed for calculation muscle excitation
        % │   ├── lCEN                      muscle fiber lengths, normalized
        % │   └── [DEPRECATED] vCE          velocity of lCEN
        % ├── grf                       GROUND REACTION FORCES
        % │   ├── fyr                       normal reaction force for the right leg
        % │   ├── fxr                       friction force for the right leg
        % │   ├── fyl                       normal reaction force for the left leg
        % │   ├── fxl                       friction force for the left leg
        % ├── limitForce                COORDINATE LIMIT FORCES
        % │   ├── kneeR                     knee joint limit force, right
        % │   └── kneeL                     knee joint limit force, left
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
            % obj.dy.muscle.act = nan(nMuscles, npts);    % DEPRECATED
            % obj.dy.muscle.fMTU = nan(nMuscles, npts);    % DEPRECATED
            % obj.dy.muscle.fCE = nan(nMuscles, npts);    % DEPRECATED
            obj.dy.muscle.fATN = nan(nMuscles, npts);
            obj.dy.muscle.lCEN = nan(nMuscles, npts);
            % obj.dy.muscle.vCE = nan(nMuscles, npts);    % DEPRECATED
            obj.dy.stateHistory = nan(nStates, npts);
            obj.dy.grf.fyr = nan(1, npts);
            obj.dy.grf.fxr = nan(1, npts);
            obj.dy.grf.fyl = nan(1, npts);
            obj.dy.grf.fxl = nan(1, npts);
            obj.dy.limitForce.kneeR = nan(1, npts);
            obj.dy.limitForce.kneeL = nan(1, npts);
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
            % obj.dy.muscle.act(:) = nan;    % DEPRECATED
            % obj.dy.muscle.fMTU(:) = nan;    % DEPRECATED
            % obj.dy.muscle.fCE(:) = nan;    % DEPRECATED
            obj.dy.muscle.fATN(:,2:end) = nan;
            obj.dy.muscle.lCEN(:,2:end) = nan;
            % obj.dy.muscle.vCE(:) = nan;    % DEPRECATED
            obj.dy.stateHistory(:) = nan;
            obj.dy.grf.fyr(:) = nan;
            obj.dy.grf.fxr(:) = nan;
            obj.dy.grf.fyl(:) = nan;
            obj.dy.grf.fxl(:) = nan;
            obj.dy.limitForce.kneeR(:) = nan;
            obj.dy.limitForce.kneeL(:) = nan;
            obj.dy.phase.r(:) = nan;
            obj.dy.phase.l(:) = nan;
            obj.dy.lastTime = -1;
            obj.dy.hasRun  = false;

        end

        function update(obj, state, allMuscles, kneeLimitForceR, kneeLimitForceL, idx)
            % Name: update
            % Description: update the muscle-related information in the Class modelInfo
            %   for the muscle excitation in the next loop. The recorded information
            %   includes:
            %   - obj.dy.muscle.fATN: normalized fiber force along tendon
            %   - obj.dy.muscle.lCEN: normalized fiber length
            %   - obj.dy.stateHistory: state of the model
            %   - obj.dy.limitForce.kneeR / kneeL: knee joint coordinate limit forces

            fopt = obj.st.muscle.fopt;
            map = obj.st.model.map;
            Y = state.getY.getAsMat();

            for i = 1 : numel(fopt)
                muscle = allMuscles{i};
                obj.dy.muscle.fATN(i, idx+1) = muscle.getFiberForceAlongTendon(state) / fopt(i);
                obj.dy.muscle.lCEN(i, idx+1) = muscle.getNormalizedFiberLength(state);
                % obj.dy.muscle.vCE(i, idx+1) = muscle.getFiberVelocity(state);
                % obj.dy.muscle.act(i, idx+1) = muscle.getActivation(state);
                % obj.dy.muscle.fMTU(i, idx+1) = muscle.getFiberForce(state);
                % obj.dy.muscle.fCE(i, idx+1) = muscle.getActiveFiberForce(state);

                key = ['/forceset/' char(allMuscles{i}.getName) '/activation'];
                Y(map(key)) = muscle.getActivation(state);
            end

            obj.dy.stateHistory(:, idx) = Y;
            obj.dy.limitForce.kneeR(idx) = kneeLimitForceR.getRecordValues(state).get(0);
            obj.dy.limitForce.kneeL(idx) = kneeLimitForceL.getRecordValues(state).get(0);
        end

    end

end