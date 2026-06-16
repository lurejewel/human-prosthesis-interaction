function [initPara, reflexParamMap, reflexTemplate] = ...
    load_lasso_reflex_controller(lassoMatFile)
% Name: load_lasso_reflex_controller
% Description: Load and validate a LASSO-informed sparse linear phase
%   controller from a .mat file.  Returns the flat CMA-ES initial parameter
%   vector, a parameter-layout map, and a controller template struct.
%
%   Supports two LASSO struct formats:
%     - 'grouped' (v2): controllers shared across phase groups; parameters
%        are NOT duplicated in the flat vector.  Fields: .groups(g).beta,
%        .groups(g).bias, .groups(g).mask, .groups(g).phases, .groups(g).label.
%     - 'legacy'  (v1): per-phase beta/bias/mask cell arrays.
%
% Input:
%   lassoMatFile – path to .mat file containing a 'lasso' struct (or a single
%                  struct variable)
% Output:
%   initPara       – N×1 column vector of optimizable parameters
%   reflexParamMap – struct describing the flat-vector layout
%   reflexTemplate – struct with masks, phase IDs, and metadata

%% ---- load and resolve variable ----
vars = whos('-file', lassoMatFile);
if isempty(vars)
    error('LASSO file "%s" is empty.', lassoMatFile);
end

contents = load(lassoMatFile);

if isfield(contents, 'lasso')
    lasso = contents.lasso;
elseif numel(vars) == 1 && isstruct(contents.(vars(1).name))
    lasso = contents.(vars(1).name);
else
    error(['LASSO file must contain a variable named "lasso" or a single ' ...
           'struct variable. Found: %s'], strjoin({vars.name}, ', '));
end

%% ---- detect format ----
if isfield(lasso, 'groups') && ~isempty(lasso.groups)
    useGrouped = true;
    nGroups = numel(lasso.groups);
    nPhases = lasso.nPhases;
else
    useGrouped = false;
    nPhases = validate_and_normalise(lasso);
end

%% ---- build reflexTemplate ----
reflexTemplate = struct();
reflexTemplate.controllerType  = 'lasso_linear_phase';
reflexTemplate.nFeatures       = 22;
reflexTemplate.nMusclesPerSide = 7;
reflexTemplate.nPhases         = nPhases;
reflexTemplate.phaseIds        = lasso.phaseIds(:).';

if useGrouped
    % ---- grouped format: expand group masks to per-phase masks ----
    reflexTemplate.mask = cell(nPhases, 1);
    reflexTemplate.groupOfPhase = zeros(1, nPhases);
    reflexTemplate.nGroups = nGroups;
    reflexTemplate.groupLabels = cell(1, nGroups);
    reflexTemplate.groupPhases = cell(1, nGroups);

    for g = 1:nGroups
        grp = lasso.groups(g);
        reflexTemplate.groupLabels{g} = grp.label;
        reflexTemplate.groupPhases{g} = grp.phases;
        for idx = 1:numel(grp.phases)
            ph = grp.phases(idx);
            phaseIdx = find(lasso.phaseIds == ph, 1);
            if isempty(phaseIdx)
                error('Group %d phase %d not found in lasso.phaseIds.', g, ph);
            end
            reflexTemplate.mask{phaseIdx} = grp.mask;
            reflexTemplate.groupOfPhase(phaseIdx) = g;
        end
    end
else
    % ---- legacy format ----
    if isfield(lasso, 'mask')
        reflexTemplate.mask = lasso.mask(:);  % column cell
    else
        reflexTemplate.mask = cell(nPhases, 1);
        for p = 1:nPhases
            reflexTemplate.mask{p} = abs(lasso.beta{p}) > 0;
        end
    end
end

% optional metadata
if isfield(lasso, 'featureNames')
    reflexTemplate.featureNames = lasso.featureNames;
end
if isfield(lasso, 'muscleNames')
    reflexTemplate.muscleNames = lasso.muscleNames;
end

%% ---- build reflexParamMap ----
reflexParamMap = struct();
reflexParamMap.controllerType  = 'lasso_linear_phase';
reflexParamMap.nFeatures       = 22;
reflexParamMap.nMusclesPerSide = 7;
reflexParamMap.nPhases         = nPhases;
reflexParamMap.phaseIds        = lasso.phaseIds(:).';
reflexParamMap.phases          = [];

initPara = [];

if useGrouped
    % ---- grouped: one set of parameters per group (no duplication) ----
    reflexParamMap.nGroups = nGroups;
    reflexParamMap.groups  = [];

    for g = 1:nGroups
        grp = lasso.groups(g);
        betaMat = grp.beta;         % 22×7
        biasVec = grp.bias(:)';     % 1×7
        maskMat = grp.mask;         % 22×7 logical

        % ---- validate mask consistency ----
        nonzeroOutsideMask = abs(betaMat(~maskMat)) > 0;
        if any(nonzeroOutsideMask)
            error('Group %d (%s): beta has non-zero entries where mask is false.', ...
                g, grp.label);
        end

        betaLinIdx = find(maskMat);
        nBeta = numel(betaLinIdx);

        betaStart = numel(initPara) + 1;
        betaEnd   = betaStart + nBeta - 1;
        biasStart = betaEnd + 1;
        biasEnd   = biasStart + 6;   % 7 bias values

        initPara = [initPara; betaMat(betaLinIdx); biasVec(:)];  %#ok<AGROW>

        reflexParamMap.groups(g).label      = grp.label;
        reflexParamMap.groups(g).phases     = grp.phases;
        reflexParamMap.groups(g).betaLinIdx = betaLinIdx;
        reflexParamMap.groups(g).nBeta      = nBeta;
        reflexParamMap.groups(g).betaStart  = betaStart;
        reflexParamMap.groups(g).betaEnd    = betaEnd;
        reflexParamMap.groups(g).biasStart  = biasStart;
        reflexParamMap.groups(g).biasEnd    = biasEnd;
    end

    % ---- also build per-phase map (for backward compatibility) ----
    for p = 1:nPhases
        g = reflexTemplate.groupOfPhase(p);
        gm = reflexParamMap.groups(g);
        reflexParamMap.phases(p).phaseId    = lasso.phaseIds(p);
        reflexParamMap.phases(p).groupIdx   = g;
        reflexParamMap.phases(p).betaLinIdx = gm.betaLinIdx;
        reflexParamMap.phases(p).nBeta      = gm.nBeta;
        reflexParamMap.phases(p).betaStart  = gm.betaStart;
        reflexParamMap.phases(p).betaEnd    = gm.betaEnd;
        reflexParamMap.phases(p).biasStart  = gm.biasStart;
        reflexParamMap.phases(p).biasEnd    = gm.biasEnd;
    end
else
    % ---- legacy: per-phase layout ----
    for p = 1:nPhases
        betaMat  = lasso.beta{p};   % 22×7
        biasVec  = lasso.bias{p}(:)';  % 1×7
        maskMat  = reflexTemplate.mask{p};

        % ---- validate mask consistency ----
        nonzeroOutsideMask = abs(betaMat(~maskMat)) > 0;
        if any(nonzeroOutsideMask)
            error('Phase %d: beta has non-zero entries where mask is false.', p);
        end

        % ---- column-major linear indices of sparse beta entries ----
        betaLinIdx = find(maskMat);   % column-major into 22×7
        nBeta      = numel(betaLinIdx);

        betaStart = numel(initPara) + 1;
        betaEnd   = betaStart + nBeta - 1;
        biasStart = betaEnd + 1;
        biasEnd   = biasStart + 6;   % 7 bias values

        % ---- append to flat vector ----
        initPara = [initPara; betaMat(betaLinIdx); biasVec(:)];  %#ok<AGROW>

        % ---- store per-phase map ----
        reflexParamMap.phases(p).phaseId   = lasso.phaseIds(p);
        reflexParamMap.phases(p).betaLinIdx = betaLinIdx;
        reflexParamMap.phases(p).nBeta     = nBeta;
        reflexParamMap.phases(p).betaStart = betaStart;
        reflexParamMap.phases(p).betaEnd   = betaEnd;
        reflexParamMap.phases(p).biasStart = biasStart;
        reflexParamMap.phases(p).biasEnd   = biasEnd;
    end
end

reflexParamMap.totalLen = numel(initPara);
initPara = initPara(:);  % column vector

%% ---- print summary (once, not inside optimisation) ----
if useGrouped
    fprintf('[LASSO controller] %d groups → %d phases, %d total optimizable parameters.\n', ...
        nGroups, nPhases, reflexParamMap.totalLen);
    for g = 1:nGroups
        gm = reflexParamMap.groups(g);
        fprintf('  Group %d (%s): phases [%s], %d non-zero beta entries, 7 bias entries\n', ...
            g, gm.label, strjoin(string(gm.phases), ', '), gm.nBeta);
    end
else
    fprintf('[LASSO controller] %d phases loaded, %d total optimizable parameters.\n', ...
        nPhases, reflexParamMap.totalLen);
    for p = 1:nPhases
        fprintf('  Phase %d: %d non-zero beta entries, 7 bias entries\n', ...
            lasso.phaseIds(p), reflexParamMap.phases(p).nBeta);
    end
end

end

% ========================================================================
function nPhases = validate_and_normalise(lasso)
%% ---- phase count ----
if isfield(lasso, 'nPhases')
    nPhases = lasso.nPhases;
else
    nPhases = numel(lasso.beta);
end

if ~isfield(lasso, 'beta')
    error('LASSO struct must contain a "beta" field.');
end
if ~isfield(lasso, 'bias') && ~isfield(lasso, 'b')
    error('LASSO struct must contain a "bias" (or "b") field.');
end

% ---- normalise to column cell arrays ----
lasso.beta = lasso.beta(:);
if isfield(lasso, 'bias')
    lasso.bias = lasso.bias(:);
else
    lasso.bias = lasso.b(:);
    lasso.bias = lasso.bias(:);
end

if isfield(lasso, 'mask')
    lasso.mask = lasso.mask(:);
end

%% ---- phaseIds ----
if isfield(lasso, 'phaseIds')
    lasso.phaseIds = lasso.phaseIds(:)';
else
    lasso.phaseIds = 0:(nPhases-1);
end

%% ---- validation ----
assert(numel(lasso.beta) == nPhases, ...
    'nPhases (%d) != numel(beta) (%d).', nPhases, numel(lasso.beta));
assert(numel(lasso.bias) == nPhases, ...
    'nPhases (%d) != numel(bias) (%d).', nPhases, numel(lasso.bias));

for p = 1:nPhases
    [br, bc] = size(lasso.beta{p});
    assert(br == 22 && bc == 7, ...
        'Phase %d: beta must be 22×7, got %d×%d.', p, br, bc);
    biasVec = lasso.bias{p}(:)';
    assert(numel(biasVec) == 7, ...
        'Phase %d: bias must have 7 elements, got %d.', p, numel(biasVec));
    lasso.bias{p} = biasVec;

    if isfield(lasso, 'mask')
        [mr, mc] = size(lasso.mask{p});
        assert(mr == 22 && mc == 7, ...
            'Phase %d: mask must be 22×7, got %d×%d.', p, mr, mc);
        assert(islogical(lasso.mask{p}), ...
            'Phase %d: mask must be logical.', p);
    end
end
end
