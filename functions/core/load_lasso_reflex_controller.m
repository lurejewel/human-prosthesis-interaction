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
%     - 'legacy'  (v1): per-phase beta/bias/mask cell arrays.  A `b` field is
%       accepted as an alias for `bias`, and a missing `phaseIds` defaults to
%       0:(nPhases-1).  The normalised struct is what the rest of this
%       function consumes (the old helper normalised a local copy and its
%       changes were lost).
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
elseif isscalar(vars) && isstruct(contents.(vars(1).name))
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
    [nPhases, lasso] = validate_and_normalise(lasso);
end

%% ---- grouped format: fail early with a clear message ----
if useGrouped
    requiredGroupFields = {'label', 'phases', 'beta', 'bias', 'mask'};
    for g = 1:nGroups
        grp = lasso.groups(g);
        if ~isstruct(grp)
            error('lasso.groups(%d) must be a struct.', g);
        end
        for f = 1:numel(requiredGroupFields)
            if ~isfield(grp, requiredGroupFields{f})
                error('Group %d is missing required field "%s".', ...
                    g, requiredGroupFields{f});
            end
        end
    end
end

%% ---- build reflexTemplate ----
reflexTemplate = struct();
reflexTemplate.controllerType  = 'lasso_linear_phase';
reflexTemplate.nFeatures       = 26;
reflexTemplate.nMusclesPerSide = 9;
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
reflexParamMap.nFeatures       = 26;
reflexParamMap.nMusclesPerSide = 9;
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

        % ---- validate shapes with clear error messages ----
        betaMat = grp.beta;
        if ~isnumeric(betaMat) || ~isequal(size(betaMat), [26, 9])
            error('Group %d (%s): beta must be a 26×9 numeric matrix, got %s.', ...
                g, grp.label, mat2str(size(betaMat)));
        end
        biasVec = grp.bias;
        if ~isnumeric(biasVec) || numel(biasVec) ~= 9
            error('Group %d (%s): bias must be numeric with 9 elements, got %s.', ...
                g, grp.label, mat2str(size(biasVec)));
        end
        biasVec = biasVec(:)';     % 1×9
        maskMat = grp.mask;
        if ~islogical(maskMat) || ~isequal(size(maskMat), [26, 9])
            error('Group %d (%s): mask must be a 26×9 logical matrix, got %s (%s).', ...
                g, grp.label, mat2str(size(maskMat)), class(maskMat));
        end

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
        biasEnd   = biasStart + 8;   % 9 bias values

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
        betaMat  = lasso.beta{p};   % 26×9
        biasVec  = lasso.bias{p}(:)';  % 1×9
        maskMat  = reflexTemplate.mask{p};

        % ---- validate mask consistency ----
        nonzeroOutsideMask = abs(betaMat(~maskMat)) > 0;
        if any(nonzeroOutsideMask)
            error('Phase %d: beta has non-zero entries where mask is false.', p);
        end

        % ---- column-major linear indices of sparse beta entries ----
        betaLinIdx = find(maskMat);   % column-major into 26×9
        nBeta      = numel(betaLinIdx);

        betaStart = numel(initPara) + 1;
        betaEnd   = betaStart + nBeta - 1;
        biasStart = betaEnd + 1;
        biasEnd   = biasStart + 8;   % 9 bias values

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
        fprintf('  Group %d (%s): phases [%s], %d non-zero beta entries, 9 bias entries\n', ...
            g, gm.label, strjoin(string(gm.phases), ', '), gm.nBeta);
    end
else
    fprintf('[LASSO controller] %d phases loaded, %d total optimizable parameters.\n', ...
        nPhases, reflexParamMap.totalLen);
    for p = 1:nPhases
        fprintf('  Phase %d: %d non-zero beta entries, 9 bias entries\n', ...
            lasso.phaseIds(p), reflexParamMap.phases(p).nBeta);
    end
end

end

% ========================================================================
function [nPhases, lasso] = validate_and_normalise(lasso)
% Normalise a legacy (v1) per-phase LASSO struct and validate it.
% Returns the normalised struct as a second output: MATLAB passes structs by
% value, so normalising a local copy without returning it left the caller
% with the raw struct (missing phaseIds, `b` alias, column-cell layouts, ...).
%
% Normalisation:
%   - beta, bias (and mask, when present) become column cell arrays;
%   - `b` is accepted as an alias for `bias` and is copied to lasso.bias;
%   - missing phaseIds defaults to 0:(nPhases-1);
%   - each per-phase bias is stored as a 1x9 row vector.
% Validation:
%   - beta must be 26x9 numeric and bias must have 9 elements per phase;
%   - mask, when present, must be 26x9 logical per phase.

%% ---- required fields ----
if ~isfield(lasso, 'beta')
    error('Legacy LASSO struct must contain a "beta" field.');
end
if ~isfield(lasso, 'bias') && ~isfield(lasso, 'b')
    error('Legacy LASSO struct must contain a "bias" (or "b") field.');
end
if ~iscell(lasso.beta)
    error(['Legacy LASSO struct: "beta" must be a cell array of per-phase ' ...
           '26×9 matrices.']);
end
if isfield(lasso, 'bias')
    if ~iscell(lasso.bias)
        error(['Legacy LASSO struct: "bias" must be a cell array of ' ...
               'per-phase 9-element vectors.']);
    end
else
    if ~iscell(lasso.b)
        error(['Legacy LASSO struct: "b" (alias for "bias") must be a cell ' ...
               'array of per-phase 9-element vectors.']);
    end
end
if isfield(lasso, 'mask') && ~iscell(lasso.mask)
    error(['Legacy LASSO struct: "mask" must be a cell array of per-phase ' ...
           '26×9 logical matrices.']);
end

%% ---- phase count ----
if isfield(lasso, 'nPhases')
    nPhases = lasso.nPhases;
else
    nPhases = numel(lasso.beta);
end

% ---- normalise to column cell arrays ----
lasso.beta = lasso.beta(:);
if isfield(lasso, 'bias')
    lasso.bias = lasso.bias(:);
else
    % legacy v1 alias: `b` == `bias`
    lasso.bias = lasso.b(:);
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
    assert(br == 26 && bc == 9, ...
        'Phase %d: beta must be 26×9, got %d×%d.', p, br, bc);
    biasVec = lasso.bias{p}(:)';
    assert(numel(biasVec) == 9, ...
        'Phase %d: bias must have 9 elements, got %d.', p, numel(biasVec));
    lasso.bias{p} = biasVec;

    if isfield(lasso, 'mask')
        [mr, mc] = size(lasso.mask{p});
        assert(mr == 26 && mc == 9, ...
            'Phase %d: mask must be 26×9, got %d×%d.', p, mr, mc);
        assert(islogical(lasso.mask{p}), ...
            'Phase %d: mask must be logical.', p);
    end
end
end
