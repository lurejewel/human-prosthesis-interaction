function D = load_lasso_design_matrices(stoFilePath, options)
%LOAD_LASSO_DESIGN_MATRICES  Build the regression system for LASSO fitting.
%
%   D = load_lasso_design_matrices(stoFilePath, options) reads a LASSO-ready
%   .sto once and returns the 26-feature design matrix, the 9 right-leg
%   excitation targets, the sample weights and the stance/swing row indices
%   (already delay-aligned).  This is the expensive, data-only half of
%   fit_adaptive_lasso_phase_model.m: it is independent of any resampling and
%   is therefore computed once and reused by every replicate of the
%   structure-identification ensemble.
%
%   Input
%   -----
%   stoFilePath : LASSO-ready .sto (column naming as produced by
%                 write_lasso_dataset_sto.m, or any SCONE-style file with
%                 <muscle>_r.activation / .excitation / .fiber_length_norm /
%                 .mtu_force_norm and leg1_r.state).
%   options     : struct, defaults in lasso_fit_default_options.m.  Only
%                 delaySteps, excitationThreshold and sampleWeight are used.
%
%   Output D
%   --------
%   .X            n x p   design matrix (original feature scale)
%   .X_std        n x p   population-standardised X (used by the solver)
%   .muX, .sX     1 x p   standardisation constants (kept for back-transform)
%   .Y            n x nMus excitation targets (right leg, canonical order)
%   .w            n x 1   sample weights already trimmed by delaySteps
%   .phaseIdxMap  .stance / .swing  row indices (column vectors)
%   .phaseMask    .stance / .swing  logical n x 1 masks
%   .muscleNames  1 x nMus cellstr
%   .featureNames p x 1 cellstr (row order of X)
%   .nSamples, .nRawSamples, .p, .nMus, .delaySteps, .stoFile
%
%   See also fit_adaptive_lasso_core, fit_adaptive_lasso_phase_model.

if nargin < 1 || isempty(stoFilePath)
    error('load_lasso_design_matrices:badInput', 'stoFilePath is required.');
end
if nargin < 2 || isempty(options)
    options = struct();
end
opts = lasso_fit_default_options(options);
stoFilePath = char(stoFilePath);
if ~isfile(stoFilePath)
    error('load_lasso_design_matrices:notFound', ...
        'STO file not found: %s', stoFilePath);
end

delaySteps = opts.delaySteps;

%% ===== STEP 1a: READ COLUMN LABELS, DETECT RIGHT-LEG MUSCLES ============
% The old pipeline assumed the 9 right-leg muscles appear in canonical order.
% Here we still fit the right leg only and assume bilateral symmetry (same
% assumption as before).  Detecting via endsWith('_r.activation') keeps both
% the SCONE UN.sto and the Moco interface output working unchanged.
fid = fopen(stoFilePath, 'r');
if fid == -1
    error('load_lasso_design_matrices:cannotOpen', ...
        'Cannot open STO file: %s', stoFilePath);
end
ln = '';
while ischar(ln) && ~contains(ln, 'endheader')
    ln = fgetl(fid);
end
if ~ischar(ln)
    fclose(fid);
    error('load_lasso_design_matrices:noEndHeader', ...
        'endheader not found in STO file: %s', stoFilePath);
end
headerLine = fgetl(fid);
while ischar(headerLine) && isempty(strtrim(headerLine))
    headerLine = fgetl(fid);
end
fclose(fid);
if ~ischar(headerLine)
    error('load_lasso_design_matrices:noHeader', ...
        'No column-label line found in %s.', stoFilePath);
end
colNames = regexp(strtrim(headerLine), '\s+', 'split');
colNames = colNames(~cellfun(@isempty, colNames));

isRightActivation = endsWith(colNames, '_r.activation');
muscleIdx = find(isRightActivation);
if isempty(muscleIdx)
    error('load_lasso_design_matrices:noMuscles', ...
        'No "<muscle>_r.activation" columns found in %s.', stoFilePath);
end
muscleNames = cell(1, numel(muscleIdx));
for i = 1:numel(muscleIdx)
    muscleNames{i} = extractBefore(colNames{muscleIdx(i)}, '.activation');
end
numTargets = numel(muscleNames);
if opts.verbose
    fprintf('[LASSO] %d right-leg muscles detected in %s\n', numTargets, stoFilePath);
    fprintf('        %s\n', strjoin(muscleNames, ', '));
end

%% ===== STEP 1b: SHARED FEATURE MATRIX + TARGETS ==========================
% extractStoMuscleFeatures reads the target excitation, the 9 same-side
% lcEN/fATN columns and the 8 pelvis/joint kinematics columns.  Every muscle
% must share the exact same design matrix; a mismatch means the .sto is
% misaligned, and failing loudly is safer than silently using the wrong
% features.
[yAllBase, XAll, featureToIndexMap, phase] = ...
    extractStoMuscleFeatures(stoFilePath, muscleNames{1});
[yBase, X] = preprocessExcitationAndFeatures(yAllBase, XAll, ...
    delaySteps, opts.excitationThreshold);

Y = zeros(numel(yBase), numTargets);
Y(:, 1) = yBase;
for m = 2:numTargets
    [yAllM, XAllM, featureToIndexMapM, ~] = ...
        extractStoMuscleFeatures(stoFilePath, muscleNames{m});
    if ~isequal(size(XAllM), size(XAll)) || max(abs(XAllM(:) - XAll(:))) > 1e-12
        error('load_lasso_design_matrices:FeatureMismatch', ...
            'Feature matrix for %s is not identical to base muscle %s.', ...
            muscleNames{m}, muscleNames{1});
    end
    if ~isequal(sort(keys(featureToIndexMapM)), sort(keys(featureToIndexMap)))
        error('load_lasso_design_matrices:FeatureMapMismatch', ...
            'Feature map for %s differs from base muscle %s.', ...
            muscleNames{m}, muscleNames{1});
    end
    [yM, XM] = preprocessExcitationAndFeatures(yAllM, XAllM, ...
        delaySteps, opts.excitationThreshold);
    if ~isequal(size(XM), size(X)) || max(abs(XM(:) - X(:))) > 1e-12
        error('load_lasso_design_matrices:PreprocessMismatch', ...
            'Preprocessed feature matrix for %s differs from base muscle %s.', ...
            muscleNames{m}, muscleNames{1});
    end
    Y(:, m) = yM;
end

%% ===== STEP 1c: DESIGN MATRIX + STANDARDISATION ==========================
% The design matrix is X directly (legacy: A = X*M with identity M; removed).
% Standardisation uses the population std; coefficients are transformed back
% to the original feature scale, so this only affects conditioning.
X_orig = X;
muX = mean(X_orig, 1);
sX = std(X_orig, 1, 1);
sX(sX == 0) = 1;
X_std = (X_orig - muX) ./ sX;
p = size(X_std, 2);

nMus = (p - 8) / 2;
if nMus < 1 || abs(nMus - round(nMus)) > 1e-12
    error('load_lasso_design_matrices:InvalidPredictorCount', ...
        'size(X,2) must satisfy p = 2*#mus + 8. Got p = %d.', p);
end
if nMus ~= numTargets
    error('load_lasso_design_matrices:MuscleCountMismatch', ...
        'Design matrix implies %d muscles but %d targets were found.', ...
        nMus, numTargets);
end

%% ===== STEP 2B: MAP PHASE INDICES TO THE DELAYED TIME GRID ===============
% preprocessExcitationAndFeatures shifts the target forward by delaySteps and
% the features backward, so the raw phase indices must be shifted by the same
% amount to align with rows of X/Y/w.  Samples that fall off the front are
% dropped.
nSamples = size(X_std, 1);
rawStanceIdx = phase.stanceIdx(:);
rawSwingIdx = phase.swingIdx(:);

stanceIdx = rawStanceIdx(rawStanceIdx > delaySteps) - delaySteps;
swingIdx = rawSwingIdx(rawSwingIdx > delaySteps) - delaySteps;
stanceIdx = stanceIdx(stanceIdx >= 1 & stanceIdx <= nSamples);
swingIdx = swingIdx(swingIdx >= 1 & swingIdx <= nSamples);

stanceMask = false(nSamples, 1);
stanceMask(stanceIdx) = true;
swingMask = false(nSamples, 1);
swingMask(swingIdx) = true;

if any(stanceMask & swingMask)
    error('load_lasso_design_matrices:PhaseOverlap', ...
        'Stance and swing indices overlap after delay mapping.');
end
if any(~(stanceMask | swingMask))
    warning('fit_adaptive_lasso_phase_model:PhaseGap', ...
        'Phase mapping leaves %d samples unlabeled after delay.', ...
        nnz(~(stanceMask | swingMask)));
end

%% ===== SAMPLE WEIGHTS ===================================================
% The interface stores weights on the full (pre-delay) time grid.  The
% response/features are trimmed by delaySteps, so trim identically.  A
% caller-supplied weight vector already of length nSamples is accepted as is.
wAll = opts.sampleWeight;
if isempty(wAll)
    wAll = ones(nSamples, 1);
else
    wAll = double(wAll(:));
    if numel(wAll) == nSamples + delaySteps
        wAll = wAll((delaySteps + 1):end);
    end
    if numel(wAll) ~= nSamples
        error('load_lasso_design_matrices:BadWeights', ...
            'options.sampleWeight has %d entries, expected %d.', ...
            numel(wAll), nSamples);
    end
    if any(~isfinite(wAll)) || any(wAll < 0)
        error('load_lasso_design_matrices:BadWeights', ...
            'options.sampleWeight must be finite and nonnegative.');
    end
end

%% ===== ASSEMBLE OUTPUT ==================================================
D = struct();
D.stoFile       = stoFilePath;
D.X             = X_orig;
D.X_std         = X_std;
D.muX           = muX;
D.sX            = sX;
D.Y             = Y;
D.w             = wAll;
D.phaseIdxMap   = struct('stance', stanceIdx, 'swing', swingIdx);
D.phaseMask     = struct('stance', stanceMask, 'swing', swingMask);
D.muscleNames   = muscleNames;
D.featureNames  = buildInverseFeatureMap(featureToIndexMap, size(X, 2));
D.nSamples      = nSamples;
D.nRawSamples   = numel(yAllBase);
D.p             = p;
D.nMus          = nMus;
D.delaySteps    = delaySteps;

if opts.verbose
    fprintf('[LASSO] design matrix: %d samples x %d features, %d targets\n', ...
        D.nSamples, D.p, D.nMus);
    fprintf('[LASSO] stance n=%d, swing n=%d, valid weight sum=%.1f\n', ...
        numel(stanceIdx), numel(swingIdx), sum(wAll));
end
end

%% ========================================================================
% Local function: feature name -> row index inverse map.
% ========================================================================
% Strict bijection check: every feature maps to exactly one row, and every row
% is claimed exactly once.
function indexToFeatureName = buildInverseFeatureMap(featureToIndexMap, numFeatures)
ks = keys(featureToIndexMap);
vs = cell2mat(values(featureToIndexMap));
if numel(ks) ~= numFeatures
    error('buildInverseFeatureMap:CountMismatch', ...
        'feature map size (%d) does not match numFeatures (%d).', ...
        numel(ks), numFeatures);
end
if any(vs < 1) || any(vs > numFeatures) || any(abs(vs - round(vs)) > 0)
    error('buildInverseFeatureMap:BadIndex', ...
        'feature indices must be integers in [1, %d].', numFeatures);
end
if numel(unique(vs)) ~= numel(vs)
    error('buildInverseFeatureMap:NonBijective', 'feature map is not one-to-one.');
end
indexToFeatureName = cell(numFeatures, 1);
for t = 1:numel(ks)
    indexToFeatureName{vs(t)} = ks{t};
end
end
