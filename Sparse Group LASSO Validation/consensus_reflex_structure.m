function C = consensus_reflex_structure(masks, opts)
%CONSENSUS_REFLEX_STRUCTURE  Consensus 0-1 structure from K replicate masks.
%
%   C = consensus_reflex_structure(masks, opts)
%
%   `masks` is a K x nGroups cell array; masks{k,g} is the p x nMus logical
%   structure (26 x 9 for the current 9-muscle controller) recovered by
%   replicate k for phase group g.  This function turns the ensemble into:
%
%     * Pi{g}          per-entry selection frequency (the "confidence map"),
%     * the modal (most frequent) exact mask and its frequency,
%     * the medoid replicate (highest mean Jaccard to the others),
%     * a Jaccard similarity matrix + single-linkage clusters at
%       opts.clusterJaccard, and the dominant-cluster fraction,
%     * a consensus mask selected by opts.rule.
%
%   Selection rules
%   ---------------
%   'majority'  (default) C = (Pi >= 0.5), ties at exactly 0.5 resolved by
%               opts.tieBreak ('include' | 'exclude').  This is the
%               coordinate-wise median of the masks, i.e. the mask with the
%               smallest total Hamming distance to all replicates.
%               Property: if m of the K replicates produce the *identical*
%               mask M, then Pi equals 1 or 0 on every entry as long as
%               m/K > 0.5, so the majority vote returns M exactly.  The
%               "13 out of 20 seeds agree -> take that framework" rule is
%               therefore a special case, not extra logic.
%   'mode'      the most frequent exact mask (ties broken by mean Jaccard).
%   'medoid'    the mask of the replicate with the highest mean Jaccard.
%   'stability' C = (Pi >= opts.piThreshold), a Meinshausen-Buhlmann style
%               stricter threshold (default 0.7).
%
%   opts
%   ----
%   rule             'majority' (default)
%   piThreshold      0.7   threshold used by the 'stability' rule
%   tieBreak         'include' (default) | 'exclude' for Pi == 0.5 exactly
%   stableOnPi       0.9   diagnostic: entries with Pi >= this are "stable on"
%   stableOffPi      0.1   diagnostic: entries with Pi <= this are "stable off"
%   clusterJaccard   0.9   edge threshold of the mask-similarity graph
%   conformJaccard   0.8   a replicate "supports" the consensus when its
%                          Jaccard similarity to the consensus mask is >= this
%   valid            K x nGroups logical (default all true); invalid
%                    replicates are excluded from that group's statistics.
%                    Invalid masks may be empty (e.g. failed replicates are
%                    recorded as []) -- they are ignored, and only valid
%                    masks are size-checked against each other.
%
%   Output C
%   --------
%   .rule, .piThreshold, .tieBreak, .clusterJaccard, .nSeeds
%   .mask{g}                   p x nMus logical consensus structure
%   .Pi{g}                     p x nMus selection frequency
%   .nValid(g)                 number of replicates used for group g
%   .modeMask{g}, .modeCount(g), .modeFraction(g), .modeTie(g)
%   .medoidIdx(g), .medoidMask{g}
%   .jaccardMean(g), .jaccardMatrix{g}
%   .jaccardToConsensus{g}     Jaccard of each valid replicate to C.mask{g}
%   .coverage{g}, .meanCoverage(g)
%                              share of each replicate's structure retained by
%                              the consensus (for a majority-vote superset this
%                              is the meaningful agreement number)
%   .nConforming(g), .conformingFraction(g)
%                              replicates within conformJaccard of the chosen
%                              framework (strict "same framework" count)
%   .clusters{g}, .clusterSizes{g}, .dominantCluster(g),
%   .dominantClusterFraction(g)
%   .nStableOn(g), .nStableOff(g), .nAmbiguous(g)
%   .supportSizes{g}           nnz of each valid replicate mask
%   .exactMaskCounts{g}        frequency spectrum of the distinct exact masks
%                              (descending), i.e. "how many seeds share what"
%   .nGroups, .nElem
%
%   See also lasso_stability_identification, bootstrap_phase_runs.

if nargin < 2 || isempty(opts)
    opts = struct();
end
opts = fill_defaults(opts);

if ~iscell(masks) || isempty(masks)
    error('consensus_reflex_structure:BadInput', 'masks must be a non-empty cell array.');
end
K = size(masks, 1);
nGroups = size(masks, 2);
if nGroups < 1
    error('consensus_reflex_structure:BadInput', 'masks must have at least one column.');
end

% Validate opts.valid BEFORE the mask-size check: failed replicates are
% recorded with empty masks and must be ignored here.  Only valid masks are
% size-checked; the reference size comes from the first valid, non-empty mask.
if isfield(opts, 'valid') && ~isempty(opts.valid)
    valid = logical(opts.valid);
    if ~isequal(size(valid), [K, nGroups])
        error('consensus_reflex_structure:BadValid', ...
            'opts.valid must be %d x %d.', K, nGroups);
    end
else
    valid = true(K, nGroups);
end

refSize = [];
for k = 1:K
    for g = 1:nGroups
        if valid(k, g) && ~isempty(masks{k, g})
            refSize = size(masks{k, g});
            break;
        end
    end
    if ~isempty(refSize)
        break;
    end
end
if isempty(refSize)
    error('consensus_reflex_structure:NoValidReplicate', ...
        'No valid, non-empty replicate mask found (opts.valid has no usable entries).');
end

for k = 1:K
    for g = 1:nGroups
        if ~valid(k, g)
            continue;   % failed/invalid replicate: mask may be empty, ignore it
        end
        if isempty(masks{k, g}) || ~isequal(size(masks{k, g}), refSize)
            error('consensus_reflex_structure:SizeMismatch', ...
                'masks{%d,%d} is %s but the reference mask is %s.', ...
                k, g, mat2str(size(masks{k, g})), mat2str(refSize));
        end
    end
end

C = struct();
C.rule           = opts.rule;
C.piThreshold    = opts.piThreshold;
C.tieBreak       = opts.tieBreak;
C.clusterJaccard = opts.clusterJaccard;
C.nSeeds         = K;
C.nGroups        = nGroups;
C.nElem          = prod(refSize);
C.valid          = valid;

for g = 1:nGroups
    vk = find(valid(:, g));
    n = numel(vk);
    if n == 0
        error('consensus_reflex_structure:NoValidReplicate', ...
            'Group %d has no valid replicate.', g);
    end

    % ---- stack the valid masks as an n x nElem logical matrix -----------
    M = false(n, C.nElem);
    for i = 1:n
        M(i, :) = masks{vk(i), g}(:)';
    end

    Pi = reshape(mean(M, 1), refSize);

    % ---- exact-mask frequencies ----------------------------------------
    keys = cell(n, 1);
    for i = 1:n
        keys{i} = char('0' + M(i, :));
    end
    [~, ~, ic] = unique(keys);
    counts = accumarray(ic, 1);
    maxCount = max(counts);

    % ---- Jaccard similarity + clusters ---------------------------------
    X = double(M);
    inter = X * X';
    sz = sum(X, 2);
    uni = sz + sz' - inter;
    J = inter ./ uni;
    J(uni == 0) = 1;                       % two empty masks are identical
    J(1:n + 1:end) = 1;

    adj = (J >= opts.clusterJaccard);
    adj = adj | adj';
    labels = local_components(adj);
    clusterSizes = accumarray(labels, 1);
    [domSize, domLabel] = max(clusterSizes);
    domFraction = domSize / n;

    % ---- modal mask (ties -> highest mean Jaccard) ---------------------
    cand = find(counts == maxCount);
    modeTie = numel(cand) > 1;
    bestScore = -inf;
    modeMask = reshape(M(find(ic == cand(1), 1), :), refSize);
    for c = cand(:)'
        rows = find(ic == c);
        score = mean(mean(J(rows, :), 2));
        if score > bestScore + 1e-12
            bestScore = score;
            modeMask = reshape(M(rows(1), :), refSize);
        end
    end

    % ---- medoid --------------------------------------------------------
    meanJ = mean(J, 2);
    [~, imed] = max(meanJ);
    medoidMask = reshape(M(imed, :), refSize);

    % ---- consensus mask ------------------------------------------------
    switch lower(opts.rule)
        case 'majority'
            if strcmpi(opts.tieBreak, 'exclude')
                consensusMask = Pi > 0.5;
            else
                consensusMask = Pi >= 0.5 - 1e-12;
            end
        case 'stability'
            consensusMask = Pi >= opts.piThreshold - 1e-12;
        case 'mode'
            consensusMask = modeMask;
        case 'medoid'
            consensusMask = medoidMask;
        otherwise
            error('consensus_reflex_structure:BadRule', ...
                'Unknown rule "%s" (majority|mode|medoid|stability).', opts.rule);
    end
    consensusMask = logical(consensusMask);

    % ---- how many replicates actually support the CHOSEN framework ------
    % Exact-mask agreement (modeFraction) is brittle: with ~230 binary entries
    % two runs practically never coincide bit-for-bit.  The meaningful
    % question is "how many seeds are close to the framework we selected",
    % measured by the Jaccard distance to the consensus mask itself.
    cm = consensusMask(:).';
    jacToCons = zeros(n, 1);
    covToCons = zeros(n, 1);
    for i = 1:n
        inter = nnz(M(i, :) & cm);
        uni = nnz(M(i, :) | cm);
        jacToCons(i) = inter / max(uni, 1);
        % coverage: share of THIS replicate's structure retained by the
        % consensus.  For a majority-vote superset, Jaccard to a single seed is
        % low by construction; coverage is the meaningful "does the consensus
        % still contain what this seed needed" number.
        covToCons(i) = inter / max(nnz(M(i, :)), 1);
    end
    nConforming = nnz(jacToCons >= opts.conformJaccard - 1e-12);

    % ---- diagnostics ---------------------------------------------------
    nOn  = nnz(Pi >= opts.stableOnPi - 1e-12);
    nOff = nnz(Pi <= opts.stableOffPi + 1e-12);

    C.mask{g}                   = consensusMask;
    C.Pi{g}                     = Pi;
    C.nValid(g)                 = n;
    C.modeMask{g}               = modeMask;
    C.modeCount(g)              = maxCount;
    C.modeFraction(g)           = maxCount / n;
    C.modeTie(g)                = modeTie;
    C.medoidIdx(g)              = vk(imed);
    C.medoidMask{g}             = medoidMask;
    C.jaccardMean(g)            = mean(J(~eye(n, 'logical')));
    C.jaccardMatrix{g}          = J;
    C.clusters{g}               = labels;
    C.clusterSizes{g}           = clusterSizes;
    C.dominantCluster(g)        = domLabel;
    C.dominantClusterFraction(g) = domFraction;
    C.jaccardToConsensus{g}      = jacToCons;
    C.coverage{g}                = covToCons;
    C.meanCoverage(g)            = mean(covToCons);
    C.nConforming(g)             = nConforming;
    C.conformingFraction(g)      = nConforming / n;
    C.nStableOn(g)              = nOn;
    C.nStableOff(g)             = nOff;
    C.nAmbiguous(g)             = C.nElem - nOn - nOff;
    C.supportSizes{g}           = sz;
    C.exactMaskCounts{g}        = sort(counts, 'descend');

    if n == 1
        C.jaccardMean(g) = NaN;   % no pair exists
    end
end
end

%% ========================================================================
% Local: fill option defaults.
% ========================================================================
function opts = fill_defaults(opts)
d = struct('rule', 'majority', 'piThreshold', 0.7, 'tieBreak', 'include', ...
    'stableOnPi', 0.9, 'stableOffPi', 0.1, 'clusterJaccard', 0.9, ...
    'conformJaccard', 0.8, 'valid', []);
fn = fieldnames(d);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}) || isempty(opts.(fn{i}))
        opts.(fn{i}) = d.(fn{i});
    end
end
if ~any(strcmpi(opts.rule, {'majority', 'mode', 'medoid', 'stability'}))
    error('consensus_reflex_structure:BadRule', 'Unknown rule "%s".', opts.rule);
end
if ~any(strcmpi(opts.tieBreak, {'include', 'exclude'}))
    error('consensus_reflex_structure:BadTieBreak', ...
        'tieBreak must be ''include'' or ''exclude''.');
end
if ~(isscalar(opts.piThreshold) && opts.piThreshold > 0 && opts.piThreshold <= 1)
    error('consensus_reflex_structure:BadThreshold', 'piThreshold must be in (0, 1].');
end
end

%% ========================================================================
% Local: connected components of a symmetric adjacency matrix (BFS).
% ========================================================================
function labels = local_components(adj)
n = size(adj, 1);
labels = zeros(n, 1);
cur = 0;
for i = 1:n
    if labels(i) ~= 0
        continue;
    end
    cur = cur + 1;
    stack = i;
    while ~isempty(stack)
        v = stack(end);
        stack(end) = [];
        if labels(v) ~= 0
            continue;
        end
        labels(v) = cur;
        nb = find(adj(v, :) & labels(:).' == 0);
        stack = [stack; nb(:)]; %#ok<AGROW>
    end
end
end
