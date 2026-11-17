function [idxOut, info] = bootstrap_phase_runs(idxIn, mode, opts)
%BOOTSTRAP_PHASE_RUNS  Resample the samples of one gait-phase group by runs.
%
%   [idxOut, info] = bootstrap_phase_runs(idxIn, mode, opts)
%
%   idxIn is the vector of (already delay-aligned) row indices belonging to one
%   phase group, e.g. D.phaseIdxMap.stance.  Because walking data are strongly
%   autocorrelated, resampling individual samples would understate the true
%   variability of the identified structure.  Instead, idxIn is split into its
%   natural contiguous runs (one stance / one swing run per stride) and whole
%   runs are resampled.  This is the "same experiment, different draw" that
%   each random seed of the stability identification stands for.
%
%   Modes
%   -----
%   'run_bootstrap'  (default) draw nBlocks runs WITH replacement and
%                    concatenate them.  The output length is therefore close to
%                    the input length (sum of nBlocks run lengths) but not
%                    fixed; no run is ever cut.
%   'run_subsample'  draw ceil(subsampleFraction * nBlocks) runs WITHOUT
%                    replacement.  Fewer samples than the original, so it
%                    conflates "different data" with "less data" -- useful as an
%                    ablation, not as the default.
%   'iid_bootstrap'  draw numel(idxIn) samples with replacement (classic
%                    bootstrap).  Kept for comparison only: it breaks the
%                    temporal correlation and is expected to look too stable.
%
%   opts
%   ----
%   randStream        RandStream object (required unless opts.seed is given).
%                     Passing an explicit stream keeps parfor workers
%                     reproducible and leaves the global RNG untouched.
%   seed              scalar seed; used to build a Threefry stream if
%                     opts.randStream is absent.
%   maxBlockLen       0 (default) keeps each run as one block; > 0 splits every
%                     run into consecutive blocks of at most that many samples
%                     (finer granularity for short records).
%   subsampleFraction fraction of blocks kept by 'run_subsample' (default 0.8).
%
%   info
%   ----
%   .mode, .nRuns, .nBlocks, .nOut, .nIn, .runLengths, .blockLengths,
%   .drawIdx (indices of the drawn blocks), .seed
%
%   See also lasso_stability_identification, consensus_reflex_structure.

if nargin < 3 || isempty(opts)
    opts = struct();
end
if nargin < 2 || isempty(mode)
    mode = 'run_bootstrap';
end
if ~isfield(opts, 'maxBlockLen') || isempty(opts.maxBlockLen)
    opts.maxBlockLen = 0;
end
if ~isfield(opts, 'subsampleFraction') || isempty(opts.subsampleFraction)
    opts.subsampleFraction = 0.8;
end
if ~isfield(opts, 'seed') || isempty(opts.seed)
    opts.seed = [];
end

idxIn = idxIn(:);
if isempty(idxIn)
    error('bootstrap_phase_runs:EmptyInput', 'idxIn must not be empty.');
end
if any(idxIn < 1) || any(abs(idxIn - round(idxIn)) > 0)
    error('bootstrap_phase_runs:BadIndex', 'idxIn must contain positive integers.');
end
if numel(unique(idxIn)) ~= numel(idxIn)
    idxIn = unique(idxIn);
    warning('bootstrap_phase_runs:DuplicateIndex', ...
        'Duplicate indices removed before resampling.');
end
idxIn = sort(idxIn);

% ---- resolve the random stream (explicit stream > seed > error) ----------
if isfield(opts, 'randStream') && ~isempty(opts.randStream)
    rs = opts.randStream;
    seedUsed = NaN;
elseif ~isempty(opts.seed)
    rs = RandStream('Threefry', 'Seed', opts.seed);
    seedUsed = opts.seed;
else
    error('bootstrap_phase_runs:NoRng', ...
        'Provide opts.randStream or opts.seed so the draw is reproducible.');
end

% ---- split into contiguous runs -----------------------------------------
d = diff(idxIn);
runStarts = [1; find(d > 1) + 1];
runEnds   = [find(d > 1); numel(idxIn)];
nRuns = numel(runStarts);
runLengths = runEnds - runStarts + 1;

% ---- optionally split runs into shorter blocks --------------------------
blocks = cell(0, 1);
if opts.maxBlockLen > 0
    for r = 1:nRuns
        s = idxIn(runStarts(r):runEnds(r));
        nChunk = ceil(numel(s) / opts.maxBlockLen);
        for c = 1:nChunk
            a = (c - 1) * opts.maxBlockLen + 1;
            b = min(c * opts.maxBlockLen, numel(s));
            blocks{end + 1, 1} = s(a:b); %#ok<AGROW>
        end
    end
else
    for r = 1:nRuns
        blocks{end + 1, 1} = idxIn(runStarts(r):runEnds(r)); %#ok<AGROW>
    end
end
nBlocks = numel(blocks);
blockLengths = cellfun(@numel, blocks);

% ---- draw ---------------------------------------------------------------
switch lower(mode)
    case 'run_bootstrap'
        % NOTE: randi(s, imax, m, n) -- the (m, n) arguments must be separate
        % scalars.  randi(rs, nBlocks, 1) would be parsed as imax=nBlocks, n=1
        % and silently return a 1x1 scalar (a single block!), so write out all
        % three dimensions explicitly.
        drawIdx = randi(rs, nBlocks, nBlocks, 1);
        idxOut = vertcat(blocks{drawIdx});
    case 'run_subsample'
        k = max(1, ceil(opts.subsampleFraction * nBlocks));
        drawIdx = randperm(rs, nBlocks, k).';
        drawIdx = sort(drawIdx);
        idxOut = vertcat(blocks{drawIdx});
    case 'iid_bootstrap'
        drawIdx = randi(rs, numel(idxIn), numel(idxIn), 1);
        idxOut = idxIn(drawIdx);
    otherwise
        error('bootstrap_phase_runs:BadMode', ...
            'Unknown resampling mode "%s" (expected run_bootstrap, run_subsample or iid_bootstrap).', ...
            mode);
end

idxOut = idxOut(:);

% Sanity guard: a run-level bootstrap must return roughly the original number
% of samples.  A much shorter draw means the block indices were built wrong.
nOut = numel(idxOut);
if strcmpi(mode, 'run_bootstrap') && nOut < 0.6 * numel(idxIn)
    error('bootstrap_phase_runs:DegenerateDraw', ...
        ['The bootstrap draw returned %d of %d samples (%.0f%%); expected ' ...
         'roughly the full sample count.  Check the block-index drawing.'], ...
        nOut, numel(idxIn), 100 * nOut / numel(idxIn));
end
if strcmpi(mode, 'run_subsample') && nOut < 0.3 * numel(idxIn)
    error('bootstrap_phase_runs:DegenerateDraw', ...
        'The subsample draw returned only %d of %d samples.', nOut, numel(idxIn));
end

info = struct();
info.mode         = mode;
info.nRuns        = nRuns;
info.nBlocks      = nBlocks;
info.nIn          = numel(idxIn);
info.nOut         = numel(idxOut);
info.runLengths   = runLengths;
info.blockLengths = blockLengths;
info.drawIdx      = drawIdx;
info.seed         = seedUsed;

if nBlocks == 1
    warning('bootstrap_phase_runs:SingleBlock', ...
        'The phase group contains a single contiguous run; the draw is degenerate.');
end
end
