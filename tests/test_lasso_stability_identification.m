function test_lasso_stability_identification()
% Name: test_lasso_stability_identification
% Description: Pure-MATLAB tests for the multi-seed consensus structure
%   identification (run-level bootstrap + Pi/majority consensus + frozen-mask
%   refit).  No OpenSim required: the design system is synthetic.
%
%   Covered
%   -------
%   1. majority vote returns the framework when 13/20 seeds agree;
%   2. Pi / modal-frequency / cluster statistics are consistent;
%   3. Pi == 0.5 tie-break (even K) is deterministic and configurable;
%   4. bootstrap is reproducible per seed and resamples whole runs;
%   5. end-to-end synthetic identification: deterministic per baseSeed,
%      structure invariants (beta == 0 outside mask), loader round-trip,
%      and Pi separates true from null entries.

fprintf('=== LASSO stability identification tests ===\n');

thisDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(thisDir);
lassoDir = fullfile(repoRoot, 'Sparse Group LASSO Validation');
addpath(lassoDir);
addpath(fullfile(repoRoot, 'functions', 'core'));

try
    %% ---- Test 1: 13/20 agreement -> that exact framework ----------------
    fprintf('Test 1: majority vote with 13/20 agreement ... ');
    K = 20;
    nFeat = 26;
    nMus = 9;
    rng(11, 'threefry');
    baseMask = (rand(nFeat, nMus) > 0.75);
    masks = cell(K, 1);
    for k = 1:13
        masks{k, 1} = baseMask;
    end
    for k = 14:K
        masks{k, 1} = (rand(nFeat, nMus) > 0.5);
    end
    C = consensus_reflex_structure(masks, struct('rule', 'majority'));
    assert(isequal(C.mask{1}, baseMask), ...
        'Majority vote did not return the 13/20 framework.');
    assert(C.modeCount(1) == 13, 'Modal count %d ~= 13.', C.modeCount(1));
    assert(abs(C.modeFraction(1) - 13/20) < 1e-12, 'Modal fraction wrong.');
    assert(C.nValid(1) == K, 'nValid wrong.');
    assert(C.dominantClusterFraction(1) >= 13/K - 1e-12, ...
        'Dominant cluster should contain the 13 identical masks.');
    assert(C.nStableOff(1) == nnz(C.Pi{1} <= 0.1 + 1e-12), 'stable-off count wrong.');
    fprintf('PASSED\n');

    %% ---- Test 2: Pi separates the agreed entries ------------------------
    fprintf('Test 2: Pi of agreed entries ... ');
    piBase = C.Pi{1}(baseMask);
    assert(min(piBase) >= 13/K - 1e-12, 'Pi of agreed entries < 13/20.');
    assert(max(C.Pi{1}(:)) <= 1, 'Pi above 1.');
    fprintf('PASSED\n');

    %% ---- Test 3: Pi == 0.5 tie-break ------------------------------------
    fprintf('Test 3: Pi == 0.5 tie-break ... ');
    m2 = cell(2, 1);
    m2{1, 1} = false(4, 2);
    m2{1, 1}(1, 1) = true;
    m2{2, 1} = false(4, 2);
    Cinc = consensus_reflex_structure(m2, struct('tieBreak', 'include'));
    Cexc = consensus_reflex_structure(m2, struct('tieBreak', 'exclude'));
    assert(Cinc.mask{1}(1, 1) == true, 'include tie-break failed.');
    assert(Cexc.mask{1}(1, 1) == false, 'exclude tie-break failed.');
    fprintf('PASSED\n');

    %% ---- Test 4: bootstrap reproducibility and run structure -----------
    fprintf('Test 4: run-level bootstrap ... ');
    % 10 runs of 40 separated by 10-sample gaps (with only 3 runs two seeds
    % collide with probability ~4%, which made this test flaky).
    nRun = 10;
    runLen = 40;
    gap = 10;
    idxRuns = zeros(nRun * runLen, 1);
    for r = 1:nRun
        idxRuns((r - 1) * runLen + (1:runLen)) = ...
            (r - 1) * (runLen + gap) + (1:runLen);
    end
    o1 = struct('maxBlockLen', 0, 'randStream', RandStream('Threefry', 'Seed', 7));
    o2 = struct('maxBlockLen', 0, 'randStream', RandStream('Threefry', 'Seed', 7));
    o3 = struct('maxBlockLen', 0, 'randStream', RandStream('Threefry', 'Seed', 8));
    [i1, info1] = bootstrap_phase_runs(idxRuns, 'run_bootstrap', o1);
    i2 = bootstrap_phase_runs(idxRuns, 'run_bootstrap', o2);
    i3 = bootstrap_phase_runs(idxRuns, 'run_bootstrap', o3);
    assert(isequal(i1, i2), 'Same seed must give the same draw.');
    assert(~isequal(i1, i3), 'Different seeds should give different draws.');
    assert(info1.nRuns == nRun, 'Expected %d contiguous runs, got %d.', ...
        nRun, info1.nRuns);
    % Regression: randi(rs, nBlocks, 1) used to return a single block, which
    % silently shrank every replicate to one gait run.
    assert(info1.nOut == info1.nIn, ...
        'Run bootstrap must return ~nIn samples, got %d of %d.', ...
        info1.nOut, info1.nIn);
    assert(numel(i1) == sum(info1.blockLengths(info1.drawIdx)), ...
        'Output length must equal the sum of the drawn block lengths.');
    assert(all(ismember(i1, idxRuns)), 'Draw contains foreign indices.');

    % every drawn block must be a contiguous run of the input
    % (rebuild the blocks and verify each is contiguous in idxRuns)
    pos = cumsum([1, info1.blockLengths(:).']);
    okContig = true;
    for b = unique(info1.drawIdx).'
        s = idxRuns(pos(b):pos(b+1) - 1);
        if numel(s) > 1 && any(diff(s) ~= 1)
            okContig = false;
        end
    end
    assert(okContig, 'A drawn block is not contiguous.');
    fprintf('PASSED\n');

    %% ---- Test 5: end-to-end synthetic identification --------------------
    fprintf('Test 5: synthetic end-to-end identification ... ');
    D = make_synthetic_design(2026);

    fitOpts = struct();
    fitOpts.gamma = 1;
    fitOpts.c = 0.01;
    fitOpts.rho = 50;
    fitOpts.delaySteps = 1;
    fitOpts.robustInitIter = 1;
    fitOpts.robustIter = 1;
    fitOpts.verbose = false;
    fitOpts.storeLambdaPath = false;
    fitOpts.groupLabels = {'stance', 'swing'};
    fitOpts.groupPhases = {[0 1], [2 3 4]};
    fitOpts.source = 'synthetic';
    fitOpts.lambdaPathConfig = struct('nLambda', 6, 'lambdaMinRatio', 1e-3, ...
        'tolSelect', 0.05, 'selectionMetric', 'post_refit_train_censored_rmse');

    opt = struct();
    opt.nSeeds = 4;
    opt.baseSeed = 2026;
    opt.resampleMode = 'run_bootstrap';
    opt.consensusRule = 'majority';
    opt.useParallel = false;
    opt.saveResults = false;
    opt.verbose = false;
    opt.nRetry = 2;
    opt.fitOptions = fitOpts;

    % The synthetic design leaves a few samples unlabeled on purpose (gaps
    % between runs), so the stitched-prediction NaN warning is expected here.
    prevWarn = warning('off', 'fit_adaptive_lasso_phase_model:PredictionGap');
    R1 = lasso_stability_from_design(D, opt);
    R2 = lasso_stability_from_design(D, opt);
    warning(prevWarn);

    assert(isequal(R1.consensus.mask, R2.consensus.mask), ...
        'Same baseSeed must give the same consensus mask.');
    assert(isequal(R1.lasso.groups(1).mask, R2.lasso.groups(1).mask), ...
        'Same baseSeed must give the same controller mask.');

    % --- structural invariants ---
    for g = 1:2
        beta = R1.lasso.groups(g).beta;
        mask = R1.lasso.groups(g).mask;
        assert(all(beta(~mask) == 0), ...
            'Group %d: beta is non-zero outside the frozen mask.', g);
        assert(nnz(mask) > 0, 'Group %d: empty consensus mask.', g);
    end
    assert(R1.totalLen == nnz(R1.lasso.groups(1).mask) + 9 + ...
        nnz(R1.lasso.groups(2).mask) + 9, 'totalLen mismatch.');

    % --- Pi separates true from null entries ---
    trueMask = D.trueMask;
    piAll = [R1.consensus.Pi{1}(:); R1.consensus.Pi{2}(:)];
    tmAll = [trueMask(:); trueMask(:)];
    piTrue = mean(piAll(tmAll));
    piNull = mean(piAll(~tmAll));
    assert(piTrue > piNull, ...
        'Pi should be higher on the true support (%.3f) than off it (%.3f).', ...
        piTrue, piNull);
    fprintf('(pi true %.2f vs null %.2f) ', piTrue, piNull);

    % --- loader / unpacker round-trip ---
    tmpFile = [tempname, '.mat'];
    cleanupTmp = onCleanup(@() delete_if_exists(tmpFile));
    lasso = R1.lasso;
    save(tmpFile, 'lasso');
    [initPara, paramMap, template] = load_lasso_reflex_controller(tmpFile);
    assert(paramMap.totalLen == R1.totalLen, ...
        'Loader totalLen %d ~= %d.', paramMap.totalLen, R1.totalLen);
    params = unpack_lasso_reflex_params(initPara, paramMap, template);
    for g = 1:2
        ph = template.groupPhases{g};
        phaseIdx = find(template.phaseIds == ph(1), 1);
        d = max(abs(params.beta{phaseIdx}(:) - R1.lasso.groups(g).beta(:)));
        assert(d < 1e-12, 'Group %d beta round-trip mismatch (%.3g).', g, d);
        assert(all(params.beta{phaseIdx}(~template.mask{phaseIdx}) == 0), ...
            'Group %d: structural zeros violated after unpack.', g);
    end
    fprintf('PASSED\n');

    fprintf('\n=== ALL STABILITY TESTS PASSED ===\n');
catch ME
    fprintf('\n*** TEST FAILED ***\n%s\n', getReport(ME));
    rethrow(ME);
end
end

%% ========================================================================
% Synthetic design system: 26 features, 9 muscles, plant-like sparse truth,
% censored at c, with stance/swing runs so the run bootstrap has structure.
% ========================================================================
function D = make_synthetic_design(seed)
rng(seed, 'threefry');
p = 26;
nMus = 9;
nStance = 500;
nSwing = 500;
n = nStance + nSwing;

% features: 18 muscle-ish (lcEN-like positive, fATN-like 0..1) + 8 kinematics
lcEN = 0.6 + 0.4 * rand(n, nMus);
fATN = 0.8 * rand(n, nMus);
kin = 0.5 * randn(n, 8);
X = [lcEN, fATN, kin];

% sparse truth: 4 entries per muscle per group, well separated in magnitude
Btrue = zeros(p, nMus, 2);
for g = 1:2
    for m = 1:nMus
        rows = randperm(p, 4).';
        Btrue(rows, m, g) = 0.6 + 0.8 * rand(4, 1);
    end
end

bias = 0.05 * ones(1, nMus);
noise = 0.01 * randn(n, nMus);
Y = zeros(n, nMus);
stanceIdx = (1:nStance).';
swingIdx = ((nStance + 1):n).';
Y(stanceIdx, :) = X(stanceIdx, :) * Btrue(:, :, 1) + bias + noise(stanceIdx, :);
Y(swingIdx, :) = X(swingIdx, :) * Btrue(:, :, 2) + bias + noise(swingIdx, :);
c = 0.01;
Y = max(Y, c);

% contiguous runs so that run-level bootstrap does something sensible
% 10 stance runs of 40 samples separated by 10-sample gaps, inside 1..500;
% swing runs are the same pattern shifted by 500 so the two groups are disjoint.
runLen = 40;
nRuns = 10;
gap = 10;
baseRuns = arrayfun(@(r) ((r - 1) * (runLen + gap) + (1:runLen)).', ...
    1:nRuns, 'UniformOutput', false);
stanceUse = vertcat(baseRuns{:});
swingUse = stanceUse + nStance;

muX = mean(X, 1);
sX = std(X, 1, 1);
sX(sX == 0) = 1;
X_std = (X - muX) ./ sX;

D = struct();
D.X = X;
D.X_std = X_std;
D.muX = muX;
D.sX = sX;
D.Y = Y;
D.w = ones(n, 1);
D.phaseIdxMap = struct('stance', stanceUse, 'swing', swingUse);
D.phaseMask = struct('stance', ismember((1:n).', stanceUse), ...
    'swing', ismember((1:n).', swingUse));
D.muscleNames = arrayfun(@(m) sprintf('mus%d_r', m), 1:nMus, 'UniformOutput', false);
D.featureNames = [arrayfun(@(i) sprintf('f%d_lCEN', i), 1:nMus, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('m%d_fATN', i), 1:nMus, 'UniformOutput', false), ...
    {'pelvis_tilt', 'pelvis_tilt_u', 'hip_flexion_r', 'hip_flexion_r_u', ...
     'knee_extension_r', 'knee_extension_r_u', 'ankle_dorsiflexion_r', ...
     'ankle_dorsiflexion_r_u'}].';
D.nSamples = n;
D.nRawSamples = n;
D.p = p;
D.nMus = nMus;
D.delaySteps = 1;
D.stoFile = 'synthetic';
D.trueMask = any(abs(Btrue) > 0, 3);   % 26 x 9 union over both groups
end

function delete_if_exists(f)
if isfile(f)
    delete(f);
end
end
