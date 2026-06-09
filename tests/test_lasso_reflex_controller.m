function test_lasso_reflex_controller()
% Name: test_lasso_reflex_controller
% Description: Pure-MATLAB tests for LASSO reflex controller
%   loader / unpacker round-trip.  Does NOT require OpenSim.

fprintf('=== LASSO reflex controller tests ===\n');
addpath(genpath(fullfile('..', 'functions', 'core')));

tempFile = [tempname, '.mat'];

try
    %% ---- build synthetic LASSO struct ----
    rng(2026, 'threefry');
    nPhases = 5;
    nFeatures = 22;
    nMusclesPerSide = 7;

    lasso = struct();
    lasso.nPhases = nPhases;
    lasso.phaseIds = 0:(nPhases-1);
    lasso.beta = cell(nPhases, 1);
    lasso.bias = cell(nPhases, 1);
    lasso.mask = cell(nPhases, 1);

    origBeta = cell(nPhases, 1);
    origBias = cell(nPhases, 1);
    for p = 1:nPhases
        % random sparse mask (~30% non-zero)
        lasso.mask{p} = rand(22, 7) > 0.7;
        betaMat = zeros(22, 7);
        betaMat(lasso.mask{p}) = randn(nnz(lasso.mask{p}), 1);
        lasso.beta{p} = betaMat;
        lasso.bias{p} = randn(1, 7);
        origBeta{p} = betaMat;
        origBias{p} = lasso.bias{p};
    end

    save(tempFile, 'lasso');

    %% ---- Test 1: loader ----
    fprintf('Test 1: loader ... ');
    [initPara, reflexParamMap, reflexTemplate] = ...
        load_lasso_reflex_controller(tempFile);

    expectedLen = 0;
    for p = 1:nPhases
        expectedLen = expectedLen + nnz(lasso.mask{p}) + 7;
    end
    assert(numel(initPara) == expectedLen, ...
        'initPara length %d != expected %d', numel(initPara), expectedLen);
    assert(reflexParamMap.totalLen == expectedLen);

    % verify map contents
    assert(reflexParamMap.nPhases == nPhases);
    assert(reflexTemplate.nPhases == nPhases);
    assert(isequal(reflexTemplate.phaseIds, 0:4));
    for p = 1:nPhases
        ph = reflexParamMap.phases(p);
        assert(ph.nBeta == nnz(lasso.mask{p}));
    end
    fprintf('PASSED\n');

    %% ---- Test 2: round-trip unpack ----
    fprintf('Test 2: round-trip unpack ... ');
    reflexParams = unpack_lasso_reflex_params( ...
        initPara, reflexParamMap, reflexTemplate);

    for p = 1:nPhases
        betaReconstructed = reflexParams.beta{p};
        biasReconstructed = reflexParams.bias{p};
        assert(max(abs(betaReconstructed(:) - origBeta{p}(:))) < 1e-12, ...
            'Phase %d: beta mismatch.', p);
        assert(max(abs(biasReconstructed(:) - origBias{p}(:))) < 1e-12, ...
            'Phase %d: bias mismatch.', p);
        % structural zeros
        assert(all(betaReconstructed(~lasso.mask{p}) == 0), ...
            'Phase %d: non-zero beta outside mask.', p);
    end
    fprintf('PASSED\n');

    %% ---- Test 3: inconsistent mask error ----
    fprintf('Test 3: inconsistent mask error ... ');
    badLasso = lasso;
    badLasso.mask{1}(1, 1) = false;
    badLasso.beta{1}(1, 1) = 5.0;  % non-zero where mask is false
    badFile = [tempname, '.mat'];
    save(badFile, 'badLasso');
    cleanupBad = onCleanup(@() delete(badFile));

    try
        load_lasso_reflex_controller(badFile);
        error('Expected error was not thrown.');
    catch ME
        assert(contains(ME.message, 'non-zero'), ...
            'Unexpected error: %s', ME.message);
    end
    fprintf('PASSED\n');

    %% ---- Test 4: single-struct variable (not named "lasso") ----
    fprintf('Test 4: single-struct variable ... ');
    myCtrl = lasso;  %#ok<NASGU>
    singleFile = [tempname, '.mat'];
    save(singleFile, 'myCtrl');
    cleanupSingle = onCleanup(@() delete(singleFile));
    [initPara2, ~, ~] = load_lasso_reflex_controller(singleFile);
    assert(numel(initPara2) == expectedLen);
    fprintf('PASSED\n');

    fprintf('\n=== ALL TESTS PASSED ===\n');

catch ME
    fprintf('\n*** TEST FAILED ***\n%s\n', getReport(ME));
    rethrow(ME);
end

% cleanup temp file
delete(tempFile);

end
