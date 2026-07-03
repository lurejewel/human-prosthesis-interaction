function test_lasso_reflex_controller()
% Name: test_lasso_reflex_controller
% Description: Pure-MATLAB tests for LASSO reflex controller
%   loader / unpacker round-trip.  Tests both legacy per-phase and
%   grouped controller formats.  Does NOT require OpenSim.

fprintf('=== LASSO reflex controller tests ===\n');
addpath(genpath(fullfile('..', 'functions', 'core')));

tempFile = [tempname, '.mat'];

try
    %% ---- build synthetic LASSO struct (legacy per-phase format) ----
    rng(2026, 'threefry');
    nPhases = 5;
    nFeatures = 26;
    nMusclesPerSide = 9;

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
        lasso.mask{p} = rand(26, 9) > 0.7;
        betaMat = zeros(26, 9);
        betaMat(lasso.mask{p}) = randn(nnz(lasso.mask{p}), 1);
        lasso.beta{p} = betaMat;
        lasso.bias{p} = randn(1, 9);
        origBeta{p} = betaMat;
        origBias{p} = lasso.bias{p};
    end

    save(tempFile, 'lasso');

    %% ---- Test 1: legacy loader ----
    fprintf('Test 1: legacy loader ... ');
    [initPara, reflexParamMap, reflexTemplate] = ...
        load_lasso_reflex_controller(tempFile);

    expectedLen = 0;
    for p = 1:nPhases
        expectedLen = expectedLen + nnz(lasso.mask{p}) + 9;
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

    %% ---- Test 2: legacy round-trip unpack ----
    fprintf('Test 2: legacy round-trip unpack ... ');
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

    %% ---- Test 5: grouped format loader ----
    fprintf('Test 5: grouped format loader ... ');

    % Build a grouped LASSO struct: stance (phases 0,1) + swing (phases 2,3,4)
    nGroups = 2;
    groupLabels = {'stance', 'swing'};
    groupPhases = {[0 1], [2 3 4]};

    grpLasso = struct();
    grpLasso.format = 'grouped';
    grpLasso.nPhases = 5;
    grpLasso.phaseIds = 0:4;
    grpLasso.nGroups = nGroups;
    grpLasso.groups = struct();

    groupOrigBeta = cell(nGroups, 1);
    groupOrigBias = cell(nGroups, 1);
    groupOrigMask = cell(nGroups, 1);

    for g = 1:nGroups
        maskG = rand(26, 9) > 0.7;
        betaG = zeros(26, 9);
        betaG(maskG) = randn(nnz(maskG), 1);
        biasG = randn(1, 9);

        grpLasso.groups(g).label  = groupLabels{g};
        grpLasso.groups(g).phases = groupPhases{g};
        grpLasso.groups(g).beta   = betaG;
        grpLasso.groups(g).bias   = biasG;
        grpLasso.groups(g).mask   = maskG;

        groupOrigBeta{g} = betaG;
        groupOrigBias{g} = biasG;
        groupOrigMask{g} = maskG;
    end

    grpFile = [tempname, '.mat'];
    save(grpFile, 'grpLasso');
    cleanupGrp = onCleanup(@() delete(grpFile));

    [initParaGrp, reflexParamMapGrp, reflexTemplateGrp] = ...
        load_lasso_reflex_controller(grpFile);

    % verify no parameter duplication: totalLen = sum_g (nnz(mask_g) + 7)
    expectedGrpLen = 0;
    for g = 1:nGroups
        expectedGrpLen = expectedGrpLen + nnz(groupOrigMask{g}) + 9;
    end
    assert(numel(initParaGrp) == expectedGrpLen, ...
        'Grouped initPara length %d != expected %d.', numel(initParaGrp), expectedGrpLen);
    assert(reflexParamMapGrp.totalLen == expectedGrpLen);

    % verify map metadata
    assert(isfield(reflexParamMapGrp, 'nGroups'));
    assert(reflexParamMapGrp.nGroups == nGroups);
    assert(isfield(reflexTemplateGrp, 'groupOfPhase'));
    assert(isfield(reflexTemplateGrp, 'groupLabels'));
    assert(isfield(reflexTemplateGrp, 'groupPhases'));

    % verify groupOfPhase mapping
    assert(reflexTemplateGrp.groupOfPhase(1) == 1);  % phase 0 → stance
    assert(reflexTemplateGrp.groupOfPhase(2) == 1);  % phase 1 → stance
    assert(reflexTemplateGrp.groupOfPhase(3) == 2);  % phase 2 → swing
    assert(reflexTemplateGrp.groupOfPhase(4) == 2);  % phase 3 → swing
    assert(reflexTemplateGrp.groupOfPhase(5) == 2);  % phase 4 → swing

    % verify group params in map
    for g = 1:nGroups
        gm = reflexParamMapGrp.groups(g);
        assert(strcmp(gm.label, groupLabels{g}));
        assert(isequal(gm.phases, groupPhases{g}));
        assert(gm.nBeta == nnz(groupOrigMask{g}));
    end

    fprintf('PASSED\n');

    %% ---- Test 6: grouped format round-trip unpack ----
    fprintf('Test 6: grouped round-trip unpack ... ');
    reflexParamsGrp = unpack_lasso_reflex_params( ...
        initParaGrp, reflexParamMapGrp, reflexTemplateGrp);

    % phases 0,1 must share the stance beta/bias
    for pIdx = [1 2]
        betaP = reflexParamsGrp.beta{pIdx};
        biasP = reflexParamsGrp.bias{pIdx};
        assert(max(abs(betaP(:) - groupOrigBeta{1}(:))) < 1e-12, ...
            'Phase %d: beta does not match stance group.', pIdx-1);
        assert(max(abs(biasP(:) - groupOrigBias{1}(:))) < 1e-12, ...
            'Phase %d: bias does not match stance group.', pIdx-1);
        assert(all(betaP(~groupOrigMask{1}) == 0), ...
            'Phase %d: non-zero beta outside stance mask.', pIdx-1);
    end

    % phases 2,3,4 must share the swing beta/bias
    for pIdx = [3 4 5]
        betaP = reflexParamsGrp.beta{pIdx};
        biasP = reflexParamsGrp.bias{pIdx};
        assert(max(abs(betaP(:) - groupOrigBeta{2}(:))) < 1e-12, ...
            'Phase %d: beta does not match swing group.', pIdx-1);
        assert(max(abs(biasP(:) - groupOrigBias{2}(:))) < 1e-12, ...
            'Phase %d: bias does not match swing group.', pIdx-1);
        assert(all(betaP(~groupOrigMask{2}) == 0), ...
            'Phase %d: non-zero beta outside swing mask.', pIdx-1);
    end

    % verify reflexParams has per-phase masks (not group masks)
    assert(numel(reflexParamsGrp.mask) == 5);
    assert(isequal(reflexParamsGrp.mask{1}, groupOrigMask{1}));
    assert(isequal(reflexParamsGrp.mask{2}, groupOrigMask{1}));
    assert(isequal(reflexParamsGrp.mask{3}, groupOrigMask{2}));
    assert(isequal(reflexParamsGrp.mask{4}, groupOrigMask{2}));
    assert(isequal(reflexParamsGrp.mask{5}, groupOrigMask{2}));

    fprintf('PASSED\n');

    %% ---- Test 7: grouped format — perturb & verify no crosstalk ----
    fprintf('Test 7: grouped format — parameter independence ... ');

    % Clone the grouped initPara and perturb stance bias only
    initParaMod = initParaGrp;
    gmStance = reflexParamMapGrp.groups(1);
    % add 0.1 to stance bias entries
    initParaMod(gmStance.biasStart : gmStance.biasEnd) = ...
        initParaMod(gmStance.biasStart : gmStance.biasEnd) + 0.1;

    reflexParamsMod = unpack_lasso_reflex_params( ...
        initParaMod, reflexParamMapGrp, reflexTemplateGrp);

    % stance phases should reflect the change
    for pIdx = [1 2]
        biasDiff = reflexParamsMod.bias{pIdx} - reflexParamsGrp.bias{pIdx};
        assert(max(abs(biasDiff - 0.1)) < 1e-12, ...
            'Phase %d: stance bias not perturbed correctly.', pIdx-1);
        % beta should be unchanged
        betaDiff = reflexParamsMod.beta{pIdx} - reflexParamsGrp.beta{pIdx};
        assert(max(abs(betaDiff(:))) < 1e-12, ...
            'Phase %d: stance beta changed unexpectedly.', pIdx-1);
    end

    % swing phases should be completely unchanged
    for pIdx = [3 4 5]
        biasDiff = reflexParamsMod.bias{pIdx} - reflexParamsGrp.bias{pIdx};
        assert(max(abs(biasDiff)) < 1e-12, ...
            'Phase %d: swing bias changed unexpectedly.', pIdx-1);
        betaDiff = reflexParamsMod.beta{pIdx} - reflexParamsGrp.beta{pIdx};
        assert(max(abs(betaDiff(:))) < 1e-12, ...
            'Phase %d: swing beta changed unexpectedly.', pIdx-1);
    end

    fprintf('PASSED\n');

    fprintf('\n=== ALL TESTS PASSED ===\n');

catch ME
    fprintf('\n*** TEST FAILED ***\n%s\n', getReport(ME));
    rethrow(ME);
end

% cleanup temp file
delete(tempFile);

end
