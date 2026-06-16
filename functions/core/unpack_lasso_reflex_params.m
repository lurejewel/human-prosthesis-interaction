function reflexParams = unpack_lasso_reflex_params( ...
    paramVec, reflexParamMap, reflexTemplate)
% Name: unpack_lasso_reflex_params
% Description: Reconstruct structured controller parameters (beta cells and
%   bias cells) from a flat CMA-ES parameter vector, using the layout
%   described by reflexParamMap and the sparse mask from reflexTemplate.
%
%   Supports both 'grouped' and 'legacy' reflexParamMap layouts.
%   In grouped mode, the flat vector contains one set of parameters per
%   controller group (no duplication); these are expanded to per-phase
%   beta/bias cells for runtime use.
%
% Input:
%   paramVec       – N×1 flat parameter vector (or 1×N, automatically reshaped)
%   reflexParamMap – struct from load_lasso_reflex_controller
%   reflexTemplate – struct from load_lasso_reflex_controller
% Output:
%   reflexParams   – struct with fields: controllerType, nFeatures,
%                    nMusclesPerSide, nPhases, phaseIds, beta, bias, mask

paramVec = paramVec(:);
assert(numel(paramVec) == reflexParamMap.totalLen, ...
    'Parameter vector length (%d) does not match expected total (%d).', ...
    numel(paramVec), reflexParamMap.totalLen);

nPhases = reflexParamMap.nPhases;

reflexParams = struct();
reflexParams.controllerType  = 'lasso_linear_phase';
reflexParams.nFeatures       = reflexParamMap.nFeatures;
reflexParams.nMusclesPerSide = reflexParamMap.nMusclesPerSide;
reflexParams.nPhases         = nPhases;
reflexParams.phaseIds        = reflexParamMap.phaseIds;
reflexParams.beta            = cell(nPhases, 1);
reflexParams.bias            = cell(nPhases, 1);
reflexParams.mask            = reflexTemplate.mask;

if isfield(reflexParamMap, 'nGroups') && ~isempty(reflexParamMap.groups)
    % ======== GROUPED format ========
    % Unpack one beta/bias per group, then expand to phases.
    nGroups = reflexParamMap.nGroups;
    groupBeta = cell(nGroups, 1);
    groupBias = cell(nGroups, 1);

    for g = 1:nGroups
        gm = reflexParamMap.groups(g);

        % ---- reconstruct group beta ----
        betaMat = zeros(22, 7);
        betaMat(gm.betaLinIdx) = paramVec(gm.betaStart : gm.betaEnd);
        groupBeta{g} = betaMat;

        % ---- reconstruct group bias ----
        groupBias{g} = paramVec(gm.biasStart : gm.biasEnd)';  % 1×7
    end

    % ---- expand groups to per-phase cells ----
    for p = 1:nPhases
        g = reflexTemplate.groupOfPhase(p);
        reflexParams.beta{p} = groupBeta{g};
        reflexParams.bias{p} = groupBias{g};

        % ---- structural-zero enforcement ----
        maskMat = reflexTemplate.mask{p};
        betaMat = reflexParams.beta{p};
        betaMat(~maskMat) = 0;
        reflexParams.beta{p} = betaMat;
    end
else
    % ======== LEGACY per-phase format ========
    for p = 1:nPhases
        ph = reflexParamMap.phases(p);

        % ---- reconstruct beta ----
        betaMat = zeros(22, 7);
        betaMat(ph.betaLinIdx) = paramVec(ph.betaStart : ph.betaEnd);
        reflexParams.beta{p} = betaMat;

        % ---- reconstruct bias ----
        reflexParams.bias{p} = paramVec(ph.biasStart : ph.biasEnd)';  % 1×7

        % ---- structural-zero enforcement ----
        maskMat = reflexTemplate.mask{p};
        betaMat(~maskMat) = 0;
        reflexParams.beta{p} = betaMat;
    end
end

end
