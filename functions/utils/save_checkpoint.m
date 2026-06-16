function save_checkpoint(checkpointFile, g, gAtRestartStart, lastImprovementGen, ...
                        bestFitEver, softBoostTriggered, restartCount, stop, ...
                        bestFits, computationLoad, optConfig, optCfg, projName, ...
                        reflexParamMap, reflexTemplate, initPara, sigma, nWorkers, ...
                        simConfig, initPose, lassoFile)
% Name: save_checkpoint
% Description: Save all CMA-ES optimization state for crash recovery.
%   Call this periodically (e.g. every 10 generations) and at the end of
%   the optimization. Use the companion load_checkpoint() to resume.
%
% Saved state includes:
%   - CMA-ES internal state (population, step-size, covariance, paths, etc.)
%   - Loop bookkeeping (generation counter, stall tracking, restart counter)
%   - Immutable configuration (controller metadata, problem setup)
%
% Usage:
%   save_checkpoint('results/checkpoint.mat', g, gAtRestartStart, ...
%       lastImprovementGen, bestFitEver, softBoostTriggered, restartCount, ...
%       stop, bestFits, computationLoad, optConfig, optCfg, projName, ...
%       reflexParamMap, reflexTemplate, initPara, sigma, nWorkers, lassoFile);

% ---- trim history arrays to actual length ----
bestFits        = bestFits(1:g);
computationLoad = computationLoad(1:g);

% ---- assemble checkpoint struct ----
ckp = struct();

% loop bookkeeping
ckp.g                  = g;
ckp.gAtRestartStart    = gAtRestartStart;
ckp.lastImprovementGen = lastImprovementGen;
ckp.bestFitEver        = bestFitEver;
ckp.softBoostTriggered = softBoostTriggered;
ckp.restartCount       = restartCount;
ckp.stop               = stop;

% history
ckp.bestFits         = bestFits;
ckp.computationLoad  = computationLoad;

% CMA-ES object (full handle object — includes core, records, bookkeeping)
ckp.optConfig = optConfig;

% immutable configuration
ckp.optCfg          = optCfg;
ckp.projName        = projName;
ckp.simConfig       = simConfig;        % needed to recompute modelStaticProp on resume
ckp.initPose        = initPose;
ckp.reflexParamMap  = reflexParamMap;   % LASSO controller metadata
ckp.reflexTemplate  = reflexTemplate;
ckp.initPara        = initPara;
ckp.sigma           = sigma;
ckp.nWorkers        = nWorkers;

% optional LASSO file path (may be empty string if warm-start)
if exist('lassoFile', 'var') && ~isempty(lassoFile)
    ckp.lassoFile = lassoFile;
else
    ckp.lassoFile = '';
end

ckp.timestamp = datetime;

% ---- atomic save (write to temp first, then rename) ----
tmpFile = [checkpointFile '.tmp'];
save(tmpFile, '-struct', 'ckp', '-v7.3');
movefile(tmpFile, checkpointFile, 'f');

fprintf('[%s] Checkpoint saved to %s (gen %d, bestFit=%.6g).\n', ...
    char(datetime), checkpointFile, g, optConfig.recordForBestParticle.fit);

end
