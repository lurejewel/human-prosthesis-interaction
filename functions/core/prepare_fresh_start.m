function s = prepare_fresh_start(projName)
% Name: prepare_fresh_start
% Description: Set up a fresh CMA-ES optimisation run for the muscle-reflex
%   predictive simulation.  Reads initial kinematics and joint moments from
%   UN.sto, runs iterative static optimisation to obtain initial muscle
%   activations (a_opt), loads the LASSO controller, and initialises the
%   CMA-ES optimiser.  Returns a single struct s with all fields needed by
%   the main optimisation loop.
%
% Input:
%   projName  - model base name (e.g. 'human0714')
%
% Output struct fields:
%   optCfg, simConfig, initPose, a_opt, modelStaticProp, ...
%   reflexParamMap, reflexTemplate, initPara, sigma, ...
%   optConfig, nWorkers, histLen, bestFits, computationLoad, ...
%   g, gAtRestartStart, lastImprovementGen, bestFitEver, ...
%   softBoostTriggered, restartCount, stop, lassoFile

%% ---- reproducibility configuration ----
optCfg.rngSeed         = 2026;
optCfg.softPatience    = 30;
optCfg.patience        = 100;
optCfg.maxRestarts     = 3;
optCfg.sigmaBoostFactor = 2.0;
optCfg.minImprovement  = 1e-12;
rng(optCfg.rngSeed, 'threefry');

%% ---- activate parallel computing ----
nWorkers = feature('numcores');
p = gcp('nocreate');
if isempty(p) || p.NumWorkers ~= nWorkers
    delete(p);
    parpool('local', nWorkers);
end

%% ---- simulation configuration ----
simConfig.endTime  = 10;
simConfig.stepTime = 0.005;
simConfig.saveSTO  = 1;
simConfig.speed    = 1.0;
simConfig.slope    = 0;

%% ---- read initial kinematics and joint moments from UN.sto ----
stoPath = 'Sparse Group LASSO Validation\UN.sto';
assert(isfile(stoPath), 'UN.sto not found: %s', stoPath);

fid = fopen(stoPath, 'r');
for k = 1:6, fgetl(fid); end
headerLine = fgetl(fid);
dataLine   = str2double(strsplit(fgetl(fid), '\t'));
fclose(fid);
colNames = strsplit(headerLine, '\t');
colMap = containers.Map('KeyType', 'char', 'ValueType', 'int32');
for i = 1:numel(colNames), colMap(colNames{i}) = i; end

% 9 DOF (Q then U → 18 elements)
dofNames = {'pelvis_tilt','pelvis_tx','pelvis_ty', ...
            'hip_flexion_r','knee_flexion_r','ankle_dorsiflexion_r', ...
            'hip_flexion_l','knee_flexion_l','ankle_dorsiflexion_l'};
nDof = numel(dofNames);
initPose = zeros(2 * nDof, 1);
for j = 1:nDof
    initPose(j)        = dataLine(colMap(dofNames{j}));
    initPose(nDof + j) = dataLine(colMap([dofNames{j} '_u']));
end

% Joint moments + knee limit torques for static optimisation
soCoordNames = {'hip_flexion_r','hip_flexion_l', ...
                'knee_flexion_r','knee_flexion_l', ...
                'ankle_dorsiflexion_r','ankle_dorsiflexion_l'};
nSoCoords = numel(soCoordNames);
tauTarget = zeros(nSoCoords, 1);
for j = 1:nSoCoords
    tauTarget(j) = dataLine(colMap([soCoordNames{j} '.moment']));
end
tauLimit = zeros(nSoCoords, 1);
tauLimit(strcmp(soCoordNames, 'knee_flexion_r')) = dataLine(colMap('knee_r.torque'));
tauLimit(strcmp(soCoordNames, 'knee_flexion_l')) = dataLine(colMap('knee_l.torque'));

%% ---- model static properties ----
modelStaticProp = read_muscle_static_prop(projName, simConfig, initPose);

%% ---- iterative static optimisation for initial muscle activations ----
fprintf('[%s] Running iterative static optimisation ...\n', char(datetime));
a_opt = iterative_static_optimization(projName, initPose, dofNames, ...
    tauTarget, tauLimit, soCoordNames);

%% ---- LASSO linear-phase controller initialisation ----
lassoFile = fullfile('results', 'lasso_controller_result.mat');
[initPara, reflexParamMap, reflexTemplate] = ...
    load_lasso_reflex_controller(lassoFile);

sigma = 0.02;
optConfig = CMAES_optimization(initPara, sigma, nWorkers);

%% ---- history arrays & loop bookkeeping ----
histLen = (optCfg.maxRestarts + 1) * (optCfg.patience + 50) + 200;
bestFits        = nan(histLen, 1);
computationLoad = nan(histLen, 1);

g                  = 0;
gAtRestartStart    = 0;
lastImprovementGen = 0;
bestFitEver        = inf;
softBoostTriggered = false;
restartCount       = 0;
stop               = false;

%% ---- assemble output struct ----
s.optCfg              = optCfg;
s.simConfig           = simConfig;
s.initPose            = initPose;
s.a_opt               = a_opt;
s.modelStaticProp     = modelStaticProp;
s.reflexParamMap      = reflexParamMap;
s.reflexTemplate      = reflexTemplate;
s.initPara            = initPara;
s.sigma               = sigma;
s.optConfig           = optConfig;
s.nWorkers            = nWorkers;
s.histLen             = histLen;
s.bestFits            = bestFits;
s.computationLoad     = computationLoad;
s.g                   = g;
s.gAtRestartStart     = gAtRestartStart;
s.lastImprovementGen  = lastImprovementGen;
s.bestFitEver         = bestFitEver;
s.softBoostTriggered  = softBoostTriggered;
s.restartCount        = restartCount;
s.stop                = stop;
s.lassoFile           = lassoFile;

end
