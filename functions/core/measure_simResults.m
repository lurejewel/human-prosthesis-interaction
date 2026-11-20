function [fit, parts] = measure_simResults(modelInfo)
% Name: measure_simResults
% Description: evaluate the simulation results of the model from the
% aspects of walking speed, metabolic expenditure, joint hyperextension,
% etc.
% Author(s): Jin Wei, Peking U. wjin24@stu.pku.edu.cn
%
% [fit, parts] = measure_simResults(modelInfo) also returns the weighted
% per-term contributions (parts.completeness, parts.velocity, parts.knee_limit,
% parts.ankle_limit, parts.grf_limit, parts.effort).  Their sum is fit; use
% this to diagnose a suspicious score (e.g. a fall that still looks "good").
if ~modelInfo.dy.hasRun
    error('simulation did not operate normally.')
end

measureObj = Measurements(modelInfo);

% if modelInfo.stage == 1
% fit = measureObj.gait_completeness_measure();

% elseif modelInfo.stage == 2
% weight of each term
weight_time = 1;
weight_vel = 100;
weight_knee = 0.01;
weight_ankle = 0.01;
weight_grf = 10;
weight_effort = 1;

% calculate fitness (per-term values are kept for diagnostics)
parts = struct();
parts.completeness = measureObj.gait_completeness_measure(weight_time);
parts.velocity     = measureObj.gait_velocity_measure(weight_vel);
parts.knee_limit   = measureObj.knee_limit_measure(weight_knee);
parts.ankle_limit  = measureObj.ankle_limit_measure(weight_ankle);
parts.grf_limit    = measureObj.grf_limit_measure(weight_grf);
parts.effort       = measureObj.effort_measure(weight_effort);

fit = parts.completeness + parts.velocity + parts.knee_limit ...
    + parts.ankle_limit + parts.grf_limit + parts.effort;

% Defensive check: every term must be non-negative, otherwise the objective
% can be gamed.  Regression: dividing the effort by the SIGNED pelvis
% displacement made a 1.9 s backward fall score -55 (effort = -221), which
% CMA-ES read as better than walking.
termNames = {'completeness', 'velocity', 'knee_limit', ...
             'ankle_limit', 'grf_limit', 'effort'};
negSum = 0;
for it = 1:numel(termNames)
    negSum = negSum + min(0, parts.(termNames{it}));
end
if negSum < -1e-9
    warning('measure_simResults:NegativeTerm', ...
        'A fitness term is negative (sum = %.3f); the objective may be exploitable.', ...
        negSum);
end

parts.total = fit;
parts.pelvisX = measureObj.distance;
parts.lastTime = measureObj.lastTime;

% else
% error('[CUSTOMIZED ERROR] invalid optimization stage information.')

% end

end