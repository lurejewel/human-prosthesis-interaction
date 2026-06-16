function fit = measure_simResults(modelInfo)
% Name: measure_simResults
% Description: evaluate the simulation results of the model from the
% aspects of walking speed, metabolic expenditure, joint hyperextension,
% etc.
% Author(s): Jin Wei, Peking U. wjin24@stu.pku.edu.cn
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

% calculate fitness
fit = measureObj.gait_completeness_measure(weight_time) ...
    + measureObj.gait_velocity_measure(weight_vel) ...
    + measureObj.knee_limit_measure(weight_knee) ...
    + measureObj.ankle_limit_measure(weight_ankle) ...
    + measureObj.grf_limit_measure(weight_grf) ...
    + measureObj.effort_measure(weight_effort);

% else
% error('[CUSTOMIZED ERROR] invalid optimization stage information.')

% end

end