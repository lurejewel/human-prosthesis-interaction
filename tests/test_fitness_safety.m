% -------------------------------------------------------------------------
% Name: test_fitness_safety.m
% Description: Regression test for the effort term (metabolic cost of
%   transport).  A model that falls early while drifting BACKWARD must be
%   charged a positive, large cost -- not a negative one.
%
%   Motivation: the run results/opt_result_2026-11-20_10-37-17.mat finished
%   with fitness -55 although the model fell at t = 1.9 s (pelvisX = -0.15 m).
%   The old effort term was
%       fit = weight * totalMetabolic / mass / pelvisX      % signed!
%   which gave -220.8 for that failed gait and cancelled the fall penalty
%   (65.3) plus the velocity-failure floor (100).
% -------------------------------------------------------------------------
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions', 'core'));

metab = 2461.3;   % [J]   total metabolic energy of the failed 1.9 s gait
mass  = 74.3;     % [kg]  body mass of that simulation

cFell   = Measurements.cost_of_transport(metab, mass, -0.15);
cStand  = Measurements.cost_of_transport(metab, mass, 0.0);
cShort  = Measurements.cost_of_transport(metab, mass, 0.2);
cWalk1  = Measurements.cost_of_transport(metab, mass, 1.0);
cWalk10 = Measurements.cost_of_transport(metab, mass, 10.0);

% non-positive displacement -> floored, never negative
assert(cFell  > 0, 'backward displacement produced a non-positive cost');
assert(cStand > 0, 'zero displacement produced a non-positive cost');
assert(abs(cFell  - metab / mass / 0.5) < 1e-9, 'backward case must use the floor');
assert(abs(cStand - metab / mass / 0.5) < 1e-9, 'zero case must use the floor');
% any positive displacement -> EXACTLY the old expression (no behaviour change)
assert(abs(cShort - metab / mass / 0.2) < 1e-9, 'positive displacement must be unchanged');
assert(abs(cWalk1  - metab / mass / 1.0) < 1e-9, 'positive displacement must be unchanged');
assert(cWalk1 > cWalk10, 'cost must decrease as forward distance grows');
assert(cWalk10 > 0, 'forward walking cost must be positive');

fprintf('test_fitness_safety: OK  (fell=%.2f, stand=%.2f, 1m=%.2f, 10m=%.2f)\n', ...
    cFell, cStand, cWalk1, cWalk10);
