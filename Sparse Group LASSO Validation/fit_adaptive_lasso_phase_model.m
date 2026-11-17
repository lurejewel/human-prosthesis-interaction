function [fit, lasso] = fit_adaptive_lasso_phase_model(stoFilePath, options)
%FIT_ADAPTIVE_LASSO_PHASE_MODEL  Adaptive LASSO muscle-reflex model per gait phase.
%
%   按步态相位拟合 adaptive LASSO 肌肉反射模型（单次、确定性拟合）。
%
%   [fit, lasso] = fit_adaptive_lasso_phase_model(stoFilePath, options)
%
%   For every (muscle, phase group) this estimates
%
%       y_m(t) = X(t) * beta_{g,m} + b_{g,m} + noise
%
%   where X is the 26-feature matrix (9 lCEN + 9 fATN + 8 joint kinematics) and
%   y is the excitation target, using the censored adaptive-LASSO pipeline:
%   weighted (Huber) OLS -> adaptive weights -> lambda path with per-lambda
%   post-selection refit -> lambda choice -> Huber IRLS on the selected support.
%
%   This function is now a thin wrapper: the data-dependent half lives in
%   load_lasso_design_matrices.m and the fitting half in
%   fit_adaptive_lasso_core.m, so that a multi-seed bootstrap ensemble can
%   reuse one expensive STO read (see lasso_stability_identification.m).
%   Signature, defaults and outputs are unchanged.
%
%   Input
%   -----
%   stoFilePath : LASSO-ready .sto.  Column naming as produced by
%                 write_lasso_dataset_sto.m, or any SCONE-style file with
%                 <muscle>_r.activation / .excitation / .fiber_length_norm /
%                 .mtu_force_norm and leg1_r.state.
%   options     : struct, defaults in lasso_fit_default_options.m:
%                 gamma 1, epsilonOLS 1e-8, nzTol 1e-5, lambda1 0,
%                 lambdaPathConfig.nLambda 60, lambdaMinRatio 1e-4,
%                 tolSelect 0.05, c 0.01, rho 50, delaySteps 1,
%                 excitationThreshold 0.01, sampleWeight [],
%                 robustInitIter 2, robustIter 2, huberDeltaFactor 1.345,
%                 verbose true, groupLabels {'stance','swing'},
%                 groupPhases {[0 1],[2 3 4]}, source 'MocoInverse',
%                 preprocessing struct().
%
%   Output
%   ------
%   fit   : full fit result (Beta_refit / Intercept_refit / YPred_refit /
%           RMSE_refit / LambdaPath / S / options ...).  See
%           fit_adaptive_lasso_core.m for the complete field list.
%   lasso : grouped-v2 controller ready for load_lasso_reflex_controller.m.
%
%   See also load_lasso_design_matrices, fit_adaptive_lasso_core,
%            lasso_stability_identification.
%
%   数值/实现注意
%   -------------
%   * c 是 censoring floor，必须与接口 cleanOpts.floor 一致；
%   * delaySteps 在加载阶段把 target 前移、features 后移，相位索引同样平移；
%   * lambda path 的每个点都做 refit，计算量较大，这是旧流程的设定；
%   * beta 在 |beta| <= nzTol 处被强制置零，保证 mask 与结构零严格一致。

if nargin < 1 || isempty(stoFilePath)
    error('fit_adaptive_lasso_phase_model:badInput', 'stoFilePath is required.');
end
if nargin < 2 || isempty(options)
    options = struct();
end

opts = lasso_fit_default_options(options);
D = load_lasso_design_matrices(stoFilePath, opts);
[fit, lasso] = fit_adaptive_lasso_core(D, opts);
end
