function [R, D] = lasso_stability_identification(stoFilePath, options)
%LASSO_STABILITY_IDENTIFICATION  Multi-seed consensus 0-1 structure from a .sto.
%
%   [R, D] = lasso_stability_identification(stoFilePath, options)
%
%   Convenience entry point: reads the LASSO-ready .sto ONCE
%   (load_lasso_design_matrices), then runs the whole multi-seed bootstrap
%   identification (lasso_stability_from_design) and writes
%
%     results/lasso_controller_result_moco_consensus.mat   (the controller)
%     results/lasso_stability_diagnostics.mat              (everything else)
%     results/lasso_stability_checkpoint.mat               (resume point)
%
%   options : see lasso_stability_default_options.  The identification
%             settings themselves (gamma, rho, c, lambda path, robust IRLS...)
%             live in options.fitOptions and are SHARED by all replicates.
%
%   Output R (see lasso_stability_from_design) and the design system D.
%
%   See also lasso_stability_from_design, lasso_stability_default_options,
%            Adaptive_LASSO_MocoInverse_Stability.

if nargin < 1 || isempty(stoFilePath)
    error('lasso_stability_identification:badInput', 'stoFilePath is required.');
end
if nargin < 2 || isempty(options)
    options = struct();
end

opts = lasso_stability_default_options(options);
opts.fitOptions = lasso_fit_default_options(opts.fitOptions);

% Read the .sto once; every replicate reuses these matrices.
D = load_lasso_design_matrices(stoFilePath, opts.fitOptions);

R = lasso_stability_from_design(D, opts);
end
