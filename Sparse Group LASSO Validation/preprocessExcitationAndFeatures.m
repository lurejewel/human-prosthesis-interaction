function [excitationOut, featureMatrixOut] = preprocessExcitationAndFeatures(excitationIn, featureMatrixIn, delaySteps, excitationThreshold)
% preprocessExcitationAndFeatures
% Apply excitation/feature alignment and sample filtering:
% 1) shift excitation forward by delaySteps and trim feature rows to match
% 2) remove samples with excitation <= excitationThreshold
%
% Inputs:
%   excitationIn       : T x 1 excitation vector
%   featureMatrixIn    : T x N feature matrix
%   delaySteps         : nonnegative integer (default 1)
%   excitationThreshold: scalar threshold (default 0.01)
%
% Outputs:
%   excitationOut    : processed excitation vector
%   featureMatrixOut : processed feature matrix aligned with excitationOut

if nargin < 3 || isempty(delaySteps)
    delaySteps = 1;
end
if nargin < 4 || isempty(excitationThreshold)
    excitationThreshold = 0.01;
end

if ~isvector(excitationIn)
    error('excitationIn must be a vector.');
end
excitationIn = excitationIn(:);

if ~ismatrix(featureMatrixIn)
    error('featureMatrixIn must be a 2D matrix.');
end
if size(featureMatrixIn, 1) ~= numel(excitationIn)
    error('excitationIn and featureMatrixIn must have the same number of rows.');
end

if ~(isscalar(delaySteps) && isnumeric(delaySteps) && isfinite(delaySteps) && delaySteps >= 0 && floor(delaySteps) == delaySteps)
    error('delaySteps must be a nonnegative integer scalar.');
end
if delaySteps >= numel(excitationIn)
    error('delaySteps (%d) must be smaller than signal length (%d).', delaySteps, numel(excitationIn));
end

if ~(isscalar(excitationThreshold) && isnumeric(excitationThreshold) && isfinite(excitationThreshold))
    error('excitationThreshold must be a finite numeric scalar.');
end

% Align by neuromuscular transmission delay.
excitationOut = excitationIn((delaySteps + 1):end);
featureMatrixOut = featureMatrixIn(1:(end - delaySteps), :);

% Remove near-zero excitation samples and keep rows aligned.
keepMask = excitationOut > excitationThreshold;
excitationOut = excitationOut(keepMask);
featureMatrixOut = featureMatrixOut(keepMask, :);
end
