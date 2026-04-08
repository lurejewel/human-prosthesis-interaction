% Sparse Group LASSO runner script.
% Problems:
% - soleus_r在swing阶段的优化效果极差，可能需要增加excitation<0.01时的不等式约束
% - 所有的截距可能都不需要 -> 缩减mask matrix和beta的尺度
%
% TODO:
% - maskMatrix需要根据肌肉特性，导入先验（重读经典论文，寻找先验依据）
% - 研究SGL代码，看是否需要针对性改造（权重？惩罚？）
% - （soleus_r成功后）切分数据（stance-swing），对其他肌肉的控制律进行验证
% - （所有肌肉成功后）进行正动力学验证
% - （正动力学成功后）进行代理模型研究（对应Annuals Review of Biomedical Engineering中的大脑预测模型）
% - （用于论文）lambda1 lambda2的灵敏度分析

close, clear, clc

stoFilePath = 'UN.sto';
muscleName = 'soleus_r';
[yAll, XAll, featureToIndexMap, M, phase] = extractStoMuscleFeatures(stoFilePath, muscleName); % yAll: excitation of the muscle; XAll: feature matrix; featureToIndexMap: feature->column map; M: mask matrix (composed of 0 and 1)
% yAllStance = yAll(phase.stanceIdx, :);
% yAllSwing = yAll(phase.swingIdx, :);
% XAllStance = XAll(phase.stanceIdx, :);
% XAllSwing = XAll(phase.swingIdx, :);
% [y, X] = preprocessExcitationAndFeatures(yAllSwing, XAllSwing, 1, 0.01); % 1-step delay and remove excitation <= 0.01
[y, X] = preprocessExcitationAndFeatures(yAll, XAll, 1, 0.01); % 1-step delay and remove excitation <= 0.01

lambda1 = 0.5;
lambda2 = 2;

opts = struct();
opts.maxIterOuter = 2000;
opts.maxIterInner = 2000;
opts.tolOuter = 1e-6;
opts.tolInner = 1e-6;
opts.tol1D = 1e-10;
opts.maxIter1D = 50;
opts.verbose = true;

A = X * M; % feature matrix X (augmented), masked with the sparse matrix M
[beta, intercept, betaStd, stats] = sgl_fit(A, y, lambda1, lambda2, opts);

% test
% yPred = max(0.01, X * M * beta + intercept);
% figure, plot(y), hold on, plot(yPred)
% sqrt(mean((y - yPred).^2))
% beta(:) = 0; beta(26) = 0.68247; beta(27) = -0.58974; beta(37) = 0.44245; intercept = 0.097071;
% yPredAllSwing = XAllSwing * M * beta + intercept;
% figure, plot(yAllSwing), hold on, plot(yPredAllSwing)
% sqrt(mean((yAllSwing - yPredAllSwing).^2))


yPredAll = XAll * M * beta + intercept;
rmse = sqrt(mean((yAll - yPredAll).^2));

fprintf('Sparse Group LASSO algorithm run complete. Final objective = %.12g\n', stats.objHist(end));
fprintf('RMSE: %.6g\n', rmse);
% fprintf('Converged = %d, iterations = %d\n', stats.converged, stats.iters);

%% results

fprintf('\nNonzero muscle reflex controls:\n');
indexToFeatureName = buildInverseFeatureMap(featureToIndexMap, size(X, 2));
nzTol = 1e-5;
nzIdx = find(abs(beta) > nzTol);
if isempty(nzIdx)
    fprintf('(none)\n');
else
    for k = 1:numel(nzIdx)
        i = nzIdx(k);
        featureName = indexToFeatureName{i};
        fprintf('%s -> %s: %.5g\n', featureName, muscleName, beta(i));
    end
end
fprintf('prestimulation: %.5g\n', intercept);

figure, plot(yAll), hold on
plot(max(0.01, yPredAll))
xlim([0 numel(yAll)])
xlabel('#timesteps'), ylabel('excitation')
legend('SCONE muscle excitation', 'SGL fitting results')
title('SGL fitting results, compared with SCONE muscle excitation')

%% functions
function indexToFeatureName = buildInverseFeatureMap(featureToIndexMap, numFeatures)
ks = keys(featureToIndexMap);
vs = cell2mat(values(featureToIndexMap));

if numel(ks) ~= numFeatures
    error('buildInverseFeatureMap:CountMismatch', ...
        'feature map size (%d) does not match numFeatures (%d).', numel(ks), numFeatures);
end

if any(vs < 1) || any(vs > numFeatures) || any(abs(vs - round(vs)) > 0)
    error('buildInverseFeatureMap:BadIndex', 'feature indices must be integers in [1, %d].', numFeatures);
end

if numel(unique(vs)) ~= numel(vs)
    error('buildInverseFeatureMap:NonBijective', 'feature map is not one-to-one.');
end

indexToFeatureName = cell(numFeatures, 1);
for t = 1:numel(ks)
    indexToFeatureName{vs(t)} = ks{t};
end
end