% Sparse Group LASSO runner script.
% Problems:
% - 目前的代码未考虑reflex delay
% - 目前的beta参数未考虑L0, P0等；并非单纯的sigma = X * beta
% - 目前的数据是全步态周期的数据。后续需要根据步态阶段，对数据进行切割和重新拼接
%
% TODO:
% - 首先验证soleus_r的真实肌肉策略
% - 重塑sigma = X * beta公式，使之包含L0 P0
% - 研究SGL代码，看是否需要针对性改造（权重？惩罚？）
% - （soleus_r成功后）切分数据，对其他肌肉的控制律进行验证
% - （所有肌肉成功后）进行正动力学验证
% - （正动力学成功后）进行代理模型研究（对应Annuals Review of Biomedical Engineering中的大脑预测模型）

close, clear, clc

stoFilePath = 'UN.sto';
muscleName = 'soleus_r';
[y, X, XMap] = extractStoMuscleFeatures(stoFilePath, muscleName); % y: excitation of the muscle; X: feature matrix; XMap: feature map

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

[beta, intercept, betaStd, stats] = sgl_fit(X, y, lambda1, lambda2, opts);

fprintf('\nScript run complete. Final objective = %.12g\n', stats.objHist(end));
fprintf('Converged = %d, iterations = %d\n', stats.converged, stats.iters);
fprintf('Intercept = %.12g\n', intercept);

disp('beta (original scale):');
disp(beta);
disp('betaStd (standardized scale):');
disp(betaStd);
