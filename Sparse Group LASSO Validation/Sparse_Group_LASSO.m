% Sparse Group LASSO runner script.
% Problems:
% - 目前的数据是全步态周期的数据。后续需要根据步态阶段，对数据进行切割和重新拼接
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
[yAll, XAll, XMap, M] = extractStoMuscleFeatures(stoFilePath, muscleName); % yAll: excitation of the muscle; XAll: feature matrix; XMap: feature map; M: mask matrix (composed of 0 and 1)
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

A = [X, ones(size(X,1),1)] * M; % feature matrix X (augmented), masked with the sparse matrix M
[beta, intercept, betaStd, stats] = sgl_fit(A, y, lambda1, lambda2, opts);

yPredAll = [XAll, ones(size(XAll,1),1)] * M * beta + intercept;
rmse = sqrt(mean((yAll - yPredAll).^2));

fprintf('\Sparse Group LASSO algorithm run complete. Final objective = %.12g\n', stats.objHist(end));
fprintf('RMSE: %.6g\n', rmse);
% fprintf('Converged = %d, iterations = %d\n', stats.converged, stats.iters);

fprintf('\nNonzero muscle reflex controls:\n');
nzTol = 1e-5;
nzIdx = find(abs(beta) > nzTol);
if isempty(nzIdx)
	fprintf('(none)\n');
else
	for k = 1:numel(nzIdx)
		i = nzIdx(k);
		[sourceName, param] = beta_index_to_fields(i, numel(beta), XMap, muscleName);
		fprintf('%s from %s to %s: %.5g\n', param, sourceName, muscleName, beta(i));
	end
end
fprintf('prestimulation: %.5g\n', intercept);

figure, plot(yAll), hold on
plot(max(0.01, yPredAll))
xlim([0 numel(yAll)])
xlabel('#timesteps'), ylabel('excitation')
legend('SCONE muscle excitation', 'SGL fitting results')
title('SGL fitting results, compared with SCONE muscle excitation')

function key = reverseMap(mp, val)
	ks = keys(mp);
	vs = values(mp);
	idx = find([vs{:}] == val, 1);
	if isempty(idx)
		error('reverseMap:KeyNotFound', 'Cannot find value %d in XMap.', val);
	end
	key = ks{idx};
	dotPos = strfind(key, '.');
	if ~isempty(dotPos)
		key = key(1:dotPos(1)-1);
	end
end

function [sourceName, param] = beta_index_to_fields(i, p, XMap, muscleName)
	nMus = (p - 12) / 3;
	if nMus < 1 || abs(nMus - round(nMus)) > eps(max(1, nMus))
		error('beta_index_to_fields:InvalidLength', ...
			'beta length must satisfy p = 3*#mus + 12. Got p = %d.', p);
	end
	nMus = round(nMus);

	if i < 1 || i > p
		error('beta_index_to_fields:IndexOutOfRange', 'Index %d out of range [1, %d].', i, p);
	end

	side = 'r';
	if endsWith(muscleName, '_l')
		side = 'l';
	elseif endsWith(muscleName, '_r')
		side = 'r';
	end

	if i <= 2 * nMus
		n = ceil(i / 2);
		sourceName = reverseMap(XMap, n);
		if mod(i, 2) == 1
			param = 'KL';
		else
			param = 'L0';
		end
	elseif i <= 3 * nMus
		n = i - 2 * nMus;
		sourceName = reverseMap(XMap, n);
		param = 'KF';
	else
		n = i - 3 * nMus;
		paramList = {'KP', 'P0', 'KV', 'KP', 'P0', 'KV', 'KP', 'P0', 'KV', 'KP', 'P0', 'KV'};
		param = paramList{n};

		if n <= 3
			sourceName = 'pelvis_tilt';
		elseif n <= 6
			sourceName = ['hip_flexion_' side];
		elseif n <= 9
			sourceName = ['knee_angle_' side];
		else
			sourceName = ['ankle_angle_' side];
		end
	end
end