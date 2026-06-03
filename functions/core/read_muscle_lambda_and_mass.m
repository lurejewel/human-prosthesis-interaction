function [muscleLambda, muscleMass, muscleMap] = read_muscle_lambda_and_mass(filename)

% 读取 CSV 文件
T = readtable(filename, 'Delimiter', ',');

% 肌肉数量
nMuscles = height(T);

muscleLambda = nan(nMuscles,1); % initialize
muscleMass= nan(nMuscles,1); % initialize
muscleMap = containers.Map('KeyType', 'char', 'ValueType', 'int32');

for i = 1 : nMuscles
    muscleLambda(i) = T.type_i_ratio(i);
    muscleMass(i) = T.mass(i);
    muscleMap(T.abbr{i}) = i;
end


% % 预分配结构体数组
% muscleProp(nMuscles) = struct( ...
%     'abbr',   '', ...
%     'name',   '', ...
%     'lambda', []); % , ...
%     % 'fOpt',   [], ...
%     % 'lOpt',   []  );
% 
% % 创建 Map：key = abbr, value = index
% muscleMap = containers.Map('KeyType', 'char', 'ValueType', 'int32');
% 
% % 填充数据
% for i = 1:nMuscles
%     muscleProp(i).abbr   = T.abbr{i};
%     muscleProp(i).name   = T.muscle_name{i};
%     muscleProp(i).lambda = T.type_i_ratio(i);
%     % muscleProp(i).fOpt   = T.optimal_force(i);
%     % muscleProp(i).lOpt   = T.optimal_fiber_length(i);
% 
%     muscleMap(muscleProp(i).abbr) = i;
% end


end