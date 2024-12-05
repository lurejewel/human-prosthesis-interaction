function mat = vec_2_mat(vec)
% Name: vec2mat.m
% Description: Convert OpenSim Vector to MATLAB 1xn matrix.

mat = nan(vec.size, 1);
for i = 1 : vec.size
    mat(i) = vec.get(i-1);
end

end

