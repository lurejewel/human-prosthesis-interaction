function vec = mat_2_vec(mat)
% Name: mat2vec.m
% Description: Convert MATLAB 1xn matrix to OpenSim Vector.

vec = org.opensim.modeling.Vector(length(mat), nan);
for i = 0 : length(mat)-1
    vec.set(i, mat(i+1));
end

end

