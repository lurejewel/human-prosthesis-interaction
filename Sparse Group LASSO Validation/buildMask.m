function M = buildMask(n)
%BUILDM Construct an identity mapping matrix for linear terms.
%   M has size (2*n+8)-by-(2*n+8) and equals eye(2*n+8).

    validateattributes(n, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}, mfilename, 'n', 1);

    dim = 2*n + 8;
    M = speye(dim);

    % Safety checks.
    assert(all(M(:) == 0 | M(:) == 1), 'M must be binary.');
    assert(isequal(size(M), [dim, dim]), 'Size mismatch.');
end
