function M = buildM(n)
%BUILDM Construct the binary mapping matrix M for augmented linear terms.
%   M has size (2*n+9)-by-(3*n+12) and contains only 0/1 entries.

    validateattributes(n, {'numeric'}, {'scalar', 'real', 'finite', 'positive', 'integer'}, mfilename, 'n', 1);

    numRows = 2*n + 9;
    numCols = 3*n + 12;

    % Preallocate index lists for all nonzero entries.
    nnzCount = n + n + 8 + n + 4;
    I = zeros(nnzCount, 1);
    J = zeros(nnzCount, 1);
    k = 0;

    % Block A, Segment I: row i -> odd column (2*i-1).
    rowsA = (1:n).';
    colsA = (2*(1:n) - 1).';
    idx = k + (1:n);
    I(idx) = rowsA;
    J(idx) = colsA;
    k = k + n;

    % Block B, Segment II: identity matrix.
    rowsB = (n + (1:n)).';
    colsB = (2*n + (1:n)).';
    idx = k + (1:n);
    I(idx) = rowsB;
    J(idx) = colsB;
    k = k + n;

    % Block C, Segment III: triplet endpoints for four groups.
    rowsC = (2*n + (1:8)).';
    colsC_local = [1; 3; 4; 6; 7; 9; 10; 12];
    colsC = 3*n + colsC_local;
    idx = k + (1:8);
    I(idx) = rowsC;
    J(idx) = colsC;
    k = k + 8;

    % Block D, Segment I: even columns only.
    rowD = 2*n + 9;
    colsD_I = (2:2:2*n).';
    idx = k + (1:n);
    I(idx) = rowD;
    J(idx) = colsD_I;
    k = k + n;

    % Block D, Segment III: middle element of each triplet.
    colsD_III = (3*n + [2; 5; 8; 11]);
    idx = k + (1:4);
    I(idx) = rowD;
    J(idx) = colsD_III;

    M = sparse(I, J, 1, numRows, numCols);

    % Optional safety checks from the specification.
    assert(all(M(:) == 0 | M(:) == 1), 'M must be binary.');
    assert(isequal(size(M), [2*n+9, 3*n+12]), 'Size mismatch.');
end
