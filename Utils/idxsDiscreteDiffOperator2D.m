function [i, j, v] = idxsDiscreteDiffOperator2D(n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    returns indices for a sparse discrete differential operator 
%    anistropic, isotropic
%
%    see 
%    A WEIGHTED DIFFERENCE OF ANISOTROPIC AND ISOTROPIC TOTAL VARIATION MODEL FOR IMAGE PROCESSING
%    for more details.
%
%    Contact:
%       Suyeon Choi (suyeon@stanford.edu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % x axis diff
    i1 = [1:(n - 1) * n 1:(n - 1) * n];
    j1 = zeros(1, (n - 1) * n);
    for x = 1:n
        for y = 1:n - 1
            idx = y + (x - 1) * (n - 1);
            j1(idx) = y + (x - 1) * n;                
        end
    end
    j1 = [j1 j1 + 1];

    i2 = i1 + (n - 1) * n;
    j2 = zeros(1, (n - 1) * n);
    for x = 1:n
        for y = 1:n - 1
            idx = y + (x - 1) * (n - 1);
            j2(idx) = x + (y - 1) * n;                
        end
    end   
    j2 = [j2 j2 + n];     

    i = [i1 i2];
    j = [j1 j2];
    v1 = ones(1, (n - 1) * n);
    v2 = -1 * ones(1, (n - 1) * n);
    v = [v1 v2 v1 v2];
    % y axis diff
%     D = sparse(i, j, v, 2 * (n - 1) * n, n * n);

end