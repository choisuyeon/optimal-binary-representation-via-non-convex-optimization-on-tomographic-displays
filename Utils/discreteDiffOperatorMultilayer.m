function D = discreteDiffOperatorMultilayer(n, M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Discrete differential operator for multi-layer images
%    anistropic, isotropic
%
%    M : number of planes
%
%    Contact:
%       Suyeon Choi (suyeon@stanford.edu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('   - Creating a matrix for TV..');
    tic;
    [i, j, v] = idxsDiscreteDiffOperator2D(n);
    I = zeros(1, M * length(i));
    J = zeros(1, M * length(j));
    V = zeros(1, M * length(v));
    for m = 0:M-1
        I(1 +    m    * length(i) : ...
              (m + 1) * length(i)) = i + m * 2 * (n - 1) * n;
        J(1 +    m    * length(j) : ...
              (m + 1) * length(j)) = j + m * n^2;
        V(1 +    m    * length(v) : ...
              (m + 1) * length(v)) = v;
    end
    
    D = sparse(I, J, V, M * 2 * (n - 1) * n, M * n^2);
    elapsedTime = toc;
    fprintf('   - Created the matrix... : %.2fs \n', elapsedTime);
end