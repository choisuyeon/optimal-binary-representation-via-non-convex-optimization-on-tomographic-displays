function Y = yFromX(X, L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    It returns Y from the given X. 
%    (can be used in the initialization)
%    in non-convex optimization for Tomographic Displays
%
%    parameters
%    X : column vector;
%    L : predefined discrete set
%
%    Contact:
%       Suyeon Choi (suyeon@stanford.edu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    X = X(:)';
    numRows = length(X);
    numCols = length(L);
    Y = zeros(numRows, numCols);
    
    % matrix version : efficient implementation
    % TODO
    
    % loop version : unefficient but intuitive and easy implementation..
    tic;
    for i = 1:numRows
        u = X(i);
        [member, j] = ismember(u, L);
        if member
            Y(i, j) = 1;
        else
            for l = 1:numCols
               z = 0;
               for j = 1:numCols
                  z = z + (u - L(l))^2 / (u - L(j))^2;
               end
               z = 1/z;
               Y(i, l) = z;
            end            
        end
    end
    elapsedTime = toc;
    fprintf('   - Y from X is calculated... elasped time : %.1fs\n', elapsedTime);
end

