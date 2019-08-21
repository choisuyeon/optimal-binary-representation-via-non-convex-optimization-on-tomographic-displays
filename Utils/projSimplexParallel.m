function x = projSimplexParallel(y)
% project an n-dim vector y to the simplex Dn
% Dn = { x : x n-dim, 1 >= x >= 0, sum(x) = 1}

% (c) Xiaojing Ye
% xyex19@gmail.com
%
% Algorithm is explained as in the linked document
% http://arxiv.org/abs/1101.6081
% or
% http://ufdc.ufl.edu/IR00000353/
% Jan. 14, 2011.
%
%
% Suyeon Choi(suyeon@stanford.edu)
% modified it to matrix version (for parallel use)
%
% Oct. 21, 2018.

numRows = size(y, 1);
bget    = false(numRows, 1);
numCol  = size(y, 2); 
s       = sort(y, 2, 'descend'); 
tempSum = 0;

for col = 1 : numCol-1
    tempSum = tempSum + s(:, col);
    tmax    = (tempSum - 1) / col;
    bget    = tmax >= s(:, col+1);
end
    
tmax = bget .* tmax + (1 - bget) .* (tempSum + s(:, numCol) - 1) / numCol;
x = max(y - repmat(tmax, 1, numCol), 0);
return;