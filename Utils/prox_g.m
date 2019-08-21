function result = prox_g(q, sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Proximal operator with g and 1/sigma in the paper
%    Using projection on simplex (projection each row on simplex)
%
%    Input 
%         1 : variable of proximal operator
%         1 / sigma : coefficient of norm-2 term (non-used)
%
%    Contact:
%       Suyeon Choi (suyeon@stanford.edu)
%               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%% non-parallel version
N = size(q, 1);
%result = q;
for i = 1 : N
    i
    result(i, :) = projsplx(q(i, :));
end
%}
%% parallel version (Suyeon)
result = projSimplexParallel(q);
end

