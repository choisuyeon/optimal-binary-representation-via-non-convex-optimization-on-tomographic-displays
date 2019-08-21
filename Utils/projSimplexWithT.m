function result = projSimplexWithT(x, numCols, t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    projection on simplex considering brightness
%
%    Contact:
%       Suyeon Choi (suyeon@stanford.edu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % check numel(x) / num_col is integer
   
   result = reshape(x, numel(x) / numCols, numCols);
   result = t * projsplxParallel(result / t);
   result = result(:);

end