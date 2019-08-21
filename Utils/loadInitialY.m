%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    a script for loading initial matrix variable Y
%
%    Contact:
%       Suyeon Choi (suyeon@stanford.edu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(methodBinary, 'PALM+PD')
    Y = yFromX(X, L);                                   % load initial Y, predefined function
else
    Y = 0;
end
