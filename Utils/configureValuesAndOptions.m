%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    a script for initialization
%    + non-convex optimization configuration
%
%    Contact:
%       Suyeon Choi (suyeon@stanford.edu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Section 4.1 Display update configuration
% load initial U and V; initialize with center pupil view
if strcmp(display.initialCondition, "RGB image")
    numCV = '0025';                                % center of 7x7 viewpoints
    centerViewpointImage = double(imresize(imread(sprintf('./Contents/Target/%s/Pupilviews/%s.png', system.content, numCV)), display.resolution)) ./ 255;
    U = centerViewpointImage(:);
else
    U = ones(display.numel * 3, 1) * display.initialCondition; % LCD
end
V = zeros(size(b));

%% Section 4.2 Binary Optimization Configuration
% load initial X and Y;
L = linspace(0, 1, 2);                              % declare a discrete set
loadInitialX;                                       % load initial X
loadInitialY;                                       % load initial Y

%%% nonconvex : alpha, etc hyperparameters
power        = 2;                   
alphaInitial = 0.1 + 1e-12;         
alphaEnd     = 80;                  

i = 1:binaryMaxIterMask;
alphas                      = alphaEnd * ones(maxIter, 1);
alphas(1:binaryMaxIterMask) = alphaInitial + (alphaEnd - alphaInitial) * ((i-1) / (binaryMaxIterMask-1)).^power;

lambdas = ones(maxIter, 1) * 3;     % set lambda 3 constant;

proxF_MaxIter = 200;                %% you may change this
proxF_Threshold = 1e-5;             %% you may change this
opts.epsilon = proxF_Threshold;                       % threshold for primal-dual gap
opts.maxIter = proxF_MaxIter;                         % max iteration number if never exceed threshold

opts.checkIter = 10;                %% you may change this; frequency to check    
opts.gamma1 = 1;                                      % coefficient for lipschitz constant
opts.gamma2 = 1;                                      % coefficient for lipschitz constant
opts.discreteSet = L;
opts.numLayers = system.numLayers;
opts.b = b; 
opts.illumination_factor = t;                         % whether to regularize
opts.useGPU = flag.useGPU;

if lambdas(1) > 0 && strcmp(methodBinary, 'PALM+PD')
    opts.D = discreteDiffOperatorMultilayer(display.resolution(1), system.numLayers);
end
