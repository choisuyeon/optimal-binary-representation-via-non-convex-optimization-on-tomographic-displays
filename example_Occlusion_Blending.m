%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    optimization code for tomographic displays (https://doi.org/10.1038/s41467-019-10451-2)
%    you may find https://doi.org/10.1364/OE.27.024362 for more implementation details.
%
%    It solves multi view-based optimization problem for tomographic
%    display systems.
%
%    Contact:
%       Suyeon Choi (suyeon@stanford.edu)
%       Seungjae Lee (seungjae1012@gmail.com)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

%% Overall optimization configuration
maxIter           = 100;       
binaryMaxIter     = 90;       
binaryMaxIterMask = 90;     
periodEPR         = 30;  

methodBinary      = 'SART';       % 'SART' or 'PALM+PD' or 'None'
methodContinuous  = 'SART';       % 'SART' or 'PD' or 'None'

fprintf('  бс Optimizing a display and the DMD backlight pattern using %s and %s, respectively... \n', methodContinuous, methodBinary);

%% System, display and backlight configuration
rad = pi / 180;
err = 1e-6;

flag.onlyCreateProjectionMatrix = 0;
flag.useGPU                     = 0;               % available only when M < 10
flag.intermediateEPR            = 1;               % use intermediate rounding
flag.cropTarget                 = 0;

system.content    = 'Road';     % б┌
system.DOF        = [0 5.5];
system.brightness = 0.3;        % beta (0.1 ~ 0.4)
system.numLayers  = 80;         % M
t = round(system.brightness * system.numLayers);

loadSceneConfiguration;

display.numel            = prod(display.resolution);
display.range            = [0 1];
display.initialCondition = 'RGB image';  
                                     %  options for display initial:
                                     % 1. 'RGB image' : pupil sample from
                                     %                  the center viewpoint
                                     % 2. real number between [0 1]

backlight.DC               = 0;      % bg values
backlight.range            = [backlight.DC / (1 + backlight.DC) 1];
backlight.resolution       = display.resolution;
backlight.initialCondition = 'binary_blending';   
                                     %  options for backlight initial:
                                     % 1. 'binary_blending'
                                     % 2. 'considerBrightness (t)'
                                     % 3. real number between [0 1]

%% load the projection matrix (A)
          
if  ~exist('T', 'var')
    % TODO : check other specifications (e.g. disp.resolution)
    disp('   - Creating a projection matrix..')
    [T, T_ori, ~] = focalPlanes2pupilViews(pupil, system, display);
else    
    disp('   - The projection matrix already exists..')
end

%% load target images (b)
targetDir = sprintf('%s/Contents/Target/%s/Pupilviews', pwd, system.content);
b = loadTargetImages(pupil, targetDir);
b = b * (t / (1 + backlight.DC) + ...
                      backlight.DC / (1 + backlight.DC) * system.numLayers); 

%% initialize & set options for non-convex optimization
configureValuesAndOptions;

%% make a bin
dirData = sprintf('./Results/latest_[%s+%s]', ...
                       methodBinary, ...
                       methodContinuous);
mkdir('./Results/');
mkdir(dirData);

SLMlayer    = reshape(U, display.resolution(1),   display.resolution(2), 3);
BLlayer     = reshape(X, backlight.resolution(1), backlight.resolution(2) * system.numLayers);
imwrite(SLMlayer, sprintf('./%s/SLM0.bmp', dirData));
imwrite(BLlayer,  sprintf('./%s/BL0.bmp',  dirData));

%% Optimal binary representation via non-convex optimization on tomographic displays (Opt. Express 2019)
elapsedTime = tic;

for j = 1:maxIter       
    disp(j)    
        
    % binary update (with fixed SLM) (Section 4.2)
    if j <= binaryMaxIter && ~strcmp(methodBinary, 'None')
        
        % calculating K = A * U
        opts.K   = T * sparse(1 : 3*display.numel * system.numLayers, ...
                              repmat(1 : display.numel * system.numLayers, 1, 3), ...
                             [repmat(U(                  1 :   display.numel)', 1, system.numLayers) ...
                              repmat(U(  display.numel + 1 : 2*display.numel)', 1, system.numLayers) ...
                              repmat(U(2*display.numel + 1 : 3*display.numel)', 1, system.numLayers)], ...
                              3*display.numel * system.numLayers, ...
                                display.numel * system.numLayers);

        [X, Y]   = binaryUpdate(X, Y, methodBinary, opts);
    end
    
    % EPR (Section 4.4) 
    if j == binaryMaxIter || ...
       (j < binaryMaxIter && flag.intermediateEPR && mod(j, periodEPR) == 0)       
        [X, energyCompensated] = EPR(X, L, backlight.resolution, system.numLayers);
        if strcmp(methodBinary, 'PALM+PD')
            Y = yFromX(X, L);  
        end
    end     
    
    % SLM update (with fixed backlight) (Section 4.1)
    if j > 0 && ~strcmp(methodContinuous, 'None')
        
        % Calculating K = A * X....
        if j <= binaryMaxIter
            opts.K = T * sparse(1 : (3 * display.numel * system.numLayers), ...
                                [repmat(1 : display.numel, 1, system.numLayers) ...
                                 repmat(1 : display.numel, 1, system.numLayers) + display.numel ...
                                 repmat(1 : display.numel, 1, system.numLayers) + 2*display.numel], ...
                                 repmat(X', 1, 3), ...
                                3*display.numel * system.numLayers, ...
                                3*display.numel);
        end
        [U, V]   = displayUpdate(U, V, methodContinuous, opts);      
    end
    
    SLMlayer    = reshape(U, display.resolution(1),   display.resolution(2), 3);
    BLlayer     = reshape(X, backlight.resolution(1), backlight.resolution(2) * system.numLayers);
    imwrite(SLMlayer, sprintf('./%s/SLM%d.bmp', dirData, j));
    imwrite(BLlayer, sprintf('./%s/BL%d.bmp', dirData, j));             
end
elapsedTime = toc(elapsedTime)

clearvars D opts;
