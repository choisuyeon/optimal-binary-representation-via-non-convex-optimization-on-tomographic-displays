%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    a script for loading initial vectorized image X
%
%    Contact:
%       Suyeon Choi (suyeon@stanford.edu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(backlight, 'initialCondition')
    disp('   -!! set the initial condition of backlight !!');
    return
end

switch backlight.initialCondition
    case 'binary_blending'
        depthmapFilename = sprintf('./Contents/Target/%s/Depthmaps/%04d.%s', ...
                                    system.content, ...
                                    ceil(prod(pupil.numViewpoints) / 2), ...
                                    target.extensionBL);
        depthmap = double(imresize(imread(depthmapFilename), backlight.resolution, 'nearest'))/255;
        X = depthmap2Backlight(depthmap, ...
                               system.numLayers, ...
                               backlight.resolution, ...
                               [target.DOF(1) + err, target.DOF(2)], ...
                               system.DOF, ...
                               t);
                           
    case 'considerBrightness'
        X = ones(backlight.resolution(1), backlight.resolution(2) * system.numLayers) * system.brightness;
        
    otherwise
        X = ones(backlight.resolution(1), backlight.resolution(2) * system.numLayers) * backlight.initialCondition;
        
end

X = X(:);