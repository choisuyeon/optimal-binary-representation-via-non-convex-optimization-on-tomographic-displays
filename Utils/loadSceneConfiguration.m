%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    a script for loading scene spec; load predefined scene configurations   
%
%    Contact:
%       Suyeon Choi (suyeon@stanford.edu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rad = pi / 180;

if ~isfield(system, 'content')
    disp('   -!! set the content of system !!');
    return
end

fprintf('   - Loading scene configuration... : %s\n', system.content);

pupil.size = [6 6];                  % mm
target.extensionBL = 'tiff';

switch system.content
    case 'Buildings'
        system.FOV = 10 * rad; % Horizontal FOV
        target.DOF = [0 5.5];
        display.resolution = [400 400];
    case 'Road'
        system.FOV = 30 * rad; % Horizontal FOV
        target.DOF = [0 5.5];
        pupil.size = [6 6] * 6 / 7;
        display.resolution = [400 400];
        target.extensionBL = 'png';
    case 'Room'
        system.FOV = 20 * rad; % Horizontal FOV
        target.DOF = [0 4.0]; % CAUTION : this target depthmap is 0D : 4D = 0 : 255
        display.resolution = [350 350];
    case 'Room30'
        system.FOV = 30 * rad; % Horizontal FOV
        target.DOF = [0 5.5]; % CAUTION : this target depthmap is 0 : 5.5D = 0 : 255
        display.resolution = [450 450];
    case 'Castle'
        system.FOV = 30 * rad; % Horizontal FOV
        target.DOF = [0 5.5];
        display.resolution = [400 400];
    case 'Circles'
        system.FOV = 30 * rad; % Horizontal FOV
        target.DOF = [0 4.0];
        display.resolution = [350 350];
    case 'Ball'
        system.FOV = 30 * rad; % Horizontal FOV
        target.DOF = [0 4.0];
        display.resolution = [350 350];
    otherwise
        
end


pupil.numViewpoints = [7 7];         % 7 x 7 lf
pupil.resolution = display.resolution;
pupil.depth = 17;                    % arb    
