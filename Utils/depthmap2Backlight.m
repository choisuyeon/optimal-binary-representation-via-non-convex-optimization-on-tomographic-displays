function [BL, depthmap2Diopt] = depthmap2Backlight(depthMap, ...
                                                   numLayers, ...
                                                   sizeBL, ...
                                                   rangeDepthmap_Diopter, ...
                                                   rangeSystem_Diopter, ...
                                                   t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    It returns the initial backlight
%    from the given depthmap and the number of backlight layers.
%    ( = binary blending)
%
%    parameters
%
%    depthMap              : depth_map (in 2D matrix)
%    numLayers             : number of BLlayers
%    sizeBL                : size of BL (resolution of 1 layer)
%    rangeDepthmap_Diopter : Diopt correspoding to 0 and 1 of depthmap value
%    rangeSystem_Diopte  r : the furthest/nearest distance of layer in Diopt 
%    t                     : intended number of illumination per pixel
%
%    Contact:
%       Suyeon Choi (suyeon@stanford.edu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    heightBL  = sizeBL(1);
    widthBL   = sizeBL(2);
    
    % 1. determine depth of layer from given diopter data
    depths        = linspace(rangeSystem_Diopter(1), rangeSystem_Diopter(2), numLayers);
    layerSpacing  = (rangeSystem_Diopter(2) - rangeSystem_Diopter(1)) / (numLayers - 1);
    
    % 2. convert depthmap to diopter
    depthmap2Diopt = depthMap * (rangeDepthmap_Diopter(2) - rangeDepthmap_Diopter(1)) + rangeDepthmap_Diopter(1);
    
    % 3. assign value to every layer
    BL = zeros(heightBL, widthBL * numLayers);
    grad = 1 / layerSpacing;
    for j = 1 : numLayers     
        new_layer = min(1, max(0, 1 + (t-1) / 2 - grad * abs(depthmap2Diopt - depths(j))))  ;
        BL(:, (j-1) * widthBL + 1 : j * widthBL) = new_layer;       
    end
    
    % 4. exception handling (too close or too far)
    % 190403 Suyeon Choi
    tooFar   = depthmap2Diopt < (rangeSystem_Diopter(1) + layerSpacing * (t-1) / 2);
    tooClose = depthmap2Diopt > (rangeSystem_Diopter(2) - layerSpacing * (t-1) / 2);

    for j = 1:t
        BL(:, (j-1) * widthBL + 1 : j * widthBL) = BL(:, (j-1) * widthBL + 1 : j * widthBL) .* (~tooFar) + tooFar;
    end
    
    for j = numLayers - t + 1:numLayers
        BL(:, (j-1) * widthBL + 1 : j * widthBL) = BL(:, (j-1) * widthBL + 1 : j * widthBL) .* (~tooClose) + tooClose;
    end
    
    
end
