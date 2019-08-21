function [A, A_ori, check] = focalPlanes2pupilViews(pupil, system, display)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    create a projection matrix.    
%
%    Contact:
%       Suyeon Choi (suyeon@stanford.edu)
%       original version by Seungjae Lee (seungjae1012@gmail.com)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ticProjectionMatrix = tic;
err = 1e-6;
diopt = linspace(system.DOF(1) + err, system.DOF(2), system.numLayers); 

%% focal plane
focalPlaneDepth = 1./ diopt * 1000 + pupil.depth;
focalPlaneSize = zeros(system.numLayers, 2);
pixelSizes = zeros(system.numLayers, 2);
for layerIdx = 1:system.numLayers
    focalPlaneSize(layerIdx, :)    = 2 * tan(system.FOV / 2) * (focalPlaneDepth(layerIdx) - pupil.depth) * [1 1]; % Dimension : mm
    pixelSizes(layerIdx, :) = focalPlaneSize(layerIdx, :) ./ display.resolution;
end

%% pupil plane     
positionX = linspace(-pupil.size(1) / 2, pupil.size(1) / 2, pupil.numViewpoints(1));
positionY = linspace(-pupil.size(2) / 2, pupil.size(2) / 2, pupil.numViewpoints(2));

[PY, PX] = meshgrid(positionY, positionX);
pupilVecRes = prod(pupil.numViewpoints);

pupilPlaneDepth        = 1 / diopt(1) * 1000 + pupil.depth;
pupilPlaneSize         = 2 * tan(system.FOV/2) * (pupilPlaneDepth - pupil.depth) * [1 1];
pupilPlanePixelSize = pupilPlaneSize ./ pupil.resolution;
pupilPlaneX          = linspace(-pupilPlaneSize(1) / 2 + pupilPlanePixelSize(1) / 2, ...
                           pupilPlaneSize(1) / 2 - pupilPlanePixelSize(1) / 2, ...
                           pupil.resolution(1));
pupilPlaneY          = linspace(-pupilPlaneSize(2) / 2 + pupilPlanePixelSize(2) / 2, ...
                           pupilPlaneSize(2) / 2 - pupilPlanePixelSize(2) / 2, ...
                           pupil.resolution(2));

%% Create Projection Matrix
projectionRows = zeros(1, prod(pupil.numViewpoints) * prod(pupil.resolution));
projectionCols = zeros(1, prod(display.resolution));

idxCount = 1;
for layerIdx = 1:system.numLayers
    for positionIdx =  1:pupilVecRes
        for pupilIdx = 1:pupil.resolution(1)
            tanY = (pupilPlaneY) / (pupilPlaneDepth - pupil.depth);
            tanX = (pupilPlaneX(pupilIdx)) / (pupilPlaneDepth - pupil.depth);

            focalPlaneY = (focalPlaneDepth(layerIdx) - pupil.depth) * tanY + PY(positionIdx);
            focalPlaneX = (focalPlaneDepth(layerIdx) - pupil.depth) * tanX + PX(positionIdx);
            
            focalPlaneYIdx = ceil((focalPlaneY + focalPlaneSize(layerIdx, 2) / 2) / pixelSizes(layerIdx, 2));
            focalPlaneXIdx = ceil((focalPlaneX + focalPlaneSize(layerIdx, 1) / 2) / pixelSizes(layerIdx, 1));
            
            filter = (focalPlaneYIdx > 0) .* (focalPlaneYIdx < display.resolution(2) + 1);
            filterIdxs = find(filter);
            focalPlaneYIdx = focalPlaneYIdx(filterIdxs);
            
            if((focalPlaneXIdx>0) && (focalPlaneXIdx < display.resolution(1) + 1))
                projectionRows(idxCount:idxCount + size(filterIdxs, 2) - 1) = ...
                        pupilIdx + pupil.resolution(1) * (filterIdxs - 1) + prod(pupil.resolution) * (positionIdx - 1);
                projectionCols(idxCount:idxCount + size(filterIdxs, 2) - 1) = ...
                        focalPlaneXIdx + display.resolution(1) * (focalPlaneYIdx - 1) + prod(display.resolution) * (layerIdx - 1);
                idxCount = idxCount + size(filterIdxs, 2);
            end
        end
    end
end
count = idxCount - 1;
projectionRowsFiltered = projectionRows(1:count);
projectionColsFiltered = projectionCols(1:count);

A = sparse(projectionRowsFiltered, ...
           projectionColsFiltered, ...
           1, ...
           prod(pupil.resolution) * prod(pupil.numViewpoints), ...
           prod(display.resolution) * system.numLayers); % 2 layer

A_ori = A;
is = projectionRowsFiltered;
js = projectionColsFiltered;
isColorAugmented = [is(1:count) is(1:count) + size(A_ori, 1) is(1:count) + 2 * size(A_ori, 1)];
jsColorAugmented = [js(1:count) js(1:count) + size(A_ori, 2) js(1:count) + 2 * size(A_ori, 2)];

A = sparse(isColorAugmented, ...
           jsColorAugmented, ...
           1, ...
           3 * prod(pupil.resolution) * prod(pupil.numViewpoints), ...
           3 * prod(display.resolution) * system.numLayers);
       
ticProjectionMatrix = toc(ticProjectionMatrix);
disp(['   - Projection Matrix created : ', num2str(ticProjectionMatrix), 's']);

%% For check (TODO)
check.numLayers = system.numLayers;
check.disp_resolution = display.resolution;
check.pupil_resolution = pupil.resolution;
check.content = system.content;

%% demalloc
clearvars I J I_color J_color 
clearvars projectionCols projectionColsFiltered;
clearvars projectionRows projectionRowsFiltered;
end