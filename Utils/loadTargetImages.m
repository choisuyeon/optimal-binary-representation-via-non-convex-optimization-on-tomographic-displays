function b = loadTargetImages(pupil, targetDir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    load target images
%
%    Contact:
%       Suyeon Choi (suyeon@stanford.edu)
%       original version by Seungjae Lee (seungjae1012@gmail.com)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ticLoadTarget = tic;
b = ones(pupil.resolution(1), pupil.resolution(2), prod(pupil.numViewpoints), 3);
ix = 1:pupil.numViewpoints(1); iy = 1:pupil.numViewpoints(2);
targetPR = pupil.numViewpoints(1);
for idxX = 1:pupil.numViewpoints(1)
    for idxY = 1:pupil.numViewpoints(2)
        i = iy(idxY) + (ix(idxX)-1) * targetPR;
        filename = sprintf('%s/%04d.png', ...
            targetDir, ...
            mod(i-1, targetPR) * targetPR + floor((i-1) / targetPR) + 1);
        
        for c = 1:3
            temp = double(imresize(imread(filename), pupil.resolution, 'nearest')) ./ 255;
            b(:, :, idxY + (idxX-1) * pupil.numViewpoints(1), c) = temp(:, :, c);
        end
        
    end
end
b = b(:);

ticLoadTarget = toc(ticLoadTarget);
disp(['   - Target loaded : ', num2str(ticLoadTarget), 's']);
end
