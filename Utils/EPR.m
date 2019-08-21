function [rB, energyCompensated] = EPR(B, L, layerResolution, numLayers)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    simple energy-preserving rounding implementation
%    (see section 4.4 in https://doi.org/10.1364/OE.27.024362)
%
%    Contact:
%       Suyeon Choi (suyeon@stanford.edu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ticEPR = tic;
    fprintf('   - Energy Preserving rounding started... \n');
    
    d        = L(2) - L(1);
    B        = reshape(B, prod(layerResolution), numLayers);
    rB       = interp1(L, L, B, 'nearest', 'extrap'); 
    delta_rB = rB - B;
    sumB     = sum(B, 2);
    rsumB    = round(sumB);  %% sum or rB should be this.
                          
    sumrB    = sum(rB, 2);
    
    delta_rsumB = rsumB - sumrB;  %% this is disparity of energy per pixel.
                                  % we will compensate this.
        
    wrongCount = 0;
    for k = 1:prod(layerResolution)
        energyCompensate = delta_rsumB(k) / d;
        deltaE = delta_rB(k, :);
        if energyCompensate > 0
            [~, idx_least_order] = sort(deltaE, 'ascend');
            for i = 1:energyCompensate
                % compensate those depths
                originalValue = rB(k, idx_least_order(i));
                if originalValue >= L(end)
                    wrongCount = wrongCount + 1;
                else
                    rB(k, idx_least_order(i)) = rB(k, idx_least_order(i)) + d;
                end
            end            
        elseif energyCompensate < 0
            [~, idx_greatest_order] = sort(deltaE, 'descend');
            for i = 1:-energyCompensate
                % compensate those depths
                originalValue = rB(k, idx_greatest_order(i));
                if originalValue <= L(1)
                    wrongCount = wrongCount + 1;
                else
                    rB(k, idx_greatest_order(i)) = originalValue - d;
                end
            end       
        end
    end
    
    after_delta_rsumB = rsumB - sum(rB, 2);
    % Verification
    if sum(after_delta_rsumB) > 0 || wrongCount > 0
        fprintf('results is not perfect! wrong_count : %d \n', wrongCount);
    end
        
    rB = rB(:);
    rB = min(max(0, rB), 1);
    energyCompensated = reshape(delta_rsumB, layerResolution);

elapsedTime = toc(ticEPR);
    fprintf('   - energy preserving rounding finished. elapsed time:%.1f \n', elapsedTime);
    
end