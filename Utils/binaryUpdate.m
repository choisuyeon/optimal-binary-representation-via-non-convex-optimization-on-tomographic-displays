function [X1, Y1] = binaryUpdate(X, Y, method, opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Binary update using PALM and SART
%    for optimizing Tomographic Displays.      
%
%    Contact:
%       Suyeon Choi (suyeon@stanford.edu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    L       = opts.discreteSet;
    b       = opts.b;
    
    switch method
        % 1. using SART (previous work)
        case 'SART'
            W = opts.K  * ones(size(opts.K, 2), 1);
            V = opts.K' * ones(size(opts.K, 1), 1);            
            W(W ~= 0) = 1 ./ W(W ~= 0);
            V(V ~= 0) = 1 ./ V(V ~= 0);
            
            X1 = X + V .* (opts.K' * (W .* (b - opts.K * X)));            
            X1(X1 < 0) = 0;
            X1(X1 > 1) = 1;            
            Y1 = Y;
        
        % 2. using PALM and Three functions splitting 
        case 'PALM+PD'            
            % simple functions description
            norm_1   = @(u, C) (repmat(u, 1, length(C)) - repmat(C, length(u), 1));
            norm_2   = @(u, C) norm_1(u, C).^2;
            grad_u_H = @(z, u, C, alpha) alpha * (z.^2 .* norm_1(u, C)) * ones(length(C), 1);
            grad_z_H = @(z, u, C, alpha) alpha * z .* norm_2(u, C);
            L1       = @(z, C, alpha) alpha * max(z.^2 * ones(length(C), 1));
            L2       = @(u, C, alpha) alpha * max(max(norm_2(u, C)));

            % parameters
            alpha   = opts.alpha;

            % update
            tau = opts.gamma1 * L1(Y, L, alpha);
            X1 = prox_f(X - (1 / tau) * grad_u_H(Y, X, L, alpha), tau, opts, 'PD3O');

            sigma = opts.gamma2 * L2(X, L, alpha);
            Y1 = prox_g(Y - (1 / sigma) * grad_z_H(Y, X, L, alpha), sigma);   
        
        otherwise
            X1 = X;
            Y1 = Y;
    end
    
    
end