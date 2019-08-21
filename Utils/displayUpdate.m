function [X1, Y1] = displayUpdate(X, Y, method, opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Display update using PD algorithm by Chambolle and Pock            
%    for optimizing Tomographic Displays.      
%
%    Contact:
%       Suyeon Choi (suyeon@stanford.edu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    b = opt.b;
    
    switch method
        % 1. using SART (previous work)
        case 'SART'
            W         = opt.K  * ones(size(opt.K, 2), 1);
            V         = opt.K' * ones(size(opt.K, 1), 1);            
            W(W ~= 0) = 1 ./ W(W ~= 0);
            V(V ~= 0) = 1 ./ V(V ~= 0);
            
            X1         = X + V .* (opt.K' * (W .* (b - opt.K * X)));            
            X1(X1 < 0) = 0;
            X1(X1 > 1) = 1;            
            Y1 = Y;
            
        case 'PD'
            obj       = chambollePock;
            obj.G     = @(x)    sum((x > 1) + (x < 0)) * 1e12;
            obj.H     = @(y)    0.5 * (y - b)' * (y - b);        
            obj.K     = @(x)    opt.K  * x;
            obj.AdjK  = @(y)    opt.K' * y;
            obj.ProxG = @(x, t) min(max(x, 0), 1); % Projection on box [0, 1]
            obj.ProxH = @(y, t) (y + b .* t) ./ (t + 1);

            obj.gamma = sum(abs(opt.K), 1)';
            obj.delta = sum(abs(opt.K), 2);
            obj.gamma = 1 ./ max(1e-3, obj.gamma);
            obj.delta = 1 ./ max(1e-3, obj.delta);

            obj.opts.x     = X;
            obj.opts.s     = zeros(size(b));
            obj.opts.relax = 0;

            [X1, Y1]  = obj.update();      
        otherwise
            X1 = X;
            Y1 = Y;
    end
end