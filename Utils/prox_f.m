function result = prox_f(p, tau, opt, method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Proximal operator with ||Ax - b||_2^2 within [0, 1] of 1/tau in paper
%    used PD3O by Min Yang, computes three functions at one iteration
%
%    Input 
%         p : variable of proximal operator
%         tau : coefficient of norm-2 term
%         opt : options regarding f (ex. A, b, lambda)
%
%    Contact:
%       Suyeon Choi (suyeon@stanford.edu)
%               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    k = size(opt.K, 1);    
    
    if ~isfield(opt, 'reg')
       opt.reg = 0;
    end
    
    if isfield(opt, 'lambda')
       lambda = opt.lambda;
    else
       lambda = 0;
    end
           
    if lambda > 0
        K = [opt.K; opt.D];
    else
        K = opt.K;
    end

    b = opt.b;
    switch method
        case 'PD3O'
            obj         = PrimalDual;  % using the class PrimalDual
            obj.myGradF = @(x) tau * (x - p);
            if opt.reg
               obj.myProxG = @(x, t) projSimplexWithT(min(max(x, 0), 1), opt.numLayers, opt.illumination_factor); % box -> simplex
            else
               obj.myProxG = @(x, t) min(max(x, 0), 1); % Projection on box [0, 1]
            end
            myProxH1    = @(y, t) (y + b .* t) ./ (t + 1);
            myProxH2    = @(y, t) sign(y) .* max(abs(y) - t * lambda, 0);
            obj.myProxH = @(y, t) [myProxH1(y(    1 : k  ), t(    1 : k  )); ...
                                   myProxH2(y(k + 1 : end), t(k + 1 : end))];                               
          
            obj.gamma = sum(abs(K), 1)';
            obj.gamma = 1 ./ max(1e-6, obj.gamma);
            obj.delta = sum(abs(K), 2);
            obj.delta = gather(obj.delta);
            obj.delta = 1 ./ max(1e-6, obj.delta);

            obj.input.iter       = opt.maxIter;
            obj.input.checkiter  = opt.checkIter;
            obj.input.gapDeltaThreshold = opt.epsilon;
            obj.input.x          = zeros(size(p));
            obj.input.s          = zeros(size(K, 1), 1);
            
            if opt.useGPU
                K = gpuArray(K);
            end
            obj.myA     = @(x) K  * x;
            obj.myAdjA  = @(y) K' * y;
            [x_PD3O, ~, ~]  = obj.minimize('PD3O');

            result = x_PD3O;
            if opt.useGPU
                result = gather(x_PD3O);
            end

           clearvars K;

    end

end

