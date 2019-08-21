classdef PrimalDual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Primal-Dual algorithms for minimizing                            
%    F(x) + G(x) + H(Ax)                                              
%    F is differentiable, and both G and H are proximable             
%
%    Three methods are implemented: Condat-Vu, PDFP, and PD3O        
%    Reference: 
%       M. Yan, "A Primal-Dual Three Operator Splitting Scheme", arxiv
%       1611.????, Nov., 2016.
%
%    Contact:
%       Ming Yan yanm @ math.msu.edu
%       Downloadable from https://github.com/mingyan08/PD3O                                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   properties
      A;
      gamma;        % parameter gamma
      delta;        % parameter delta
      input;
      myF       =  @(x)   0;            % smooth function    F: x --> F(x)
      myG       =  @(x)   0;            % proximable fcn     G: x --> G(x)
      myH       =  @(y)   0;            % proximable fcn     H: y --> H(y)
      myHA      =  @(y)   0;            % conjugate fcn of H
      myA       =  @(x)   0;            % linear operator    A: x --> A(x)
      myGradF   =  @(x)   0;            % gradient of F      gradF: x    --> grad(F)(x)
      myProxG   =  @(x,t) x;            % prox of G          proxG: x,t  --> prox(t.G)(x)
      myProxH   =  @(y,t) y;            % prox of H          proxH: y,t  --> prox(t.H)(y)
      myAdjA    =  @(y)   0;            % adjoint of A       adjA:  y    --> A'*(y)
   end
   methods
       function p = E(this, x)  
           % return the function value for given x: F(x) + G(x) + H(Ax) 
           p = this.myF(x) + this.myG(x) + this.myH(this.myA(x));
       end       
      
       function [x, s, E] = minimize(this, method)       
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % s^+ = (I - Prox_{gamma/lambda H})(s + A * x_bar)
           % x^+ = Prox_{gamma G}(x - gamma GradF - lambda A' * s^+)
           %     = Prox_{gramm G}(xh - lambda A' * s^+)
           % x_bar^+ = 2 x^+ - x + gamma GradF - gamma GradF'
           %         = 2 x^+ - xh - gamma GradF'
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           if ~isfield(this.input, 'x')
                error('Error. Put the initial primal x in this.input.x')               
           end
           x     = this.input.x;
           if ~isfield(this.input, 's')
                error('Error. Put the initial dual s in this.input.s')               
           end
           s     = this.input.s;
           if isfield(this.input, 'iter')
               iter  = this.input.iter;
           else
               iter = 1000;
           end
           
           if isfield(this.input, 'checkiter')
               checkiter  = this.input.checkiter;
           else
               checkiter = 50;
           end
           
           if isfield(this.input, 'gapDeltaThreshold')
               gapDeltaThreshold  = this.input.gapDeltaThreshold;
           else
               gapDeltaThreshold = 1e-5;
           end
           
           x_bar = x;
           E     = zeros(iter, 1);
           gradF = this.myGradF(x);
           
           processCmd = [];
           PDGap0 = 1000;
           for i = 1:iter
               x0       = x;
               s0       = s;
               xh       = x - this.gamma .* gradF;  % x - \gamma \nabla F(x)
               sh       = s + this.delta .* this.myA(x_bar);     % s + A * x_bar time cost step
               s        = sh - this.delta .* this.myProxH(sh ./ this.delta, 1 ./ this.delta);
               x        = this.myProxG(xh - this.gamma .* this.myAdjA(s), this.gamma);
               
               gradF    = this.myGradF(x);
               
               switch method
                   case 'PD3O'  % Primal-Dual Three Operator
                       x_bar    = 2*x - xh - this.gamma .* gradF;                   
               end
               
               if i >= 0
                   if mod(i, checkiter) == 0 
                       PDGap =   norm((x - x0) ./ this.gamma - this.myAdjA(s - s0)) / numel(x) ...
                               + norm((s - s0) ./ this.delta - this.myA(x - x0)) / numel(s);
                       processCmd = print2cmd(processCmd, sprintf(['Calculating prox_f... Iter: #%d : PD Gap %.', num2str(ceil(-log10(gapDeltaThreshold))), 'f', '\n'], i, PDGap));
                       
                       if PDGap < gapDeltaThreshold
                           break;
                       end
                       if abs(PDGap0 - PDGap) < gapDeltaThreshold / 10000 % abs val
                           disp('abs');
                           break;
                       end
                       if abs(PDGap0 - PDGap) / PDGap0 < 1e-3 % 2nd diff
                           disp('2nd diff');
                           break;
                       end
                       PDGap0 = PDGap;
                   end
               end
           end
       end
   end
end