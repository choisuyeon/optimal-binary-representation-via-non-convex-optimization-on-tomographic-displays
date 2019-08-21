classdef chambollePock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    1 step of Primal-Dual algorithms by Chambolle_Pock                   
%    for minimizing G(x) + H(Ax)          
%
%    Contact:
%       Suyeon Choi (suyeon@stanford.edu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   properties
      gamma;
      delta;
      opts;
      G     = @(x)    0;          % proximable function 
      H     = @(y)    0;          % proximable function   
      K     = @(x)    0;          % linear operator   
      AdjK  = @(y)    0;          % adjoint of K
      ProxG = @(x, t) x;         % proximal operator of G         
      ProxH = @(y, t) y;         % proximal operator of H    
   end
   methods      
       function [x1, s] = update(this)       
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           % s <- s + Ax;
           % s <- s - delta prox_{H / delta} (s / delta)
           % x <- prox_{gamma G} (x - gamma A' s)
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           x = this.opts.x;
           s = this.opts.s;
           
           iter = 1;
           for i=1:iter
               s_ = s  + this.delta .* this.K(x);
               s  = s_ - this.delta .* this.ProxH(s_ ./ this.delta, 1 ./ this.delta);
               x1 = this.ProxG(x - this.gamma .* this.AdjK(s), this.gamma);
               if this.opts.relax
                   x1 = x1 + 1 * (x1 - x);
               end
           end
       end
   end
end