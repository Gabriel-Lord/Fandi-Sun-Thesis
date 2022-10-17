function outcome=diffTerms(W,diff,Ddiff,dt)
% to generate the terms from diffusion by Milstein method
outcome=diff*W+0.5*Ddiff*diff*(W.^2-dt);