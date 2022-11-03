function outcome=diffTerms(W,diff,Ddiff,dt)
% to generate the terms from diffusion by Milstein method
% % INPUTS:
% % 'W' -- (scalar)  Wiener increment for this step  
% % 'Diff' -- (scalar)  diffusion  
% % 'Ddiff' -- (scalar)  derivative of diffusion  
% %           with respect to input  
% % 'dt' -- (scalar)  step size  
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts  
% % OUTPUTS:
% % 'outcome' -- (scalar) contains all terms from diffusion
outcome=diff*W+0.5*Ddiff*diff*(W.^2-dt);