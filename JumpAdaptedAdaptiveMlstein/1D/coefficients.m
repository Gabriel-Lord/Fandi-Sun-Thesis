function [dft,Diff,Ddiff]=coefficients(u,parms)
% function to generate 
% drift, diffusion and the derivative of diffusion
% % INPUTS: 
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts  
% % 'u' -- (scalar)  input value  
% % OUTPUTS:
% % 'dft' -- (scalar)  of drift 
% % 'Diff' -- (scalar) of diffusion
% % 'Ddiff'  -- (scalar) derivatives of diffusion
dft=DFT(parms.nonlinearity, u);
Diff=DIFF(parms.sigma, u, parms.ModelNoise);
Ddiff=DDIFF(parms.sigma, u, parms.ModelNoise);
