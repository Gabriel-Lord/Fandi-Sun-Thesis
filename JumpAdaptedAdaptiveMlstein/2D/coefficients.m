function [dft,Diff,Ddiff1,Ddiff2]=coefficients(u,parms)
% function to generate 
% drift, diffusion and the derivative of diffusion
% % INPUTS: 
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts  
% % 'u' -- (2-by-1 vector)  input value vector 
% % OUTPUTS:
% % 'dft' -- (2-by-1 vector)  of drift 
% % 'Diff' -- (2-by-2 matrix) of diffusion
% % 'Ddiff1' 'Ddiff2' -- (2-by-2 matrix) derivatives of diffusion
dft=DFT(parms.nonlinearity, u);
Diff=DIFF(parms, u, parms.ModelNoise);
[Ddiff1, Ddiff2]=DDIFF(parms, u, parms.ModelNoise);
