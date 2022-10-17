function [dft,Diff,Ddiff]=coefficients(u,parms)
% function to generate 
% drift, diffusion and the derivative of diffusion
dft=DFT(parms.nonlinearity, u);
Diff=DIFF(parms.sigma, u, parms.ModelNoise);
Ddiff=DDIFF(parms.sigma, u, parms.ModelNoise);
