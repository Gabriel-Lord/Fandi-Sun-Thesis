function [Ddiff1, Ddiff2]=DDIFF(parms,u,ModelNoise)
% function to generate derivative of diffusion 
% according to the noise that the user choose
% % INPUTS: 
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts  
% % 'u' -- (2-by-1 vector)  input value vector
% % 'ModelNoise' -- (string) type of the noise
% % OUTPUTS:
% % 'Ddiff1' 'Ddiff1'-- (2-by-2 matrix) derivative of diffusion with respect of u
if strcmp('diagonal',ModelNoise)
    Ddiff1=parms.sigma*[2*u(1)   0;  0  0];
    Ddiff2=parms.sigma*[0   0;  0  2*u(2)];
elseif strcmp('commutative',ModelNoise)
    Ddiff1=parms.sigma*[0   2*u(2);  2*u(1)  0];
    Ddiff2=parms.sigma*[0   2*u(2);  2*u(1)  0];
else  % non-commutative
    Ddiff1=parms.sigma*[2*parms.NonliDiff*u(1) 0; 0 2*u(2)];
    Ddiff2=parms.sigma*[0 1; parms.NonliDiff 0];
end