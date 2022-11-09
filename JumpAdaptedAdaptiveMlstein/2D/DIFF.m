function outcome=DIFF(parms,u,ModelNoise)
% function to generate diffusion 
% according to the noise that the user choose
% % INPUTS: 
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts  
% % 'u' -- (2-by-1 vector)  input value vector
% % 'ModelNoise' -- (string) type of the noise
% % OUTPUTS:
% % 'outcome' -- (2-by-2 matrix) of diffusion
if strcmp('diagonal',ModelNoise)
    outcome=parms.sigma*[u(1)^2 0; 0 u(2)^2];
elseif strcmp('commutative',ModelNoise)
    outcome=parms.sigma*[u(2)^2 u(2)^2; u(1)^2 u(1)^2];
else  % non-commutative
    outcome=parms.sigma*...
        [parms.NonliDiff*u(1)^2 u(2); u(2)^2 parms.NonliDiff*u(1)];
end