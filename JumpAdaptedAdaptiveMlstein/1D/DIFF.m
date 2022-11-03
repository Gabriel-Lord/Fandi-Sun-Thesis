function outcome=DIFF(sigma,u,ModelNoise)
% function to generate diffusion 
% according to the noise that the user choose
% % INPUTS: 
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts  
% % 'u' -- (scalar)  input value 
% % 'ModelNoise' -- (string) type of the noise
% % OUTPUTS:
% % 'outcome' -- (scalar) of diffusion
if strcmp('additive',ModelNoise)
    outcome=sigma;
elseif strcmp('mul1',ModelNoise)
    outcome=sigma*u;
else
    outcome=sigma*(1-u^2);
end