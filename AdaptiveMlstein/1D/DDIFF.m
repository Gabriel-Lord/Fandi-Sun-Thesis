function outcome=DDIFF(sigma,u,ModelNoise)
% function to generate derivative of diffusion 
% according to the noise that the user choose
% % INPUTS: 
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts  
% % 'u' -- (scaler)  input value  
% % 'ModelNoise' -- (string) type of the noise
% % OUTPUTS:
% % 'Ddiff' -- (scaler) derivative of diffusion with respect of u
if strcmp('additive',ModelNoise)
    outcome=0;
elseif strcmp('mul1',ModelNoise)
    outcome=sigma;
else
    outcome=-2*sigma*u;
end