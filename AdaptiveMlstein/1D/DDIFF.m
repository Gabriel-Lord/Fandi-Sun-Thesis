function outcome=DDIFF(sigma,u,ModelNoise)
% function to generate derivative of diffusion 
% according to the noise that the user choose
if strcmp('additive',ModelNoise)
    outcome=0;
elseif strcmp('mul1',ModelNoise)
    outcome=sigma;
else
    outcome=-2*sigma*u;
end