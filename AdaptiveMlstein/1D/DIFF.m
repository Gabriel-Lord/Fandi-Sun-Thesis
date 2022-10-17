function outcome=DIFF(sigma,u,ModelNoise)
% function to generate diffusion 
% according to the noise that the user choose
if strcmp('additive',ModelNoise)
    outcome=sigma;
elseif strcmp('mul1',ModelNoise)
    outcome=sigma*u;
else
    outcome=sigma*(1-u^2);
end