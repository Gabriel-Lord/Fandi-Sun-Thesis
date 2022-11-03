function outcome=schemeJAdptSSBM(parms, StartValue, Winc, dt, funcCubicDrift, JumpSize)
%  fixed-step jump-adapted Split-Step Backward Milstein
% % INPUTS:
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts
% % 'StartValue' -- (scalar) initial value of the process for this step
% % 'Winc' -- (scalar)  Wiener increment from the reference solution for this step  
% % 'dt' -- (scalar) step size
% % 'funcCubicDrift' -- function name needed for solving drift
% % 'JumpSize' -- (scalar)  jump size occurred at the end of this step
% % OUTPUTS:
% % 'outcome' -- (scalar) terminal value of the process for this step
 
% following is to solve the cubic equation for drift
p=(dt-1)/(dt*parms.nonlinearity);
q= StartValue/(dt*parms.nonlinearity);
s = sqrt(q^2/4 - p^3/27);
input =  funcCubicDrift(q/2 + s)+funcCubicDrift(q/2-s);
% substituting 'input' into the diffusion terms
Diff=DIFF(parms.sigma, input, parms.ModelNoise);
Ddiff=DDIFF(parms.sigma, input, parms.ModelNoise);
outcomeMid=input+diffTerms(Winc,Diff,Ddiff,dt);
outcome=outcomeMid+AMPLTD(JumpSize,outcomeMid);

