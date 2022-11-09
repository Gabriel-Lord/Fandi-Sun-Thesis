function outcome=schemeJAdptSSBM2D(parms, StartValue, Winc, dt, JumpSize)
%  fixed-step jump-adapted Split-Step Backward Milstein 
% % INPUTS:
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts
% % 'StartValue' -- (2-by-1 vector) initial vector of the process for this step
% % 'Winc' -- (2-by-1 vector)  Wiener increment from the reference solution for this step  
% % 'dt' -- (scalar) step size for this step
% % 'JumpSize' -- (scalar)  
% % OUTPUTS:
% % 'outcome' -- (2-by-1 vector) terminal value of the process for this step
 
input=multinewton( StartValue, dt, StartValue, parms);
Diff=DIFF(parms, input, parms.ModelNoise);
[Ddiff1, Ddiff2]=DDIFF(parms, input, parms.ModelNoise); 
outcomeMid=input+diffTerms(Winc,Diff,Ddiff1,Ddiff2,dt,parms);
outcome=outcomeMid+AMPLTD(JumpSize,outcomeMid);


