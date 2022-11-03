function outcome=schemeJAdptPMil(parms, StartValue, Winc, dt, JumpSize)
% fixed-step jump-adapted Projected Milstein method
% % INPUTS:
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts
% % 'StartValue' -- (scalar) initial value of the process for this step
% % 'Winc' -- (scalar)  Wiener increment from the reference solution for this step  
% % 'dt' -- (scalar) step size
% % 'JumpSize' -- (scalar)  jump size occurred at the end of this step
% % OUTPUTS:
% % 'outcome' -- (scalar) terminal value of the process for this step
input= min(1 , dt^(-parms.alpha)*abs(StartValue)^(-1)) *StartValue;
[dft,Diff,Ddiff]=coefficients(input, parms);
outcomeMid=input+dft*dt+diffTerms(Winc,Diff,Ddiff,dt);
outcome=outcomeMid+AMPLTD(JumpSize,outcomeMid);
