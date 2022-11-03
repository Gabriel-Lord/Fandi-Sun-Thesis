function outcome=schemeJAdptTMil(parms, StartValue, Winc, dt, JumpSize)
%  fixed-step jump-adapted Tamed Milstein method
% % INPUTS:
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts
% % 'StartValue' -- (scalar) initial value of the process for this step
% % 'Winc' -- (scalar)  Wiener increment from the reference solution for this step  
% % 'dt' -- (scalar) step size
% % 'JumpSize' -- (scalar)  jump size occurred at the end of this step
% % OUTPUTS:
% % 'outcome' -- (scalar) terminal value of the process for this step
 

scale=1+dt*abs(StartValue)^4;
[dft,Diff,Ddiff]=coefficients(StartValue, parms);
outcomeMid=StartValue+(dft*dt+diffTerms(Winc,Diff,Ddiff,dt))/scale;
outcome=outcomeMid+AMPLTD(JumpSize,outcomeMid);
