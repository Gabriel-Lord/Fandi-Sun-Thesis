function outcome=schemePMil(parms, StartValue, Winc, dt)
% fixed-step Projected Milstein method
% % INPUTS:
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts
% % 'StartValue' -- (scalar) initial value of the process for this step
% % 'Winc' -- (scalar)  Wiener increment from the reference solution for this step  
% % 'dt' -- (scalar) step size
% % OUTPUTS:
% % 'outcome' -- (scalar) terminal value of the process for this step
input= min(1 , dt^(-parms.alpha)*abs(StartValue)^(-1)) *StartValue;
[dft,Diff,Ddiff]=coefficients(input, parms);
outcome=input+dft*dt+diffTerms(Winc,Diff,Ddiff,dt);
