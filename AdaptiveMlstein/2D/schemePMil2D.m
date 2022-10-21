function outcome=schemePMil2D(parms, StartValue, Winc, dt)
% fixed-step Projected Milstein method
% % INPUTS:
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts
% % 'StartValue' -- (2-by-1 vector) initial vector of the process for this step
% % 'Winc' -- (2-by-1 vector)  Wiener increment from the reference solution for this step  
% % 'dt' -- (scalar) step size for this step
% % OUTPUTS:
% % 'outcome' -- (2-by-1 vector) terminal value of the process for this step
input= min(1 , dt^(-parms.alpha)*norm(StartValue).^(-1)).*StartValue;
[dft,Diff,Ddiff1,Ddiff2]=coefficients(input, parms);
outcome=input+dft*dt+diffTerms(Winc,Diff,Ddiff1,Ddiff2,dt,parms);
