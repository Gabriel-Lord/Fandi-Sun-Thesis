function outcome=schemeTMil(parms, StartValue, Winc, dtType)
%  fixed-step Tamed Milstein method
% % INPUTS:
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts
% % 'StartValue' -- (scalar) initial value of the process for this step
% % 'Winc' -- (scalar)  Wiener increment from the reference solution for this step  
% % 'dtType' -- (string) type of step, either 'fixed' or 'last'
% % OUTPUTS:
% % 'outcome' -- (scalar) terminal value of the process for this step

if strcmp('fixed',dtType)
    dt=parms.dtuse;
else  % last step
    dt=parms.dtlast;
end

scale=1+dt*abs(StartValue)^4;
[dft,Diff,Ddiff]=coefficients(StartValue, parms);
outcome=StartValue+(dft*dt+diffTerms(Winc,Diff,Ddiff,dt))/scale;
