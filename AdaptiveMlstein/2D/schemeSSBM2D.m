function outcome=schemeSSBM2D(parms, StartValue, Winc, dtType)
%  fixed-step Split-Step Backward Milstein 
% % INPUTS:
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts
% % 'StartValue' -- (2-by-1 vector) initial vector of the process for this step
% % 'Winc' -- (2-by-1 vector)  Wiener increment from the reference solution for this step  
% % 'dtType' -- (string) type of step, either 'fixed' or 'last'
% % OUTPUTS:
% % 'outcome' -- (2-by-1 vector) terminal value of the process for this step
if strcmp('fixed',dtType)
    dt=parms.dtuse;
else  % last step
    dt=parms.dtlast;
end

input=multinewton( StartValue, dt, StartValue, parms);
Diff=DIFF(parms, input, parms.ModelNoise);
[Ddiff1, Ddiff2]=DDIFF(parms, input, parms.ModelNoise); 
outcome=input+diffTerms(Winc,Diff,Ddiff1,Ddiff2,dt,parms);


