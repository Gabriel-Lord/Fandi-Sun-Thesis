function outcome=schemeAMilCPU(parms, StartValue)
% % To generate the terminal value of the proposed Adaptive Milstein method
% %             without a reference solution
% % INPUTS:
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts
% % 'StartValue' -- (scalar) initial value of the process 
% % OUTPUTS:
% % 'outcome' -- (scalar) terminal value of the process  

t=0;
while abs(t-parms.T)>parms.Tol && t<=parms.T
    HcoreFunc=AMilCoreH(StartValue,parms);
    huse=max(parms.hmin, min(parms.hmax, HcoreFunc)); % choice of h
    if t+huse > parms.T  % last tep
        huse=parms.T-t;
    end
    WincAMil=sqrt(huse)*randn;  % Wiener increment
    if huse <= parms.hmin  % backstop case
        StartValue=schemePMil(parms, StartValue, WincAMil, huse);
    else  % explicit Milstein
        [dft,Diff,Ddiff]=coefficients(StartValue, parms);
        StartValue=StartValue+dft*huse+diffTerms(WincAMil,Diff,Ddiff,huse);
    end
    t=t+huse;
end
outcome=StartValue;