function outcome=schemeAMilCPU2D(parms, StartValue)
% % To generate the terminal value of the proposed Adaptive Milstein method
% %             without a reference solution
% % INPUTS:
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts
% % 'StartValue' -- (2-by-1 vector) initial vector of the process 
% % OUTPUTS:
% % 'outcome' -- (2-by-1 vector) terminal value of the process  

t=0;
while abs(t-parms.T)>parms.Tol && t<=parms.T
    HcoreFunc=AMilCoreH(StartValue,parms);
    huse=max(parms.hmin, min(parms.hmax, HcoreFunc)); % choice of h
    if t+huse > parms.T  % last tep
        huse=parms.T-t;
    end
    WincAMil=sqrt(huse)*randn(parms.rank,1);  % Wiener increment
    if huse <= parms.hmin  % backstop case
        StartValue=schemePMil2D(parms, StartValue, WincAMil, huse);
    else  % explicit Milstein
        [dft,Diff,Ddiff1,Ddiff2]=coefficients(StartValue, parms);
        StartValue=StartValue+dft*huse+diffTerms(WincAMil,Diff,Ddiff1,Ddiff2,huse,parms);
    end
    t=t+huse;
end
outcome=StartValue;