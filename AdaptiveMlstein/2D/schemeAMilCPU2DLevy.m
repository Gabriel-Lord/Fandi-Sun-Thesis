function outcome=schemeAMilCPU2DLevy(parms, StartValue)
% % To generate the terminal value of the proposed Adaptive Milstein method
% %             without a reference solution
% % INPUTS:
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts
% % 'StartValue' -- (2-by-1 vector) initial vector of the process 
% % OUTPUTS:
% % 'outcome' -- (2-by-1 vector) terminal value of the process  

t=0;
% Y2 and Y4 are needed for Levy area computation
Y2=0; Y4=0;
while abs(t-parms.T)>parms.Tol && t<=parms.T
    HcoreFunc=AMilCoreH(StartValue,parms);
    h=max(parms.hmin, min(parms.hmax, HcoreFunc)); % choice of h
    % Levy area setting
    delta=h^2; % size of ministep
    % number of ministep within one h
    K=floor(h/delta);
    % resize h to make it a multiple of ministeps 
    huse=K*delta;
    % mini Wiener increments 
    % for mini step within one h
    dWla=randn(parms.rank,K)*sqrt(delta);
    % approximating Levy area
    [parms.LA,Y2,Y4]=funcLevy(K,Y2,Y4,dWla);
    % Wiener increment to use 
    WincAMil=sum(dWla,2);

    if t+huse > parms.T  % last tep
        huse=parms.T-t;
        % we remove Levy area approximation in the last step
        parms.LA=0;
    end 
   if huse <= parms.hmin  % backstop case
        StartValue=schemePMil2D(parms, StartValue, WincAMil, huse);
    else  % explicit Milstein
        [dft,Diff,Ddiff1,Ddiff2]=coefficients(StartValue, parms);
        StartValue=StartValue+dft*huse+diffTerms(WincAMil,Diff,Ddiff1,Ddiff2,huse,parms);
    end
    t=t+huse;
end
outcome=StartValue;