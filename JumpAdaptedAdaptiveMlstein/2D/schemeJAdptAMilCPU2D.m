function outcome=schemeJAdptAMilCPU2D(parms, i, StartValue)
% % To generate the terminal value of the proposed
% %           jump-adapted Adaptive Milstein method
% %            without a reference solution
% % INPUTS:
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts
% % 'StartValue' -- (2-by-1 vector) initial vector of the process 
% % 'i' -- (scalar) the index of MC realisation
% % OUTPUTS:
% % 'outcome' -- (2-by-1 vector) terminal value of the process  

% pull out jump sizes and times 
zetaTemp=parms.zeta(parms.rank*i-(parms.rank-1):parms.rank*i,:);
JumpTimesTemp=[nonzeros(parms.JumpTimes(i,:))', parms.T+parms.Tol]; 

t=0;
while abs(t-parms.T)>parms.Tol && t<=parms.T
    % core function of choosing 'h'
    HcoreFunc=AMilCoreH(StartValue,parms);
    % bounding the function value by hmin and hmax
    h=max(parms.hmin, min(parms.hmax, HcoreFunc)); % choice of h
    
    if t+h > parms.T  % last tep
        h=parms.T-t;
    end 
    
    % check if there is a jump in this step and look for its jump size 
    [h, ~, JumpSize, JumpTimesTemp, zetaTemp]=...
                        JumpSizeFinder(t, h, JumpTimesTemp, zetaTemp);
                    
    % Wiener increment
    WincAMil=sqrt(h)*randn(parms.rank,1);  
    if h <= parms.hmin  % backstop case
        StartValue=schemeJAdptPMil2D(parms, StartValue, WincAMil, h, JumpSize);
    else  % explicit Milstein
        [dft,Diff,Ddiff1,Ddiff2]=coefficients(StartValue, parms);
        StartValueMid=StartValue+dft*h+diffTerms(WincAMil,Diff,Ddiff1,Ddiff2,h,parms);
        StartValue=StartValueMid+AMPLTD(JumpSize, StartValueMid);
    end
    t=t+h;
end
outcome=StartValue;