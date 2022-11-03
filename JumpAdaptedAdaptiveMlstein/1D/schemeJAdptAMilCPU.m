function outcome=schemeJAdptAMilCPU(parms, i, StartValue)
% % To generate the terminal value of the proposed Adaptive Milstein method
% %             without a reference solution
% % INPUTS:
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts
% % 'i' -- (scalar) the index of MC realisation
% % 'StartValue' -- (scalar) initial value of the process 
% % OUTPUTS:
% % 'outcome' -- (scalar) terminal value of the process  

% pull out jump sizes and times 
zetaTemp=parms.zeta(i,:);
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
    WincAMil=sqrt(h)*randn;  
    
    if h <= parms.hmin % backstop case
        % we use jump-adapted PMil with timestep 'h' for the backstop
        StartValue=schemeJAdptPMil(parms, StartValue, WincAMil, h, JumpSize);
    else
        % explicit jump-adapted Milstein
        [dft,Diff,Ddiff]=coefficients(StartValue, parms);
        StartValueMid=StartValue+dft*h+diffTerms(WincAMil,Diff,Ddiff,h);
        StartValue=StartValueMid+AMPLTD(JumpSize, StartValueMid); 
    end
    t=t+h;
end
outcome=StartValue;