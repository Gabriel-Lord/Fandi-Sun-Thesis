function outcome=schemeJAdptAMilCPU2DLevy(parms, i, StartValue)
% % To generate the terminal value of the proposed Adaptive Milstein method
% %             without a reference solution
% % INPUTS:
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts
% % 'i' -- (scalar) index of the MC realisation
% % 'StartValue' -- (2-by-1 vector) initial vector of the process 
% % OUTPUTS:
% % 'outcome' -- (2-by-1 vector) terminal value of the process  

t=0;
% pull out jump sizes, jump times, reference times
% and reference time increments
zetaTemp=parms.zeta(parms.rank*i-(parms.rank-1):parms.rank*i,:);
JumpTimesTemp=[nonzeros(parms.JumpTimes(i,:))', 2*parms.T];

% Y2 and Y4 are needed for Levy area computation
Y2=0; Y4=0;
while abs(t-parms.T)>parms.Tol && t<=parms.T
    HcoreFunc=AMilCoreH(StartValue,parms);
    h=max(parms.hmin, min(parms.hmax, HcoreFunc)); % choice of h
    % check if there is a jump in this step and look for its jump size
    [h, ~, JumpSize, JumpTimesTemp, zetaTemp]=...
        JumpSizeFinder(t, h, JumpTimesTemp, zetaTemp);
    
    if t+h > parms.T  % last tep 
        h=parms.T-t;
        parms.LA=0;
        WincAMil=randn(parms.rank,1).*sqrt(h); 
    else % all other steps, could have jump or not
        LevyStepNo=round(1/h);   
        LevyStepSize=h/LevyStepNo;   
        WincLevy=randn(parms.rank, LevyStepNo).*sqrt(LevyStepSize); 
        % approximating Levy area
        [parms.LA,Y2,Y4]=funcLevy(LevyStepNo,Y2,Y4,WincLevy);
        WincAMil=sum(WincLevy,2);
    end
   if h <= parms.hmin % backstop case
        % we use PMil with timestep 'huse' for the backstop
        StartValue=schemeJAdptPMil2D(parms, StartValue, WincAMil, h, JumpSize);
    else 
        % explicit Milstein
        [dft,Diff,Ddiff1,Ddiff2]=coefficients(StartValue, parms);
        StartValueMid=StartValue+dft*h+diffTerms(WincAMil,Diff,Ddiff1,Ddiff2,h,parms);
        StartValue=StartValueMid+AMPLTD(JumpSize, StartValueMid); 
    end
    t=t+h;
end
outcome=StartValue;