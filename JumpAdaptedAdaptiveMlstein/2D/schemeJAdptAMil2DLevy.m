function [outcome, j, hfix]=schemeJAdptAMil2DLevy(parms, StartValue, Winc1MC, i, j, hfix)
% % To generate the terminal value of the proposed Adaptive Milstein method
% %           based on a reference solution
% % INPUTS:
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts
% % 'StartValue' -- (2-by-1 vector) initial vector of the process
% % 'Winc1MC' -- (2-by-n vector) reference Levy Wiener increments
% % 'j' -- (scalar) number of adaptive steps taken on all MC realisations
% % 'hfix' -- (scalar) sum length of adaptive steps taken on all MC realisations
% % OUTPUTS:
% % 'outcome' -- (2-by-1 vector) terminal value of the process
% % 'j' -- (scalar) updated value
% % 'hfix' -- (scalar) updated value

% Since each step 'h' is adaptive and recalculated after each mesh grid
%           we do not know how many reference Wiener increaments
%           they will take for each step
% Therefore, 'Marker' below acts like a bookmark
%           which will record the starting position of each adaptive step
%           on the corresponding reference line
% 'Marker' is updated after each adaptive step is done
Marker=0;

% we add each adaptive step to 't'
% until it reaches terminal time 'parms.T'
t=0;

% pull out reference time points with Levy steps
times1MC=parms.times1MC;
timesIncLevy=parms.timesIncLevy;
% pull out jump sizes, jump times, reference times
% and reference time increments
zetaTemp=parms.zeta(parms.rank*i-(parms.rank-1):parms.rank*i,:);
JumpTimesTemp=[nonzeros(parms.JumpTimes(i,:))', 2*parms.T];

% Y2 and Y4 are needed for Levy area computation
Y2=0; Y4=0;
while abs(t-parms.T)>parms.Tol && t<=parms.T
    
    % core function of choosing 'h'
    HcoreFunc=AMilCoreH(StartValue,parms);
    % bounding the function value by hmin and hmax
    h=max(parms.hmin, min(parms.hmax, HcoreFunc));
    
    % tau=T+Tol < t+dt > T ?
    % check if there is a jump in this step and look for its jump size
    [h, ~, JumpSize, JumpTimesTemp, zetaTemp]=...
        JumpSizeFinder(t, h, JumpTimesTemp, zetaTemp);
    
    % separate the last step (no Levy) and the rest
    if t+h > parms.T  % last tep with no jump
        h=parms.T-t;
        WincAMil=sum(Winc1MC(:, Marker+1 : end),2);
        % we don't apply Levy approximation in the last step
        % in case the Levy steps are too small
        parms.LA=0;
    else % all other steps, could have jump or not
        % number of Levy steps
        LevyStepNo=round(1/h);
        % size of each Levy step
        LevyStepSize=h/LevyStepNo;
        % create empty vector for storing Levy Wiener increments
        WincLevy=zeros(parms.rank,LevyStepNo);
        for k=1:LevyStepNo
            times1MCDiff=t+k*LevyStepSize-times1MC;
            % index of the right mesh point of ''t+huse''
            hLocCeil=find(times1MCDiff<=0,1);
            % Wiener increment up to the left mesh point of ''t+h''
            WincLevy(:,k)=sum(Winc1MC(:, Marker+1 : hLocCeil-2),2);
            % left portion of this step within two reference mesh points
            bbDt1= t+k*LevyStepSize - times1MC(hLocCeil-1);
            % whole Wiener increment of these two reference mesh points
            bbDW=Winc1MC(:,hLocCeil-1);
            % whole time increment of these two reference mesh points
            bbDt=timesIncLevy(hLocCeil-1);
            % run Brownian bridge
            [dW1, dW2]=bbridge(bbDW, bbDt, bbDt1);
            % add left portion of reference W to current W
            WincLevy(:,k)=WincLevy(:,k)+dW1;
            % update the right reference W
            Winc1MC(:,hLocCeil-1)=dW2;
            % update marker
            Marker= hLocCeil-2;
        end
        % approximating Levy area
        [parms.LA,Y2,Y4]=funcLevy(LevyStepNo,Y2,Y4,WincLevy);
        % Wiener increment to use for the scheme
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
        hfix=hfix+h;  % adaptive stepsizes
        j=j+1; % number of adaptive steps
    end
    t=t+h; % push 't'
end
outcome=StartValue;