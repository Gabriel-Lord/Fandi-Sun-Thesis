function [outcome, j, hfix]=schemeJAdptAMil(parms, StartValue, Winc, i, j, hfix)
% % To generate the terminal value of the proposed Adaptive Milstein method
% %           based on a reference solution
% % INPUTS:
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts
% % 'StartValue' -- (scalar) initial value of the process
% % 'Winc' -- (1-by-n vector) Wiener increment from the reference solution
% % 'i' -- (scalar) the index of MC realisation
% % 'j' -- (scalar) number of adaptive steps taken on all MC realisations
% % 'hfix' -- (scalar) sum length of adaptive steps taken on all MC realisations
% % OUTPUTS:
% % 'outcome' -- (scalar) terminal value of the process
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

WincTemp=Winc; 

% pull out jump sizes, jump times, reference times 
% and reference time increments
zetaTemp=parms.zeta(i,:);
JumpTimesTemp=[nonzeros(parms.JumpTimes(i,:))', parms.T+parms.Tol];
AllTimesI=parms.AllTimes(i,:);
timesIncI=parms.timesInc(i,:);

% we add each adaptive step to 't'
% until it reaches terminal time 'parms.T'
t=0;
while abs(t-parms.T)>parms.Tol && t<=parms.T
    
    % core function of choosing 'h'
    HcoreFunc=AMilCoreH(StartValue,parms);
    % bounding the function value by hmin and hmax
    h=max(parms.hmin, min(parms.hmax, HcoreFunc));
 
    if t+h > parms.T  % last tep
        h=parms.T-t;
    end
    
    % check if there is a jump in this step and look for its jump size 
    [h, jumpIndct, JumpSize, JumpTimesTemp, zetaTemp]=...
                        JumpSizeFinder(t, h, JumpTimesTemp, zetaTemp);
     
    % check location of h in t line
    AllTimesDiff=t+h-AllTimesI;  % from reference line
    hLocCeil=find(AllTimesDiff<0,1); % index of the right mesh point of ''t+h''
    
    % Wiener increment up to the left mesh point of ''t+h''
    WincAMil=sum(WincTemp(Marker+1 : hLocCeil-2));
    
    % check if h lands between 2 times, if so, need to Brownian bridge
    if jumpIndct == 0  % true if it is not landing at jump mesh point
        if  isempty(hLocCeil) ~=  true  % if found
            % left portion of this step within two reference mesh points
            bbDt1= t+h - AllTimesI(hLocCeil-1);
            % whole Wiener increment of these two reference mesh points
            bbDW=WincTemp(hLocCeil-1);  
            % whole time increment of these two reference mesh points
            bbDt=timesIncI(hLocCeil-1);  
            % run Brownian bridge 
            [dW1, dW2]=bbridge(bbDW, bbDt, bbDt1);
            % add left portion of reference W to current W
            WincAMil=WincAMil+dW1; 
            % update the right reference W 
            WincTemp(hLocCeil-1)=dW2; % add right portion to next W, so it will be smaller
        else  % last step
            WincAMil=sum(WincTemp(Marker+1 : end));
        end
    end 

    % update Marker
    Marker= hLocCeil-2;
     
    if h <= parms.hmin % backstop case
        % we use jump-adapted PMil with timestep 'h' for the backstop
        StartValue=schemeJAdptPMil(parms, StartValue, WincAMil, h, JumpSize);
    else
        % explicit jump-adapted Milstein
        [dft,Diff,Ddiff]=coefficients(StartValue, parms);
        StartValueMid=StartValue+dft*h+diffTerms(WincAMil,Diff,Ddiff,h);
        StartValue=StartValueMid+AMPLTD(JumpSize, StartValueMid);
        hfix=hfix+h;  % adaptive stepsizes
        j=j+1; % number of adaptive steps
    end
    t=t+h; % push 't'
end
outcome=StartValue;