% Main script for generating:
%     Figure 1. root-mean-sqaure error plot with timestep vs error;
%     Figure 2. efficientcy plot with error vs CPU;
% Five methods will be tested and shown:
%     1. JA-AMil -- jump adapted Adaptive Milstein  (proposed)
%     2. JA-PMil -- fixed-step jump adapted Projected Milstein
%     3. JA-TMil -- fixed-step jump adapted Tamed Milstein
%     4. JA-SSBM -- fixed-step jump adapted Split-Step Backward Milstein
%     5. JA-TSRK -- fixed-step jump adapted Tamed Stochastic Runge-Kutta of order 1.0
% Model will be in Two-Dimensional with
%     1. nonlinear drift with form in 'DFT.m'
%     2. diffusion with additive, multiplicative
%                     and nonlinear multiplicative noises in 'DIFF.m'
%     3. linear jump amplitude in 'AMPLTD.m'
clc;clear all; close all;

% list of model parameters called from 'ModelParameters.m'
parms=ModelParameters;

% initial value
Xzero=[5,7]';
parms.rank=length(Xzero);

% parameter for controlling the adaptive step
% could be bigger when e.g. intial value is big
% see 'MainTelomerePath.m'
parms.kappa=1;

% list of hmax we use to see convergence
hmaxList=pow2([  8 7 6 5 4 ])*parms.dtref;
hmaxListLength=length(hmaxList);

% number of Monte Carlo realisation which is to approximate expectation
M=5;
parms.M=M;

% the noise of the model with details in 'DIFF.m' and 'DDIFF.m'
% with the choices of  'diagonal' and 'commutative'. 
parms.ModelNoise='commutative';  %  'diagonal'  %  'commutative' 


% start counting the running time of entire script
% will be displayed when the simulation is done
tStart = tic;

% pull out terminal time for convenience
T=parms.T;

% jump intensity -- can choose from e.g. 4, 10, or 20.
parms.lambda=4;

% pre-compute jump times and jump sizes for all MC realiasations
[parms.JumpTimes, parms.zeta]=JumpInfo(parms);

% generate the superposition of uniform reference steps and jump times
WTimes=0:parms.dtref:T;
parms.maxLength=length(WTimes)+size(parms.JumpTimes,2);
[parms.AllTimes, parms.timesInc]=allTimePoints(parms, WTimes, parms.JumpTimes);

% creating vector of zeros for storing data
[XerrAMil, XerrTMil, XerrPMil, XerrSSBM, XerrTSRK]...
    =deal(zeros(1,hmaxListLength));

%% convergence
for m=1:hmaxListLength % go through hmaxlist list
    % update following parameters for
    parms.hmax=hmaxList(m);
    parms.hmin=parms.hmax/parms.rho;
    Tol=parms.hmin/pow2(3);
    parms.Tol=Tol;
    
    % hfix and j are for recording the adaptive steps
    % in all MC simulation
    hfix=0; j=0;
    % create empty matrix for storing Wiener increments
    WincKeep=zeros(parms.rank*M,parms.maxLength);
    
    for i=1:M  % run through all MC realisations
        
        %---------- Reference solution - jump-adapted Projected Milstein ---------- %
        % copy out values
        timesIncI=nonzeros(parms.timesInc(i,:))';
        timesIncILen=length(timesIncI);
        zetaTemp=parms.zeta(parms.rank*i-(parms.rank-1):parms.rank*i,:);
        JumpTimesTemp=[nonzeros(parms.JumpTimes(i,:))', T+Tol];
        % generate reference Wiener increments
        Winc=randn(parms.rank,timesIncILen).*sqrt(timesIncI);
        % store reference Wiener increments
        WincKeep(parms.rank*i-(parms.rank-1):parms.rank*i,1:size(Winc,2))=Winc;
        Xref=Xzero;  t=0;
        for k=1:timesIncILen  % run through reference mesh points
            Wref=Winc(:,k);
            dt=timesIncI(k);
            % check if there is a jump and pull out size if so
            [~, ~, JumpSize, JumpTimesTemp, zetaTemp]=...
                JumpSizeFinder(t, dt, JumpTimesTemp, zetaTemp);
            % main scheme
            Xref=schemeJAdptPMil2D(parms, Xref, Wref, dt, JumpSize);
            t=t+dt;
        end
        % store reference terminal value
        XrefKeep(:,i)=Xref;
        
        %-------------------- jump-adapted Adaptive Milstein -------------------- %
        XAMil = Xzero;
        % calculating AMil with detailed explaination in 'schemeAMil.m'
        [XAMil, j, hfix]=schemeJAdptAMil2D(parms, XAMil, Winc, i, j, hfix);
        % adding sqaured error together for the sample mean later
        XerrAMil(m)=XerrAMil(m)+norm(XAMil-XrefKeep(:,i))^2;
    end
    
    % the averaged adaptive steps
    hmean(m)=hfix/j;
    % the max number of hmean that T can take as a whole
    hmeanMultiple=floor(T/hmean(m));
    % uniform mesh points
    WTimesDt=[0:hmean(m):hmeanMultiple*hmean(m), T];
    % superposition of 'hmean's and jump times
    % for all MC realisations
    [~, timesIncDt]=allTimePoints(parms, WTimesDt, parms.JumpTimes);
    
    for i=1:M % run through MC realisations for all fixed-step methods
        % copy out values
        timesIncDtI=nonzeros(timesIncDt(i,:))';
        timesIncDtLen=length(timesIncDtI);
        JumpTimesTemp=[nonzeros(parms.JumpTimes(i,:))',T+Tol];
        zetaTemp=parms.zeta(parms.rank*i-(parms.rank-1):parms.rank*i,:);
        WincRefTemp=WincKeep(parms.rank*i-(parms.rank-1):parms.rank*i,:);
        AllTimesI=parms.AllTimes(i,:);
        timesIncI=nonzeros(parms.timesInc(i,:))';
        
        % create empty vectors for storing data later
        [XTMil, XPMil, XSSBM, XTSRK]=deal(Xzero);
        Marker=0; t=0;
        for k=1:timesIncDtLen % run through mesh points
            dt=timesIncDtI(k);
            % check if there is a jump and pull out size if so
            [~, jumpIndct, JumpSize, JumpTimesTemp, zetaTemp]=...
                JumpSizeFinder(t, dt, JumpTimesTemp, zetaTemp);
            % find the location of the landing point on the reference line
            diff_AllTimes=t+dt-AllTimesI;
            hLocCeil=find(diff_AllTimes<0,1); % index of the right mesh point of ''t+dt''
            % Wiener increment up to the left mesh point of ''t+dt''
            WincDt=sum(WincRefTemp(:,Marker+1 : hLocCeil-2),2);
            % check if h lands between 2 times, if so, need to Brownian bridge
            if jumpIndct == 0  % true if it is not landing at jump mesh point
                if  isempty(hLocCeil) ~=  true  % if found
                    % left portion of this step within two reference mesh points
                    bbDt1= t+dt - AllTimesI(hLocCeil-1);
                    % whole Wiener increment of these two reference mesh points
                    bbDW=WincRefTemp(:,hLocCeil-1);
                    % whole time increment of these two reference mesh points
                    bbDt=timesIncI(hLocCeil-1);
                    % run Brownian bridge
                    [dW1, dW2]=bbridge(bbDW, bbDt, bbDt1);
                    % add left portion of reference W to current W
                    WincDt=WincDt+dW1;
                    % update the right reference W
                    WincRefTemp(:,hLocCeil-1)=dW2;
                else  % last step
                    WincDt=sum(WincRefTemp(:, Marker+1 : end),2);
                end
            end
            % update marker to be the left mesh point of 't+dt'
            Marker= hLocCeil-2;
            
            % run all fixed-step methods
            XTMil=schemeJAdptTMil2D(parms, XTMil, WincDt, dt, JumpSize);
            XPMil=schemeJAdptPMil2D(parms, XPMil, WincDt, dt, JumpSize);
            XSSBM=schemeJAdptSSBM2D(parms, XSSBM, WincDt, dt, JumpSize);
            XTSRK=schemeJAdptTSRK2D(parms, XTSRK, WincDt, dt, JumpSize);
            
            t=t+dt;
        end
        % adding sqaured error together for the sample mean later
        XerrPMil(m)=XerrPMil(m)+norm(XrefKeep(:,i)-XPMil)^2;
        XerrTMil(m)=XerrTMil(m)+norm(XrefKeep(:,i)-XTMil)^2;
        XerrSSBM(m)=XerrSSBM(m)+norm(XrefKeep(:,i)-XSSBM)^2;
        XerrTSRK(m)=XerrTSRK(m)+norm(XrefKeep(:,i)-XTSRK)^2;
        
    end
end

% sample mean of root-mean-square
RMSerrAMil=sqrt(XerrAMil./M);
RMSerrTMil=sqrt(XerrTMil./M);
RMSerrPMil=sqrt(XerrPMil./M);
RMSerrSSBM=sqrt(XerrSSBM./M);
RMSerrTSRK=sqrt(XerrTSRK./M);

%%  efficiency

% creating vector of zeros for storing data
[TimeAMil, TimeTMil, TimePMil, TimeSSBM, TimeTSRK ]=deal(zeros(1,hmaxListLength));

for m=1:hmaxListLength % go through hmaxlist list
    parms.hmax=hmaxList(m);
    parms.hmin=parms.hmax/parms.rho;
    Tol=parms.hmin/pow2(3);
    parms.Tol=Tol;
    for i=1:M
        
        %-------------- jump-adapted Adaptive Milstein -------------- %
        tAMil=tic;
        XAMil = Xzero;
        XAMil=schemeJAdptAMilCPU2D(parms, i, XAMil);
        % end counting CPU time for AMil and store
        TimeAMil(m)= TimeAMil(m)+toc(tAMil);
        
        
        % take hmean as the fixed step size
        Dt=hmean(m);
        % the max number of hmean that T can take
        DtMultiple=floor(T/Dt);
        %-------------- jump-adapted Projected Milstein -------------- %
        tPMil=tic;
        % create uniform mesh point
        WTimesDt=[0:Dt:DtMultiple*Dt, T];
        % copy out jump times and sizes
        JumpTimesI=nonzeros(parms.JumpTimes(i,:))';
        JumpTimesTemp=[JumpTimesI,T+Tol];
        zetaTemp=parms.zeta(parms.rank*i-(parms.rank-1):parms.rank*i,:);
        % combine and sort uniform times and jump times
        AllTimesDt=sort([WTimesDt, JumpTimesI]);
        % remove all repeated times in t line
        AllTimesDt(diff(AllTimesDt) == 0)=[];
        % time increments of all mesh points
        timesIncDt=diff(AllTimesDt);
        timesIncDtLen=length(timesIncDt);
        % initial value
        XPMil=Xzero;
        t=0;
        for k=1:timesIncDtLen   % run through combined mesh points
            dt=timesIncDt(k);
            WincDt=randn(parms.rank,1)*sqrt(dt);
            % check if there is a jump and pull out size if so
            [~, ~, JumpSize, JumpTimesTemp, zetaTemp]=...
                JumpSizeFinder(t, dt, JumpTimesTemp, zetaTemp);
            % main scheme
            XPMil=schemeJAdptPMil2D(parms, XPMil, WincDt, dt, JumpSize);
            t=t+dt;
        end
        % end counting CPU time for PMil and store
        TimePMil(m)= TimePMil(m)+toc(tPMil);
        
        
        %-------------- jump-adapted Tamed Milstein -------------- %
        t_Tmd=tic;
        WTimesDt=[0:Dt:DtMultiple*Dt, T];
        JumpTimesI=nonzeros(parms.JumpTimes(i,:))';
        JumpTimesTemp=[JumpTimesI,T+Tol];
        zetaTemp=parms.zeta(parms.rank*i-(parms.rank-1):parms.rank*i,:);
        AllTimesDt=sort([WTimesDt, JumpTimesI]);
        AllTimesDt(diff(AllTimesDt) == 0)=[];
        timesIncDt=diff(AllTimesDt);
        timesIncDtLen=length(timesIncDt);
        XTMil=Xzero;
        t=0;
        for k=1:timesIncDtLen
            dt=timesIncDt(k);
            [~, ~, JumpSize, JumpTimesTemp, zetaTemp]=...
                JumpSizeFinder(t, dt, JumpTimesTemp, zetaTemp);
            WincDt=randn(parms.rank,1)*sqrt(dt);
            XTMil=schemeJAdptTMil2D(parms, XTMil, WincDt, dt, JumpSize);
            t=t+dt;
        end
        TimeTMil(m)= TimeTMil(m)+toc(t_Tmd);
        
        %-------------- jump-adapted SSBM -------------- %
        t_SSBM=tic;
        WTimesDt=[0:Dt:DtMultiple*Dt, T];
        JumpTimesI=nonzeros(parms.JumpTimes(i,:))';
        JumpTimesTemp=[JumpTimesI,T+Tol];
        zetaTemp=parms.zeta(parms.rank*i-(parms.rank-1):parms.rank*i,:);
        AllTimesDt=sort([WTimesDt, JumpTimesI]);
        AllTimesDt(diff(AllTimesDt) == 0)=[];
        timesIncDt=diff(AllTimesDt);
        timesIncDtLen=length(timesIncDt);
        XSSBM=Xzero;
        t=0;
        for k=1:timesIncDtLen
            dt=timesIncDt(k);
            [~, ~, JumpSize, JumpTimesTemp, zetaTemp]=...
                JumpSizeFinder(t, dt, JumpTimesTemp, zetaTemp);
            WincDt=randn(parms.rank,1)*sqrt(dt);
            XSSBM=schemeJAdptSSBM2D(parms, XSSBM, WincDt, dt, JumpSize);
            t=t+dt;
        end
        TimeSSBM(m)= TimeSSBM(m)+toc(t_SSBM);
        
        %-------------- jump-adapted TSRK -------------- %
        t_TSRK=tic;
        WTimesDt=[0:Dt:DtMultiple*Dt, T];
        JumpTimesI=nonzeros(parms.JumpTimes(i,:))';
        JumpTimesTemp=[JumpTimesI,T+Tol];
        zetaTemp=parms.zeta(parms.rank*i-(parms.rank-1):parms.rank*i,:);
        AllTimesDt=sort([WTimesDt, JumpTimesI]);
        AllTimesDt(diff(AllTimesDt) == 0)=[];
        timesIncDt=diff(AllTimesDt);
        timesIncDtLen=length(timesIncDt);
        XTSRK=Xzero;
        t=0;
        for k=1:timesIncDtLen
            dt=timesIncDt(k);
            [~, ~, JumpSize, JumpTimesTemp, zetaTemp]=...
                JumpSizeFinder(t, dt, JumpTimesTemp, zetaTemp);
            WincDt=randn(parms.rank,1)*sqrt(dt);
            XTSRK=schemeJAdptTSRK2D(parms, XTSRK, WincDt, dt, JumpSize);
            t=t+dt;
        end
        TimeTSRK(m)= TimeTSRK(m)+toc(t_TSRK);
        
    end
end

% sample mean of CPU time
TimeAMilMean=TimeAMil./M;
TimeTMilMean=TimeTMil./M;
TimePMilMean=TimePMil./M;
TimeSSBMMean=TimeSSBM./M;
TimeTSRKMean=TimeTSRK./M;


%%  generating plots
x=1;
figure(1)  % convergence  
loglog(hmean(x:end),RMSerrAMil(x:end),'k-o','LineWidth',3,'MarkerSize',17), hold on
loglog(hmean(x:end),RMSerrPMil(x:end),'b:+', 'LineWidth',3,'MarkerSize',17), hold on
loglog(hmean(x:end),RMSerrSSBM(x:end),':x','LineWidth',3,'MarkerSize',17,'Color',[0.9290, 0.6940, 0.1250]), hold on
loglog(hmean(x:end),RMSerrTMil(x:end),'m--*','LineWidth',3,'MarkerSize',17), hold on
loglog(hmean(x:end),RMSerrTSRK(x:end),'r-.s','LineWidth',3,'MarkerSize',17), hold on 
loglog(hmean(x:end),hmean(x:end)*exp(0),'g-.','LineWidth',4), hold off
set(gca,'fontsize',14)
xlabel('hmean','FontSize',14) 
ylabel('RMS error','FontSize',16)
legend({'JA-AM','JA-PMil','JA-SSBM','JA-TMil','JA-TSRK1','Rate 1'},'Location','northwest','FontSize',12) 
title('(a)')
axis tight
grid
 
figure(2) % efficiency 
loglog(TimeAMil(x:end),RMSerrAMil(x:end),'k-o','LineWidth',3,'MarkerSize',17), hold on
loglog(TimePMil(x:end),RMSerrPMil(x:end),'b:+','LineWidth',3,'MarkerSize',17), hold on
loglog(TimeSSBM(x:end),RMSerrSSBM(x:end),':x','LineWidth',3,'MarkerSize',17,'Color',[0.9290, 0.6940, 0.1250]), hold on
loglog(TimeTMil(x:end),RMSerrTMil(x:end),'m--*','LineWidth',3,'MarkerSize',17), hold on
loglog(TimeTSRK(x:end),RMSerrTSRK(x:end),'r-.s','LineWidth',3,'MarkerSize',17), hold off 
set(gca,'fontsize',14)
xlabel('CPU time','FontSize',14)
ylabel('RMS error','FontSize',14)
legend({'JA-AM','JA-PMil','JA-SSBM','JA-TMil','JA-TSRK1'},'Location','northeast','FontSize',12) 
title('(b)')
axis tight
grid

% print clock time that this script takes
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));








