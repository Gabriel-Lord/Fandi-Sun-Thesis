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
%     2. diffusion with non-commutative noises in 'DIFF.m'
%     3. linear jump amplitude in 'AMPLTD.m'
clc;clear all; close all;

% list of model parameters called from 'ModelParameters.m'
parms=ModelParameters;

% initial value
Xzero=[3;4];
rank=length(Xzero);
parms.rank=rank;

% parameter for controlling the adaptive step
% could be bigger when e.g. intial value is big
% see 'MainTelomerePath.m' in '1D' folder
parms.kappa=1;

%%% Levy area setting
% number of reference steps
parms.Nref=pow2(9);
% size of one reference step
parms.dtref=parms.T/parms.Nref;
% square root of one reference step
parms.sqrtdtref=sqrt(parms.dtref);
% size of one mini steps within one reference step
% for Levy area approximation
parms.minilength=parms.dtref^2;
% number of total mini steps within one reference step
parms.Kref=1/parms.dtref; % no. mini LA steps in 1 dtref
% square root of one mini step
parms.sqrtmini=sqrt(parms.minilength);

% list of hmax we use to see convergence
hmaxList=pow2([ 8 7 6 5 4  3])*parms.dtref; 

% number of Monte Carlo realisation which is to approximate expectation
M=50;
parms.M=M;

% the noise of the model with details in 'DIFF.m' and 'DDIFF.m'
parms.ModelNoise='Levy';

% start counting the running time of entire script
% will be displayed when the simulation is done
tStart = tic;
hmaxListLength=length(hmaxList);

% jump intensity -- can choose from e.g. 4, 10, or 20.
parms.lambda=4;

% pre-compute jump times and jump sizes for all MC realiasations
[parms.JumpTimes, parms.zeta]=JumpInfo(parms);

% pull out values
T=parms.T;


% generate the superposition of uniform reference steps and jump times
WTimes=0:parms.dtref:T;
parms.maxLength=length(WTimes)+size(parms.JumpTimes,2);
[parms.AllTimes, parms.timesInc]=allTimePoints(parms, WTimes, parms.JumpTimes);

% creating vector of zeros for storing data
[XerrAMil, XerrTMil, XerrPMil, XerrSSBM, XerrTSRK]...
    =deal(zeros(1,hmaxListLength));

% starting generating the data for convergence plot
for m=1:hmaxListLength  % run through the list of hmax
    
    % update following parameters for
    parms.hmax=hmaxList(m);
    parms.hmin=parms.hmax/parms.rho;
    Tol=parms.hmin/pow2(3);
    parms.Tol=Tol;
    
    % create empty matrix for storing Winc
    WincKeep=zeros(rank*M, parms.maxLength);
    AllTimesLevy=zeros(M, parms.maxLength);
    
    hfix=0; j=0;
    
    % for each hmax, run through all the monte carlo realizations
    for i=1:M
        
        %-------   Reference solution -- Projected Milstein  --------------%
        % i-th reference time mesh points vector without Levy steps
        timesIncI=nonzeros(parms.timesInc(i,:))'; 
        timesIncILen=length(timesIncI);      
        % copy out jump sizes and times 
        zetaTemp=parms.zeta(rank*i-(rank-1):rank*i,:);
        JumpTimesTemp=[nonzeros(parms.JumpTimes(i,:))', T+Tol];
        Winc1MC=zeros(rank,1); % reference Levy Wiener increments vector
        times1MC=0; % reference time mesh points vector with Levy steps
        Y2=0; Y4=0; % for Levy area approximation
        Xref=Xzero; % initial value
        t=0;  
        for k=1:timesIncILen
            % go through the reference time mesh points vector without Levy steps
            dt=timesIncI(k);
            % number of Levy steps in this 'dt'
            LevyStepNo=round(1/dt);   
            % size of each Levy step 
            LevyStepSize=dt/LevyStepNo;   
            % Wiener increments for the Levy steps in this 'dt'
            WincLevy=randn(rank, LevyStepNo).*sqrt(LevyStepSize); % W increments in this particular dtref
            % approximating Levy area
            [parms.LA,Y2,Y4]=funcLevy(parms.Kref,Y2,Y4,WincLevy);
            % Wiener increment to use
            Winc=sum(WincLevy,2);
            % store reference Levy Wiener increments 
            Winc1MC=[Winc1MC, WincLevy];
            % store reference time mesh points vector with Levy steps
            times1MC=[times1MC, ...
                t+LevyStepSize : LevyStepSize : t+LevyStepNo*LevyStepSize]; 
            
            % check if there is a jump and pull out size if so
            [~, ~, JumpSize, JumpTimesTemp, zetaTemp]=...
                JumpSizeFinder(t, dt, JumpTimesTemp, zetaTemp);
            % main scheme
            Xref=schemeJAdptPMil2D(parms, Xref, Winc, dt, JumpSize);
            % update t
            t=t+dt;
        end
        % store the terminal value for each MC realisaton
        % for the use of other methods
        XrefKeep(:,i)=Xref;
        % store reference Levy Wiener increments for fixed-step methods
        WincKeep(rank*i-(rank-1):rank*i,...
            1:length(Winc1MC)-1)=Winc1MC(:,2:end);  
        % store reference time points with Levy steps for fixed-step methods
        AllTimesLevy(i,1:length(times1MC))=times1MC;  
        % for AMil use
        parms.times1MC=times1MC;
        parms.timesIncLevy=diff(times1MC);
        
        %--------------------   Adaptive Milstein  -------------------------%
        XAMil = Xzero;
        % calculating AMil with detailed explaination in 'schemeAMil.m'
        [XAMil, j, hfix]=schemeJAdptAMil2DLevy(parms, XAMil, Winc1MC(:,2:end), i, j, hfix);
        % adding sqaured error together for the sample mean later
        XerrAMil(m)=XerrAMil(m)+norm(XAMil-XrefKeep(:,i))^2;
    end
    
    % the averaged adaptive steps
    % hmean(m)=hfix/j; % old method, abandoned
    % new method, removed all jump points inherited from adapative mesh
    hmean(m)=hfix/(j-length(nonzeros(parms.JumpTimes))+M); 
    % fixed stepsize to use 
    Dt=hmean(m);
    % generate the superposition of uniform reference steps and jump times
    WTimesDt=0:Dt:T; 
    [parms.AllTimesDt, parms.timesIncDt]...
        =allTimePoints(parms, WTimesDt, parms.JumpTimes);
    
    % For the same 'parms.hmax'
    % run through all MC realisations for all fixed-step methods
    for i=1:M
        % i-th reference time mesh points vector without Levy steps
        timesIncDtI=nonzeros(parms.timesIncDt(i,:))';
        timesIncDtILen=length(timesIncDtI);
        % copy out jump sizes and times 
        JumpTimesTemp=[nonzeros(parms.JumpTimes(i,:))',2*T];
        zetaTemp=parms.zeta(rank*i-(rank-1):rank*i,:);
        % i-th time points with reference Levy steps 
        times1MCDt=nonzeros(AllTimesLevy(i,:))';
        % i-th reference Levy Wiener increments  
        Winc1MCDt=WincKeep(rank*i-(rank-1):rank*i,:);
        % i-th time increments with reference Levy steps
        timesIncLevy=diff(times1MCDt);
        %-------------------   fixed-step methods  ------------------------%
        % assign initial values of all methods to 'Xzero'
        [XTMil, XPMil, XSSBM, XTSRK]=deal(Xzero);
        Y2=0; Y4=0; % for Levy area approximation
        Marker=0; % to bookmark the starting position of each step
        t=0;  
        for k=1:timesIncDtILen
            dt=timesIncDtI(k);
            % number of Levy steps in this 'dt'
            LevyStepNo=round(1/dt);   
            % size of each Levy step 
            LevyStepSize=dt/LevyStepNo;   
            % check if there is a jump in this step and look for its jump size
            [~, ~, JumpSize, JumpTimesTemp, zetaTemp]=...
                JumpSizeFinder(t, dt, JumpTimesTemp, zetaTemp);
            % create empty vector for storeing Levy Wiener increments
            WincLevy=zeros(rank,LevyStepNo);
            % following loop is for creating Levy Wiener increments
            % with Brownian bridge on the reference line
            for r=1:LevyStepNo
                % difference with the reference line 
                times1MCDiff=t+r*LevyStepSize-times1MCDt;   
                % index of the right mesh point of ''t+dt''
                hLocCeil=find(times1MCDiff<=0,1); 
                % Wiener increment up to the left mesh point of ''t+h''
                WincLevy(:,r)=sum(Winc1MCDt(:, Marker+1 : hLocCeil-2),2);
                                    % left portion of this step within two reference mesh points
                bbDt1= t+r*LevyStepSize - times1MCDt(hLocCeil-1); 
                                    % whole Wiener increment of these two reference mesh points
                bbDW=Winc1MCDt(:,hLocCeil-1);  
                                    % whole time increment of these two reference mesh points
                bbDt=timesIncLevy(hLocCeil-1);  
                                    % run Brownian bridge
                [dW1, dW2]=bbridge(bbDW, bbDt, bbDt1);
                                    % add left portion of reference W to current W
                WincLevy(:,r)=WincLevy(:,r)+dW1; 
                                    % update the right reference W
                Winc1MCDt(:,hLocCeil-1)=dW2;  
                Marker= hLocCeil-2;
            end
            % approximating Levy area
            [parms.LA,Y2,Y4]=funcLevy(LevyStepNo,Y2,Y4,WincLevy);
            % Wiener increment to use for the scheme
            WincDt=sum(WincLevy,2);
            % run all fixed-step methods
            XTMil=schemeJAdptTMil2D(parms, XTMil, WincDt, dt, JumpSize);
            XPMil=schemeJAdptPMil2D(parms, XPMil, WincDt, dt, JumpSize);
            XSSBM=schemeJAdptSSBM2D(parms, XSSBM, WincDt, dt, JumpSize);
            XTSRK=schemeJAdptTSRK2D(parms, XTSRK, WincDt, dt, JumpSize);
            t=t+dt;
        end
   
        % adding sqaured error together for the sample mean later
        XerrTMil(m)=XerrTMil(m)+norm(XTMil-XrefKeep(:,i))^2;
        XerrPMil(m)=XerrPMil(m)+norm(XPMil-XrefKeep(:,i))^2;
        XerrSSBM(m)=XerrSSBM(m)+norm(XSSBM-XrefKeep(:,i))^2;
        XerrTSRK(m)=XerrTSRK(m)+norm(XTSRK-XrefKeep(:,i))^2;
        
    end
end

% sample mean of root-mean-square
RMSerrAMil=sqrt(XerrAMil/M);
RMSerrTMil=sqrt(XerrTMil/M);
RMSerrPMil=sqrt(XerrPMil/M);
RMSerrSSBM=sqrt(XerrSSBM/M);
RMSerrTSRK=sqrt(XerrTSRK/M);

%%    starting generating the data for efficiency plot

% In this section, we are only interested in the CPU time
%               so there will be no reference solution
% Each method will be tested  separately
%               with  steps recorded in the most fair way

% creating vector of zeros for storing data
[TimeAMil, TimeTMil, TimePMil, TimeSSBM, TimeTSRK]...
    =deal(zeros(1,hmaxListLength));
for m=1:hmaxListLength  % run through the list of hmax
    parms.hmax=hmaxList(m);
    parms.hmin=parms.hmax/parms.rho;
    parms.Tol=parms.hmin/pow2(3);
    for i=1:M  % run through all the MC realisatons
        %------------   jump-adapted Adaptive Milstein  --------------%
        tAMil=tic;  % start counting CPU time for AMil
        XAMil = Xzero;
        XAMil=schemeJAdptAMilCPU2DLevy(parms, i, XAMil);
        % end counting CPU time for AMil and store
        TimeAMil(m)=TimeAMil(m)+toc(tAMil);
        
        
        % take hmean as the fixed step size
        Dt=hmean(m);
        % the max number of hmean that T can take
        DtMultiple=floor(T/Dt);
        
        %------------  jump-adapted Tamed Milstein  ----------------%
        tTMil=tic; % start counting CPU time for TMil
        % create uniform mesh point
        WTimesDt=[0:Dt:DtMultiple*Dt, T];
        % copy out jump times and sizes
        JumpTimesI=nonzeros(parms.JumpTimes(i,:))';
        JumpTimesTemp=[JumpTimesI,2*T];
        zetaTemp=parms.zeta(rank*i-(rank-1):rank*i,:);
        % combine and sort uniform times and jump times
        AllTimesDt=sort([WTimesDt, JumpTimesI]);
        % remove all repeated times in t line
        AllTimesDt(diff(AllTimesDt) == 0)=[];
        % time increments of all mesh points
        timesIncDt=diff(AllTimesDt); 
        % initial value
        XTMil=Xzero; t=0; Y2=0; Y4=0;
        for k=1:length(timesIncDt)
            dt=timesIncDt(k);
            % check if there is a jump and pull out size if so
            [~, ~, JumpSize, JumpTimesTemp, zetaTemp]=...
                JumpSizeFinder(t, dt, JumpTimesTemp, zetaTemp);
            % setting for Levy area
            LevyStepNo=round(1/dt);   
            LevyStepSize=dt/LevyStepNo;   
            WincLevy=randn(rank,LevyStepNo).*sqrt(LevyStepSize);
            WincDt=sum(WincLevy,2);
            % Levy area approximation
            [parms.LA,Y2,Y4]=funcLevy(LevyStepNo,Y2,Y4,WincLevy); 
            % main scheme 
            XTMil=schemeJAdptTMil2D(parms, XTMil, WincDt, dt, JumpSize);
            t=t+dt;
        end 
        % end counting CPU time for TMil and store
        TimeTMil(m)=TimeTMil(m)+toc(tTMil);
        
        
        %---------  jump-adapted Projected Milstein  -----------------%
        tPMil=tic; % start counting CPU time for PMil
        % create uniform mesh point
        WTimesDt=[0:Dt:DtMultiple*Dt, T];
        % copy out jump times and sizes
        JumpTimesI=nonzeros(parms.JumpTimes(i,:))';
        JumpTimesTemp=[JumpTimesI,2*T];
        zetaTemp=parms.zeta(rank*i-(rank-1):rank*i,:);
        % combine and sort uniform times and jump times
        AllTimesDt=sort([WTimesDt, JumpTimesI]);
        % remove all repeated times in t line
        AllTimesDt(diff(AllTimesDt) == 0)=[];
        % time increments of all mesh points
        timesIncDt=diff(AllTimesDt); 
        % initial value
        XPMil=Xzero; t=0; Y2=0; Y4=0;
        for k=1:length(timesIncDt)
            dt=timesIncDt(k);
            % check if there is a jump and pull out size if so
            [~, ~, JumpSize, JumpTimesTemp, zetaTemp]=...
                JumpSizeFinder(t, dt, JumpTimesTemp, zetaTemp);
            % setting for Levy area
            LevyStepNo=round(1/dt);   
            LevyStepSize=dt/LevyStepNo;   
            WincLevy=randn(rank,LevyStepNo).*sqrt(LevyStepSize);
            WincDt=sum(WincLevy,2);
            % Levy area approximation
            [parms.LA,Y2,Y4]=funcLevy(LevyStepNo,Y2,Y4,WincLevy); 
            % main scheme 
            XPMil=schemeJAdptPMil2D(parms, XPMil, WincDt, dt, JumpSize);
            t=t+dt;
        end 
        % end counting CPU time for PMil and store
        TimePMil(m)=TimePMil(m)+toc(tPMil);
        
        
        %---------------  jump-adapted SSBM  -----------------------%
        tSSBM=tic;  % start counting CPU time for SSBM
        % create uniform mesh point
        WTimesDt=[0:Dt:DtMultiple*Dt, T];
        % copy out jump times and sizes
        JumpTimesI=nonzeros(parms.JumpTimes(i,:))';
        JumpTimesTemp=[JumpTimesI,2*T];
        zetaTemp=parms.zeta(rank*i-(rank-1):rank*i,:);
        % combine and sort uniform times and jump times
        AllTimesDt=sort([WTimesDt, JumpTimesI]);
        % remove all repeated times in t line
        AllTimesDt(diff(AllTimesDt) == 0)=[];
        % time increments of all mesh points
        timesIncDt=diff(AllTimesDt); 
        % initial value
        XSSBM=Xzero; t=0; Y2=0; Y4=0;
        for k=1:length(timesIncDt)
            dt=timesIncDt(k);
            % check if there is a jump and pull out size if so
            [~, ~, JumpSize, JumpTimesTemp, zetaTemp]=...
                JumpSizeFinder(t, dt, JumpTimesTemp, zetaTemp);
            % setting for Levy area
            LevyStepNo=round(1/dt);   
            LevyStepSize=dt/LevyStepNo;   
            WincLevy=randn(rank,LevyStepNo).*sqrt(LevyStepSize);
            WincDt=sum(WincLevy,2);
            % Levy area approximation
            [parms.LA,Y2,Y4]=funcLevy(LevyStepNo,Y2,Y4,WincLevy); 
            % main scheme 
            XSSBM=schemeJAdptSSBM2D(parms, XSSBM, WincDt, dt, JumpSize);
            t=t+dt;
        end 
        % end counting CPU time for SSBM and store
        TimeSSBM(m)=TimeSSBM(m)+toc(tSSBM);
        
        
        %--------------  jump-adapted TSRK -------------------------%
        tTSRK=tic; % start counting CPU time for TSRK
        % create uniform mesh point
        WTimesDt=[0:Dt:DtMultiple*Dt, T];
        % copy out jump times and sizes
        JumpTimesI=nonzeros(parms.JumpTimes(i,:))';
        JumpTimesTemp=[JumpTimesI,2*T];
        zetaTemp=parms.zeta(rank*i-(rank-1):rank*i,:);
        % combine and sort uniform times and jump times
        AllTimesDt=sort([WTimesDt, JumpTimesI]);
        % remove all repeated times in t line
        AllTimesDt(diff(AllTimesDt) == 0)=[];
        % time increments of all mesh points
        timesIncDt=diff(AllTimesDt); 
        % initial value
        XTSRK=Xzero; t=0; Y2=0; Y4=0;
        for k=1:length(timesIncDt)
            dt=timesIncDt(k);
            % check if there is a jump and pull out size if so
            [~, ~, JumpSize, JumpTimesTemp, zetaTemp]=...
                JumpSizeFinder(t, dt, JumpTimesTemp, zetaTemp);
            % setting for Levy area
            LevyStepNo=round(1/dt);   
            LevyStepSize=dt/LevyStepNo;   
            WincLevy=randn(rank,LevyStepNo).*sqrt(LevyStepSize);
            WincDt=sum(WincLevy,2);
            % Levy area approximation
            [parms.LA,Y2,Y4]=funcLevy(LevyStepNo,Y2,Y4,WincLevy); 
            % main scheme 
            XTSRK=schemeJAdptTSRK2D(parms, XTSRK, WincDt, dt, JumpSize);
            t=t+dt;
        end 
        % end counting CPU time for TSRK and store
        TimeTSRK(m)=TimeTSRK(m)+toc(tTSRK);
        
    end
end
% sample mean of CPU time
TimeAMilMean=TimeAMil/M;
TimeTMilMean=TimeTMil/M;
TimePMilMean=TimePMil/M;
TimeSSBMMean=TimeSSBM/M;
TimeTSRKMean=TimeTSRK/M;


%end counting time for running the entire script 
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));


%% generating plots 
x=2;
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





