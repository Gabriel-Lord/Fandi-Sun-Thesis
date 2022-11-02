% Main script for generating:
%     Figure 1. root-mean-sqaure error plot with timestep vs error;
%     Figure 2. efficientcy plot with error vs CPU;
% Five methods will be tested and shown:
%     1. AMil -- Adaptive Milstein  (proposed)
%     2. PMil -- fixed-step Projected Milstein
%     3. TMil -- fixed-step Tamed Milstein
%     4. SSBM -- fixed-step Split-Step Backward Milstein
%     5. TSRK -- fixed-step Tamed Stochastic Runge-Kutta of order 1.0
% Model will be in Two-Dimensional with
%     1. nonlinear drift with form in 'DFT.m'
%     2. diffusion with nonlinear multiplicative noises in 'DIFF.m'
clc; clear all; close all;

% initial value
Xzero=[3;4];

% list of model parameters called from 'ModelParameters.m'
parms=ModelParameters;
parms.rank=length(Xzero);

% parameter for controlling the adaptive step
% could be bigger when e.g. intial value is big
% see 'MainTelomerePath.m' in '1D' folder
parms.kappa=1;

%%% Levy area setting
% number of reference steps
parms.Nref=pow2(11); 
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
hmaxList=pow2([ 7 6 5 4 3])*parms.dtref;

% number of Monte Carlo realisation which is to approximate expectation
M=10;

% the noise of the model with details in 'DIFF.m' and 'DDIFF.m'
parms.ModelNoise='Levy';  

% start counting the running time of entire script
% will be displayed when the simulation is done
tStart = tic;
hmaxListLength=length(hmaxList);

% creating vector of zeros for storing data
[XerrAMil, XerrTMil, XerrPMil, XerrSSBM, XerrTSRK]...
    =deal(zeros(1,hmaxListLength));

% starting generating the data for convergence plot
for m=1:hmaxListLength  % run through the list of hmax

    % update following parameters for
    parms.hmax=hmaxList(m);
    parms.hmin=parms.hmax/parms.rho;
    parms.Tol=parms.hmin/pow2(3);

    hfix=0; j=0;

    % for each hmax, run through all the monte carlo realizations
    for i=1:M

        %-------   Reference solution -- Projected Milstein  --------------%
        Xref=Xzero; Y2=0; Y4=0;
        for k=1:parms.Nref
            % mini Wiener increments 
            % for mini step within one reference step 
            dWla(:,1:parms.Kref)=randn(parms.rank,parms.Kref)*parms.sqrtmini;
            % approximating Levy area
            [LA,Y2,Y4]=funcLevy(parms.Kref,Y2,Y4,dWla);
            % Wiener increment to use  
            Winc=sum(dWla,2);
            
            % updating 'LA' into parameter list
            parms.LA=LA;

            % calculating reference solution by stepsize dtref
            Xref=schemePMil2D(parms, Xref, Winc, parms.dtref);

            % store all Wiener increment on each MC realisation
            % for the use of other methods
            WKEEP((2*i-1):2*i, ...
                (k-1)*parms.Kref+1 : k*parms.Kref )=dWla;

        end
        % store the terminal value for each MC realisaton
        % for the use of other methods
        XrefKeep(:,i)=Xref;

        %--------------------   Adaptive Milstein  -------------------------%
        XAMil = Xzero;
        % calculating AMil with detailed explaination in 'schemeAMil.m'
        [XAMil, j, hfix]=schemeAMil2DLevy(parms, XAMil, WKEEP, i, j, hfix);
        % adding sqaured error together for the sample mean later
        XerrAMil(m)=XerrAMil(m)+norm(XAMil-XrefKeep(:,i))^2;
    end

    % Calculating hmean which will be used as
    %          the fixed stepsize for all fixed-step methods
    % 'hfix' is now the sum of adaptive steps from all MC realisations
    %          except the last step
    % 'j' is now the sum of number of adaptive step
    %          taken from all MC realisations
    hmean(m)=hfix/j;
    % size of one ministep within one fixed step dt
    deltaDt=hmean(m)^2;
    % number of reference ministeps within one ministep
    miniDtNo=floor(deltaDt/parms.minilength); 
    % resize ministep to make it a multiple 
    % of reference ministeps
    deltaDt=miniDtNo*parms.minilength; 
    % number of ministeps within one fixed step dt
    DtNo=floor(hmean(m)/deltaDt);  
    % updating dt that we use for fixed-step methods
    parms.dtuse=DtNo*deltaDt;   
    % number of total reference ministep in one dtuse
    miniTotalNo=miniDtNo*DtNo;  
    % number of dtuse fixed steps to reach T
    N=floor(parms.Nref*parms.Kref/miniTotalNo); 
    % size of the last step
    parms.dtlast=rem(parms.Nref*parms.Kref, miniTotalNo)*parms.minilength;


    % For the same 'parms.hmax'
    % run through all MC realisations for all fixed-step methods
    for i=1:M 
        % generate list of Levy area approximations first
        LAKEEP=zeros(1,N);
        Y2=0; Y4=0; Marker=0;
        for k=1:N
            dWla=zeros(parms.rank,DtNo);
            for a=1:DtNo 
                dWla(:,a)=sum( WKEEP((2*i-1):2*i,...
                    Marker+(a-1)*miniDtNo+1:Marker+a*miniDtNo) ,2);
            end
            [LAKEEP(k),Y2,Y4]=funcLevy(DtNo,Y2,Y4,dWla);
            Marker=Marker+miniTotalNo;
        end

        %-------------------   fixed-step methods  ------------------------%
        % assign initial values of all methods to 'Xzero'
        [XTMil, XPMil, XSSBM, XTSRK]=deal(Xzero);
        for k=1:N
            % pull Wiener increments from reference solution
            W=sum(WKEEP((2*i-1):2*i, ...
                (k-1)*miniTotalNo+1:k*miniTotalNo),2);
            parms.LA=LAKEEP(k);
            % run all methods
            XTMil=schemeTMil2D(parms, XTMil, W, 'fixed');
            XPMil=schemePMil2D(parms, XPMil, W, parms.dtuse);
            XSSBM=schemeSSBM2D(parms, XSSBM, W, 'fixed');
            XTSRK=schemeTSRK2D(parms, XTSRK, W, 'fixed');
        end

        % run the last step to reach 'parms.T'
        if parms.dtlast ~= 0
            W=sum(WKEEP((2*i-1):2*i, N*miniTotalNo+1:end),2);
            parms.LA=0;
            XTMil=schemeTMil2D(parms, XTMil, W, 'last');
            XPMil=schemePMil2D(parms, XPMil, W, parms.dtlast);
            XSSBM=schemeSSBM2D(parms, XSSBM, W, 'last');
            XTSRK=schemeTSRK2D(parms, XTSRK, W, 'last');
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
        %------------------   Adaptive Milstein  ---------------------%
        tic;  % start counting CPU time for AMil
        XAMil = Xzero;
        XAMil=schemeAMilCPU2DLevy(parms, XAMil);
        % end counting CPU time for AMil and store
        TimeAMil(m)=TimeAMil(m)+toc;

        %------------------   Tamed Milstein  -----------------------%
        tic; % start counting CPU time for TMil
        XTMil=Xzero; Y2=0; Y4=0;
        % Use 'hmean(m)' as the fixed stepsize for all fixed-step methods
        %               to make it fair to compare to AMil
        % It will have 'N' steps for size 'hmean(m)'
        %               and one last step for size 'parms.dtlast'
        parms.dtuse=hmean(m); % size of fixed step
        N=floor(parms.T/parms.dtuse); % number of fixed steps
        delta=parms.dtuse^2;  % size of one mini step
        K=floor(parms.dtuse/delta);  % number of mini step within one dtuse
        parms.dtlast=parms.T-parms.dtuse*N; % size of last step
        sqrtdelta=sqrt(delta);
        for k=1:N
            dWla=randn(parms.rank,K)*sqrtdelta;
            [parms.LA,Y2,Y4]=funcLevy(K,Y2,Y4,dWla);
            W=sum(dWla,2);
            XTMil=schemeTMil2D(parms, XTMil, W, 'fixed');
        end
        if parms.dtlast ~= 0
            parms.LA=0;
            W=randn(parms.rank,1)*sqrt(parms.dtlast);
            XTMil=schemeTMil2D(parms, XTMil, W, 'last');
        end
        % end counting CPU time for TMil and store
        TimeTMil(m)=TimeTMil(m)+toc;


        %----------------   Projected Milstein  ----------------------%
        tic;
        XPMil=Xzero; Y2=0; Y4=0;
        parms.dtuse=hmean(m);
        N=floor(parms.T/parms.dtuse);
        delta=parms.dtuse^2;
        K=floor(parms.dtuse/delta); % should be integer
        parms.dtlast=parms.T-parms.dtuse*N;
        sqrtdelta=sqrt(delta);
        for k=1:N
            dWla=randn(parms.rank,K)*sqrtdelta;
            [parms.LA,Y2,Y4]=funcLevy(K,Y2,Y4,dWla);
            W=sum(dWla,2);
            XPMil=schemePMil2D(parms, XPMil, W, parms.dtuse);
        end
        if parms.dtlast ~= 0
            parms.LA=0;
            W=randn(parms.rank,1)*sqrt(parms.dtlast);
            XPMil=schemePMil2D(parms, XPMil, W, parms.dtlast);
        end
        TimePMil(m)=TimePMil(m)+toc;


        %-------------------   SSBM  ------------------------------%
        tic;
        XSSBM=Xzero; Y2=0; Y4=0;
        parms.dtuse=hmean(m);
        N=floor(parms.T/parms.dtuse);
        delta=parms.dtuse^2;
        K=floor(parms.dtuse/delta); % should be integer
        parms.dtlast=parms.T-parms.dtuse*N;
        sqrtdelta=sqrt(delta);
        for k=1:N
            dWla=randn(parms.rank,K)*sqrtdelta;
            [parms.LA,Y2,Y4]=funcLevy(K,Y2,Y4,dWla);
            W=sum(dWla,2);
            XSSBM=schemeSSBM2D(parms, XSSBM, W, 'fixed');
        end
        if parms.dtlast ~= 0
            parms.LA=0;
            W=randn(parms.rank,1)*sqrt(parms.dtlast);
            XSSBM=schemeSSBM2D(parms, XSSBM, W, 'last');
        end
        TimeSSBM(m)=TimeSSBM(m)+toc;


        %------------------   TSRK -------------------------------%
        tic;
        XTSRK=Xzero; Y2=0; Y4=0;
        parms.dtuse=hmean(m);
        N=floor(parms.T/parms.dtuse);
        parms.dtlast=parms.T-parms.dtuse*N;
        sqrtdelta=sqrt(delta);
        for k=1:N
            dWla=randn(parms.rank,K)*sqrtdelta;
            [parms.LA,Y2,Y4]=funcLevy(K,Y2,Y4,dWla);
            W=sum(dWla,2);
            XTSRK=schemeTSRK2D(parms, XTSRK, W, 'fixed');
        end
        if parms.dtlast ~= 0
            parms.LA=0; 
            W=randn(parms.rank,1)*sqrt(parms.dtlast);
            XTSRK=schemeTSRK2D(parms, XTSRK, W, 'fixed');
        end
        TimeTSRK(m)=TimeTSRK(m)+toc;

    end
end
% sample mean of CPU time
TimeAMilMean=TimeAMil/M;
TimeTMilMean=TimeTMil/M;
TimePMilMean=TimePMil/M;
TimeSSBMMean=TimeSSBM/M;
TimeTSRKMean=TimeTSRK/M;


% end counting time for running the entire script
% display result
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));


%%   generating plots
figure(1)  % efficiency plot
loglog(TimeAMilMean,RMSerrAMil,'k-o','LineWidth',3,'MarkerSize',17), hold on
loglog(TimePMilMean,RMSerrPMil,'b:+','LineWidth',3,'MarkerSize',17), hold on
loglog(TimeSSBMMean,RMSerrSSBM,':x','LineWidth',3,'MarkerSize',17,'Color',[0.9290, 0.6940, 0.1250]), hold on
loglog(TimeTMilMean,RMSerrTMil,'m--*','LineWidth',3,'MarkerSize',17), hold on
loglog(TimeTSRKMean,RMSerrTSRK,'r-.s','LineWidth',3,'MarkerSize',17), hold off
set(gca,'fontsize',14)
xlabel('CPU time','FontSize',14)
ylabel('RMS error','FontSize',14)
legend({'AMil','PMil','SSBM','TMil','TSRK1'},'Location','northeast','FontSize',12)
title('(e)')
axis tight
grid


figure(2)  % convergence plot
loglog(hmean,RMSerrAMil,'k-o','LineWidth',3,'MarkerSize',17), hold on
loglog(hmean,RMSerrPMil,'b:+', 'LineWidth',3,'MarkerSize',17), hold on
loglog(hmean,RMSerrSSBM,':x','LineWidth',3,'MarkerSize',17,'Color',[0.9290, 0.6940, 0.1250]), hold on
loglog(hmean,RMSerrTMil,'m--*','LineWidth',3,'MarkerSize',17), hold on
loglog(hmean,RMSerrTSRK,'r-.s','LineWidth',3,'MarkerSize',17), hold on
loglog(hmean,hmean*exp(-2),'g-.','LineWidth',4), hold on
loglog([0.1, 0.01, 0.002],[0.1, 0.01, 0.002].*exp(-1.2),'g-.','LineWidth',3), hold off
set(gca,'fontsize',14)
xlabel('hmean','FontSize',14)
ylabel('RMS error','FontSize',16)
legend({'AMil','PMil','SSBM','TMil','TSRK1','Rate 1'},'Location','northwest','FontSize',12)
title('(f)')
axis tight
grid


tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
