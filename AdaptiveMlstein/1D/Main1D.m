% Main script for generating: 
%     Figure 1. root-mean-sqaure error plot with timestep vs error; 
%     Figure 2. efficientcy plot with error vs CPU; 
% Five methods will be tested and shown: 
%     1. AMil -- Adaptive Milstein  (proposed)
%     2. PMil -- fixed-step Projected Milstein
%     3. TMil -- fixed-step Tamed Milstein
%     4. SSBM -- fixed-step Split-Step Backward Milstein
%     5. TSRK -- fixed-step Tamed Stochastic Runge-Kutta of order 1.0
% Model will be in One-Dimensional with 
%     1. nonlinear drift with form in 'DFT.m'
%     2. diffusion with additive, multiplicative 
%                     and nonlinear multiplicative noises in 'DIFF.m'

clc;clear all; close all;

% initial value
Xzero=11; 

% list of model parameters called from 'ModelParameters.m'
parms=ModelParameters;  

% parameter for controlling the adaptive step
% could be bigger when e.g. intial value is big
% see 'MainTelomerePath.m'
parms.kappa=1; 

% list of hmax we use to see convergence
hmaxList=pow2([12  10  8  6  4])*parms.dtref; 

% number of Monte Carlo realisation which is to approximate expectation
M=3; 

% the noise of the model with details in 'DIFF.m' and 'DDIFF.m'
% with the choices of  'additive'; 'mul1'; and 'mul2'.
parms.ModelNoise='mul2';  %  'additive'  %  'mul1'  % 'mul2'

%  function needed for TSRK, for solving the cubic function in drift
%  will need to update if drift is changed in 'DFT.m'
funcCubicDrift=@(a) sign(a)*abs(a)^(1/3);  

%%
hmaxListLength=length(hmaxList);  

% start counting the running time of entire script
% will be displayed when the simulation is done
tStart = tic;  

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
        Xref=Xzero;
        for k=1:parms.Nref
            % generating Wiener increment
            Winc=randn*parms.sqrtdtref;  
            
            % calculating reference solution by stepsize dtref
            Xref=schemePMil(parms, Xref, Winc, parms.dtref);
            
            % store all Wiener increment on each MC realisation
            % for the use of other methods
            WKEEP(i,k)=Winc; 
        end
        % store the terminal value for each MC realisaton 
        % for the use of other methods
        XrefKeep(i)=Xref;
        
        %--------------------   Adaptive Milstein  -------------------------%
        XAMil = Xzero;  
        % calculating AMil with detailed explaination in 'schemeAMil.m'
        [XAMil, j, hfix]=schemeAMil(parms, XAMil, WKEEP, i, j, hfix);
        % adding sqaured error together for the sample mean later
        XerrAMil(m)=XerrAMil(m)+abs(XAMil-XrefKeep(i))^2;
    end
    
    % Calculating hmean which will be used as
    %          the fixed stepsize for all fixed-step methods
    % 'hfix' is now the sum of adaptive steps from all MC realisations
    %          except the last step
    % 'j' is now the sum of number of adaptive step
    %          taken from all MC realisations
    hmean(m)=hfix/j;
    
    % hmean(m) now might not be a multiple of 'parms.dtref'
    % so we calculate the floor multiple
    n=floor(hmean(m)/parms.dtref);
    
    % the fixed stepsize we actually use is
    parms.dtuse=n*parms.dtref; 
    
    % number of the fixed steps taken up to 'parms.T'
    N=floor(parms.T/parms.dtuse);
    
    % calulate the length of last step
    % which should be shorter than 'parms.dtuse'
    Nlast=rem(parms.Nref,n); 
    parms.dtlast=Nlast*parms.dtref;
    
    
    % For the same 'parms.hmax'
    % run through all MC realisations for all fixed-step methods
    for i=1:M 
        %-------------------   fixed-step methods  ------------------------%
        W=0; 
        % assign initial values of all methods to 'Xzero'
        [XTMil, XPMil, XSSBM, XTSRK]=deal(Xzero);
        for k=1:N
            % pull Wiener increments from reference solution
            W=sum(WKEEP(i,(k-1)*n+1:k*n));
            
            % run all methods 
            XTMil=schemeTMil(parms, XTMil, W, 'fixed');
            XPMil=schemePMil(parms, XPMil, W, parms.dtuse);
            XSSBM=schemeSSBM(parms, XSSBM, W, 'fixed', funcCubicDrift);
            XTSRK=schemeTSRK(parms, XTSRK, W, 'fixed');
        end
        
        % run the last step to reach 'parms.T'
        if parms.dtlast ~= 0
            W=sum(WKEEP(i,N*n+1:end));
            XTMil=schemeTMil(parms, XTMil, W, 'last');
            XPMil=schemePMil(parms, XPMil, W, parms.dtlast);
            XSSBM=schemeSSBM(parms, XSSBM, W, 'last', funcCubicDrift);
            XTSRK=schemeTSRK(parms, XTSRK, W, 'last');
        end
        
        % adding sqaured error together for the sample mean later
        XerrTMil(m)=XerrTMil(m)+abs(XTMil-XrefKeep(i))^2;
        XerrPMil(m)=XerrPMil(m)+abs(XPMil-XrefKeep(i))^2;
        XerrSSBM(m)=XerrSSBM(m)+abs(XSSBM-XrefKeep(i))^2;
        XerrTSRK(m)=XerrTSRK(m)+abs(XTSRK-XrefKeep(i))^2;
        
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
[TimeAMil, TimeTMil, TimePMil, TimeSSBM, TimeTSRK ]=deal(zeros(1,hmaxListLength));
for m=1:hmaxListLength  % run through the list of hmax
    parms.hmax=hmaxList(m);
    parms.hmin=parms.hmax/parms.rho;
    parms.Tol=parms.hmin/pow2(3);
    for i=1:M  % run through all the MC realisatons 
        %------------------   Adaptive Milstein  ---------------------%
        tic;  % start counting CPU time for AMil
        XAMil = Xzero;   
        XAMil=schemeAMilCPU(parms, XAMil);
        % end counting CPU time for AMil and store  
        TimeAMil(m)=TimeAMil(m)+toc;  

        %------------------   Tamed Milstein  -----------------------%
        tic; % start counting CPU time for TMil
        XTMil=Xzero; 
        % Use 'hmean(m)' as the fixed stepsize for all fixed-step methods
        %               to make it fair to compare to AMil
        % It will have 'N' steps for size 'hmean(m)'   
        %               and one last step for size 'parms.dtlast'
        parms.dtuse=hmean(m); 
        N=floor(parms.T/parms.dtuse); 
        parms.dtlast=parms.T-parms.dtuse*N;
        for k=1:N
            W=randn*sqrt(parms.dtuse);
            XTMil=schemeTMil(parms, XTMil, W, 'fixed');
        end   
        if parms.dtlast ~= 0
            W=randn*sqrt(parms.dtlast);
            XTMil=schemeTMil(parms, XTMil, W, 'last');
        end
        % end counting CPU time for TMil and store
        TimeTMil(m)=TimeTMil(m)+toc;
        
        
        %----------------   Projected Milstein  ----------------------%
        tic;
        XPMil=Xzero; 
        parms.dtuse=hmean(m); 
        N=floor(parms.T/parms.dtuse); 
        parms.dtlast=parms.T-parms.dtuse*N;
        for k=1:N
            W=randn*sqrt(parms.dtuse);
            XPMil=schemePMil(parms, XPMil, W, parms.dtuse);
        end 
        if parms.dtlast ~= 0
            W=randn*sqrt(parms.dtlast);
            XPMil=schemePMil(parms, XPMil, W, parms.dtlast);
        end
        TimePMil(m)=TimePMil(m)+toc;       
        
        
        %-------------------   SSBM  ------------------------------%
        tic;
        XSSBM=Xzero; 
        parms.dtuse=hmean(m); 
        N=floor(parms.T/parms.dtuse); 
        parms.dtlast=parms.T-parms.dtuse*N;
        for k=1:N
            W=randn*sqrt(parms.dtuse);
            XSSBM=schemeSSBM(parms, XSSBM, W, 'fixed', funcCubicDrift);
        end   
        if parms.dtlast ~= 0
            W=randn*sqrt(parms.dtlast);
            XSSBM=schemeSSBM(parms, XSSBM, W, 'last', funcCubicDrift);
        end
        TimeSSBM(m)=TimeSSBM(m)+toc;  
        
        
        %------------------   TSRK -------------------------------%
        tic;
        XTSRK=Xzero; 
        parms.dtuse=hmean(m); 
        N=floor(parms.T/parms.dtuse); 
        parms.dtlast=parms.T-parms.dtuse*N;
        sqrtdtuse=sqrt(parms.dtuse);
        for k=1:N
            W=randn*sqrtdtuse;
            XTSRK=schemeTSRK(parms, XTSRK, W, 'fixed');
        end   
        if parms.dtlast ~= 0      
            W=randn*sqrt(parms.dtlast);
            XTSRK=schemeTSRK(parms, XTSRK, W, 'fixed');
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
figure(1)  % convergence plot
loglog(TimeAMilMean,RMSerrAMil,'k-o','LineWidth',3,'MarkerSize',17), hold on
loglog(TimePMilMean,RMSerrPMil,'b:+','LineWidth',3,'MarkerSize',17), hold on
loglog(TimeSSBMMean,RMSerrSSBM,':x','LineWidth',3,'MarkerSize',17,'Color',[0.9290, 0.6940, 0.1250]), hold on
loglog(TimeTMilMean,RMSerrTMil,'m--*','LineWidth',3,'MarkerSize',17), hold on
loglog(TimeTSRKMean,RMSerrTSRK,'r-.s','LineWidth',3,'MarkerSize',17), hold off
set(gca,'fontsize',14)
xlabel('CPU time','FontSize',14)
ylabel('RMS error','FontSize',14)
legend({'AMil','PMil','SSBM','TMil','TSRK1'},'Location','northeast','FontSize',12)
title('(d)')
axis tight
grid


figure(2)  % efficiency plot
loglog(hmean,RMSerrAMil,'k-o','LineWidth',3,'MarkerSize',17), hold on
loglog(hmean,RMSerrPMil,'b:+', 'LineWidth',3,'MarkerSize',17), hold on
loglog(hmean,RMSerrSSBM,':x','LineWidth',3,'MarkerSize',17,'Color',[0.9290, 0.6940, 0.1250]), hold on
loglog(hmean,RMSerrTMil,'m--*','LineWidth',3,'MarkerSize',17), hold on
loglog(hmean,RMSerrTSRK,'r-.s','LineWidth',3,'MarkerSize',17), hold on
loglog(hmean,hmean*exp(-2),'g-.','LineWidth',4), hold off
set(gca,'fontsize',14)
xlabel('hmean','FontSize',14)
ylabel('RMS error','FontSize',16)
legend({'AMil','PMil','SSBM','TMil','TSRK1','Rate 1'},'Location','northwest','FontSize',12)
title('(c)')
axis tight
grid

