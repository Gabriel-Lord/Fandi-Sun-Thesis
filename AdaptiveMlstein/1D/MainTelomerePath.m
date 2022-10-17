% % To generate the paths plot from five methods:
% %             AMil; PMil; TMil; SSBM; and TSRK
% % Based on the stochastic telomere model
clc;clear all;close all;
Xzero=1000; % initial value
T=25;  % terminal time in days

% list of model parameters called from 'ModelParameters.m'
parms=ModelParameters;  
parms.kappa=1000;  

% step setting for base step
NAMil=pow2(23); 
dtBase=T/NAMil;
sqrtdtBase=sqrt(dtBase);

% step setting for AMil step
parms.hmax=pow2(3)*dtBase; 
rho=pow2(2); 
hmin=parms.hmax/rho;
Tol=hmin/pow2(3);

% step setting for fixed-step method
Nuse=pow2(20);
dtuse=T/Nuse;
StepMulti=NAMil/Nuse;

% parameters for telomere model
a = 0.73*10^(-6);
c = 5.7;

% drift function
DFT=@(u) -c-0.5*a*u^2;  

% diffusion function
DIFFCons = sqrt(a/3);   % 369.6846
DIFF=@(u) DIFFCons*realsqrt(u)^3;
DDIFFCons = DIFFCons*3/2;  % 554.5268
DDIFF = @(u) DDIFFCons*realsqrt(u); 



% parameter for PMil
alpha=0.5;

% setting for SSBM
aa=@(h) a*h/2;
bb=1;
cc=@(h,x) h*c-x;

%%
tStart = tic;

WKEEP=randn(1,NAMil)*sqrtdtBase;

%----------------   TMil   ------------------------%
XTMil(1)=Xzero;
for k=1:Nuse
    W=sum(WKEEP((k-1)*StepMulti+1:k*StepMulti));
    input=XTMil(k);
    dft=DFT(input);
    diff=DIFF(input);
    Ddiff=DDIFF(input);
    scale=1+dtuse*abs(input)^(2);
    XTMil(k+1)=input+(dft*dtuse...
        +diffTerms(W,diff,Ddiff,dtuse))/scale;
end

%----------------   adpt  ------------------------%
XAMil(1) = Xzero;   Marker=0;
t=0; k=0;
while abs(t-T)>Tol && t<=T
    k=k+1; 
    HcoreFunc=AMilCoreH(XAMil(k),parms);
    h=max(hmin,min(parms.hmax,HcoreFunc)); % choice of h
    Ninc=floor(h/dtBase); % # dtBase in each LA mini step
    
    if Marker+Ninc > NAMil  % last step
        Ninc= NAMil-Marker;
        huse=Ninc*dtBase;
        W = sum(WKEEP(Marker+1:end));
        Marker=NAMil;
    else  
        huse=Ninc*dtBase;
        W = sum(WKEEP(Marker+1:Marker+Ninc));
        Marker=Marker+Ninc; 
    end
    
    if huse <= hmin  % last step
        input=(-bb+sqrt(bb^2-4*aa(huse)*cc(huse,XAMil(k))))/(2*aa(huse));
        diff=DIFF(input);
        Ddiff=DDIFF(input);
        XAMil(k+1)=input+diffTerms(W,diff,Ddiff,huse);
    else
        input=XAMil(k);
        dft=DFT(input);
        diff=DIFF(input);
        Ddiff=DDIFF(input);
        XAMil(k+1)=input+dft*huse...
            +diffTerms(W,diff,Ddiff,huse);
    end
    t=t+huse;
end

%----------------   PMil  ------------------------%
XPMil(1)=Xzero;
for k=1:Nuse
    W=sum(WKEEP((k-1)*StepMulti+1:k*StepMulti));
    input= min(1 , dtuse^(-alpha)*abs(XPMil(k))^(-1)) *XPMil(k);
    dft=DFT(input);
    diff=DIFF(input);
    Ddiff=DDIFF(input);
    XPMil(k+1)=input+dft*dtuse+diffTerms(W,diff,Ddiff,dtuse);
end

%----------------   SSBM   ------------------------%
XSSBM(1)=Xzero; W=0;
for k=1:Nuse
    W=sum(WKEEP((k-1)*StepMulti+1:k*StepMulti));
    input = (-bb+sqrt(bb^2-4*aa(dtuse)*cc(dtuse,XSSBM(k))))/(2*aa(dtuse));
    diff=DIFF(input);
    Ddiff=DDIFF(input);
    XSSBM(k+1)=input+diffTerms(W,diff,Ddiff,dtuse);
end

%----------------   TSRK1 ------------------------%
XTSRK(1)=Xzero; W=0;  sqrtdtuse=sqrt(dtuse);
for k=1:Nuse
    W=sum(WKEEP((k-1)*StepMulti+1:k*StepMulti));
    input=XTSRK(k);
    diff=DIFF(input);
    H1=input-diff*W^2/(2*sqrtdtuse)...
        +diff*sqrtdtuse/2;
    H2=input+diff*W^2/(2*sqrtdtuse)...
        -diff*sqrtdtuse/2;
    XTSRK(k+1)=input+DFT(input)*dtuse/(1+dtuse*abs(DFT(input)))...
        +diff*W/(1+dtuse*abs(diff))...
        -0.5*DIFF(H1)*sqrtdtuse/(1+dtuse*abs(DIFF(H1)))...
        +0.5*DIFF(H2)*sqrtdtuse/(1+dtuse*abs(DIFF(H2)));
end


XAMilMean=XAMil;
XTMilMean=XTMil;
XPMilMean=XPMil;
XSSBMMean=XSSBM;
XTSRKMean=XTSRK;

Tinc=T/Nuse;
Tx=0:Tinc:T;



tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

%%
Tx_adpt=0:T/length(XAMilMean):T;

figure(1)
plot(Tx_adpt(1:end-1),XAMilMean,'k-o','LineWidth',.5,'MarkerSize',8), hold on
plot(Tx,XPMilMean,'b:+','LineWidth',1.5,'MarkerSize',8), hold on
plot(Tx,XSSBMMean,':x','LineWidth',1.5,'Color',[0.9290, 0.6940, 0.1250],'MarkerSize',8), hold on
plot(Tx,XTMilMean,'m--*','LineWidth',1.5,'MarkerSize',8), hold on
plot(Tx,XTSRKMean,'r-.s','LineWidth',1.5,'MarkerSize',8), hold on
legend({'AMil','PMil','SSBM','TMil','TSRK1'},'Location','east','FontSize',10,'AutoUpdate','off')
grid
title('\fontsize{15} (a)')
xlabel('\fontsize{15} days')
ylabel('\fontsize{15} TL in base pairs')

plot([15 16.8], [1000  790.7], 'k','LineWidth',1)
plot([15 10.3], [1000  790.7], 'k','LineWidth',1)
gca1=gca;

start=floor(7/Tinc);
ends=floor(7.0003/Tinc);

style={'b:+',1, 12};
zoomPlot(Tx.',XPMilMean,style,gca1, [start*Tinc  ends*Tinc], [0.2 0.35 0.2 0.3],[0 0 0 0],[ 3 4]);
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
xlabel('PMil','FontSize',11)

start=floor(15/Tinc);
ends=floor(15.0003/Tinc);

xbounds=[start*Tinc  ends*Tinc];
style={'m--*',1, 12};
zoomPlot(Tx.',XTMilMean,style,gca1, xbounds, [0.45 0.35 0.2 0.3],[0 0 0 0],[1 2]);
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
xlabel('TMil','FontSize',11)

%%
figure(2)
plot(Tx_adpt(1:end-1),XAMilMean,'k-o','LineWidth',2,'MarkerSize',2), hold on
plot(Tx,XSSBMMean,':x','LineWidth',2,'Color',[0.9290, 0.6940, 0.1250],'MarkerSize',2), hold on
plot(Tx,XTSRKMean,'r-.s','LineWidth',2,'MarkerSize',2), hold off
grid
title('\fontsize{15} (b)')
xlabel('\fontsize{15} days')
ylabel('\fontsize{15} TL in base pairs')


start=floor(24.9997/Tinc);
ends=floor(25/Tinc);
style={'k-o',2, 12};
[~, z] = zoomPlot(Tx_adpt(2:end).',XAMilMean,style,gca, [start*Tinc  ends*Tinc], [0.22 0.22,0.3 0.2],[0 0 -0.1 0.35],[2 3]);
hold(z, 'on');
xinput=(start*Tinc):Tinc:(ends*Tinc);
plot(xinput,XSSBMMean(start:ends),':x' ,'Color',[0.9290, 0.6940, 0.1250],'LineWidth',style{2},'MarkerSize',style{3})
plot(xinput,XTSRKMean(start:ends),'r-.s' ,'LineWidth',style{2},'MarkerSize',style{3})
grid
legend({'AMil','SSBM','TSRK1'},'Location','southwest','FontSize',10,'LineWidth',1,'AutoUpdate','off')
