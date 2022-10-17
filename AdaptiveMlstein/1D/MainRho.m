% % To generate: 
% %     Figure 1: plot of the paths
% %                    of adaptive steps 'h' when 'rho' increases
% %     Figure 2: averaged probability of using 
% %                   backstop and 'hmin' when 'rho' increases
% % drift is nonlinear in 'DFT.m'
% % diffusion is in 'DIFF.m'
clc;clear all;close all;
Xzero=100;  % initial value
T=1;
Nref=pow2(20);
hmax=T/Nref; 
rholist=[2,4,6,8,10,12,14,16];  % list of rho
hminlist=hmax./rholist;

% number of Monte Carlo realisations 
% for 'the probability of using hmin'
M=10;

% list of model parameters called from 'ModelParameters.m'
parms=ModelParameters;

% the noise of the model with details in 'DIFF.m' and 'DDIFF.m'
% with the choices of  'additive'; 'mul1'; and 'mul2'.
parms.ModelNoise='mul2';  %  'additive'  %  'mul1'  % 'mul2'

%%  path of 'h'
% 'hfix' is to store all the adaptive steps for each rho
hfix=0;

% Run through the list of rho and record the paths of 'h'
% For detailed explaination of adaptive Milstein method
%           see 'schemeAMil.m'
for m=1:length(rholist) % go through rholist list
    rho=rholist(m); hmin=hmax/rho; Tol=hmin/pow2(3);   
    j=0;  XAMil = Xzero;  t=0;
    while abs(t-T)>Tol && t<=T
        j=j+1; 
        h=max(hmin,min(hmax,hmax/abs(XAMil))); % choice of h
        WincAMil = randn*sqrt(h); %  Wiener increment
        if h <= hmin    % last step -- PMil
            XAMil=schemePMil(parms, XAMil, WincAMil, h);
        else  % explicit Milstein 
            [dft,Diff,Ddiff]=coefficients(XAMil, parms);
            XAMil=XAMil+dft*h...
                +diffTerms(WincAMil,Diff,Ddiff,h);
        end   
        t=t+h;  
        hfix(m,j)=h;
    end   
end
%  generating path plot
hfix_non0_1=nonzeros(hfix(1,:));  % path of 'h' when rho = 2
hfix_non0_2=nonzeros(hfix(3,:));  % path of 'h' when rho = 8
figure(1)
x_1=1:length(hfix_non0_1)-1;
x_2=1:length(hfix_non0_2)-1;
plot(x_1,hfix_non0_1(1:end-1),'b--','LineWidth',3), hold on
plot(x_2,hfix_non0_2(1:end-1),'r-','LineWidth',3), hold on 
yline(hmax,'k--','  hmax','LineWidth',1,'FontSize',15,'LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom')
yline(hminlist(1),'k--','hmin when \rho = 2','LineWidth',1,'FontSize',15)
yline(hminlist(3),'k--','hmin when \rho = 6','LineWidth',1,'FontSize',15)
grid
xlabel('steps','FontSize',14)
ylabel(' h','FontSize',14)
title('\fontsize{15} (e)')
legend({'path of h when \rho =2 ','path of h when \rho=6'},'Location','best','FontSize',15,'LineWidth',1,'AutoUpdate','off')

%%  probability of using 'hmin'
for m=1:length(rholist) % go through rholist list
    rho=rholist(m); hmin=hmax/rho; Tol=hmin/pow2(3); 
    % to count the times when using 'hmin' for all MC realisations
    hminCount=0;  
    % the number of all adaptive steps for all MC realisations
    j=0;
    for i=1:M  % run through all MC realisations
        XAMil = Xzero;    
        t=0; 
        while abs(t-T)>Tol && t<=T
            j=j+1;
            h=max(hmin,min(hmax,hmax/abs(XAMil))); % choice of h
            WincAMil = randn*sqrt(h);         
            if h <= hmin   
                XAMil=schemePMil(parms, XAMil, WincAMil, h);
                % adding 1 whenever it uses backstop
                hminCount=hminCount+1;        
            else 
                [dft,Diff,Ddiff]=coefficients(XAMil, parms);
                XAMil=XAMil+dft*h...
                    +diffTerms(WincAMil,Diff,Ddiff,h);
            end
            t=t+h; 
        end 
    end
    % Probability of taking backstop method
    %       for this particular 'rho'
    pr_bk(m)=hminCount/M/j;
end
% generating plot
figure(2)
loglog(rholist,pr_bk,'b-o','LineWidth',3,'MarkerSize',10), hold on
loglog(rholist,rholist.^(-1)*exp(-3.5),'r--','LineWidth',3,'MarkerSize',10), hold off
set(gca,'fontsize',14)
xlabel('$\rho$','FontSize',25,'interpreter','latex')
legend({['Probability of'  newline 'using hmin'],'ref. slope -1'},'Location','southwest','FontSize',15,'LineWidth',1,'AutoUpdate','off')
xticks(rholist)
xticklabels({'2', '4', '6', '8', '10', '12', '14', '16'})
yticks(fliplr(pr_bk))
yticklabels(string(fliplr(exp(log(round(pr_bk,3))))))
axis tight
grid
title('\fontsize{15} (f)')