% Main script for generating:
%     Figure 1. Poisson counting process
%          --- jump size 1 + uniform waiting time
%     Figure 2. compound Poisson process
%          --- jump size 'standard normal distribution' + uniform waiting time
%     Figure 2. compound Poisson + Wiener process
%          --- jump size 'standard normal distribution' + Wiener waiting time

clc;clear all; close all;
% terminal time
T=1;
% jump intensity 
lambda=5; 
% create vector of jump times and sizes
% for all 3 models 
i=0; tau=0; JumpTimes=0; zeta=0;
while tau<= T
    i=i+1;
    pi=exprnd(1/lambda);
    if JumpTimes(i)+pi>T
        pi=T-JumpTimes(i);
        JumpTimes(i+1)=JumpTimes(i)+pi;
        zeta(i+1)=randn;
        break
    end
    JumpTimes(i+1)=JumpTimes(i)+pi;
    zeta(i+1)=rand;
    tau=JumpTimes(i+1);
end

%%   Poisson Counting process
figure(1)  
length_jump=length(JumpTimes);
% plotting horizontal solid line for waiting times
for i=1:length_jump-1
    plot([JumpTimes(i), JumpTimes(i+1)],[i+1,i+1],'k-','LineWidth',3), hold on
end
% plotting vertical dashed line for jumps 
for i=1:length_jump-1
    plot([JumpTimes(i), JumpTimes(i)],[i,i+1],'k--','LineWidth',1), hold on
end
set(gca,'fontsize',14)
xlabel('time','FontSize',14)
ylabel('counting number','FontSize',14)
title('Poisson Counting process')
grid

%%  Compound Poisson process
figure(2)  
length_jump=length(JumpTimes);
% plotting horizontal solid line for waiting times
for i=1:length_jump-1
    plot([JumpTimes(i), JumpTimes(i+1)],[zeta(i+1),zeta(i+1)],'k-','LineWidth',3), hold on
end
% plotting vertical dashed line for jumps 
for i=1:length_jump-1
    plot([JumpTimes(i), JumpTimes(i)],[zeta(i),zeta(i+1)],'k--','LineWidth',1), hold on
end
set(gca,'fontsize',14)
xlabel('time','FontSize',14)
ylabel('process value','FontSize',14)
title('Compound Poisson process')
grid


%%  SDE driven by Compound Poisson
figure(3)  
length_jump=length(JumpTimes);
% plotting horizontal solid line
Nref=pow2(8);
dt=1/Nref;
% plot Wiener waiting time
for i=1:length_jump-1
    waiting_time=JumpTimes(i+1)-JumpTimes(i);
    stepno=round(waiting_time/dt);
    dtuse=waiting_time/stepno;
    sqrtdtuse=sqrt(dtuse);
    W=randn(1,stepno)*sqrtdtuse;
    plot([(JumpTimes(i)+dtuse):dtuse:JumpTimes(i+1)],zeta(i+1)+W,'k-','LineWidth',2), hold on
end
% plotting vertical dashed line for jumps
for i=1:length_jump-1
    plot([JumpTimes(i), JumpTimes(i)],[zeta(i),zeta(i+1)],'k--','LineWidth',1), hold on
end
set(gca,'fontsize',14)
xlabel('time','FontSize',14)
ylabel('process value','FontSize',14)
title('Compound Poisson + Wiener process') 
grid





