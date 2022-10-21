function outcome=schemeTSRK2D(parms, StartValue, Winc, dtType)
% fixed-step Tamed Stochastic Runge-Kutta of order 1.0 
% % INPUTS:
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts
% % 'StartValue' -- (2-by-1 vector) initial vector of the process for this step
% % 'Winc' -- (2-by-1 vector)  Wiener increment from the reference solution for this step  
% % 'dtType' -- (string) type of step, either 'fixed' or 'last'
% % OUTPUTS:
% % 'outcome' -- (2-by-1 vector) terminal value of the process for this step
if strcmp('fixed',dtType)
    dt=parms.dtuse;
else  % last step
    dt=parms.dtlast;
end
sqrtdt=sqrt(dt);
dft=DFT(parms.nonlinearity, StartValue);
DiffInput=DIFF(parms, StartValue, parms.ModelNoise);
Diff1=DiffInput(:,1);
Diff2=DiffInput(:,2);
H11=TSRK1Input(1,1,StartValue,Diff1,Diff2,Winc,sqrtdt);
H12=TSRK1Input(1,2,StartValue,Diff1,Diff2,Winc,sqrtdt);
H21=TSRK1Input(2,1,StartValue,Diff1,Diff2,Winc,sqrtdt);
H22=TSRK1Input(2,2,StartValue,Diff1,Diff2,Winc,sqrtdt);
DiffH11=DIFF(parms,H11,parms.ModelNoise);
DiffH12=DIFF(parms,H12,parms.ModelNoise);
DiffH21=DIFF(parms,H21,parms.ModelNoise);
DiffH22=DIFF(parms,H22,parms.ModelNoise); 
outcome=StartValue+dft*dt*TSRK1Scale(dt,dft)...
                +Diff1*Winc(1)*TSRK1Scale(dt,Diff1)...
                +Diff2*Winc(2)*TSRK1Scale(dt,Diff2)...
                -0.5*DiffH11*sqrtdt*TSRK1Scale(dt,DiffH11)...
                -0.5*DiffH12*sqrtdt*TSRK1Scale(dt,DiffH12)...
                +0.5*DiffH21*sqrtdt*TSRK1Scale(dt,DiffH21)...
                +0.5*DiffH22*sqrtdt*TSRK1Scale(dt,DiffH22);
        


