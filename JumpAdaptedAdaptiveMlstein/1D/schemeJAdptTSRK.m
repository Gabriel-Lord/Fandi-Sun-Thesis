function outcome=schemeJAdptTSRK(parms, StartValue, Winc, dt, JumpSize)
% fixed-step jump-adapted Tamed Stochastic Runge-Kutta of order 1.0
% % INPUTS:
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts
% % 'StartValue' -- (scalar) initial value of the process for this step
% % 'Winc' -- (scalar)  Wiener increment from the reference solution for this step  
% % 'dt' -- (scalar) step size
% % 'JumpSize' -- (scalar)  jump size occurred at the end of this step
% % OUTPUTS:
% % 'outcome' -- (scalar) terminal value of the process for this step

 
sqrtdt=sqrt(dt);
dft=DFT(parms.nonlinearity, StartValue);
DiffInput=DIFF(parms.sigma, StartValue, parms.ModelNoise);
H1=StartValue-DiffInput*Winc^2/(2*sqrtdt)...
    +DiffInput*sqrtdt/2;
H2=StartValue+DiffInput*Winc^2/(2*sqrtdt)...
    -DiffInput*sqrtdt/2;
Diff_H1=DIFF(parms.sigma, H1, parms.ModelNoise);
Diff_H2=DIFF(parms.sigma, H2, parms.ModelNoise);
outcomeMid=StartValue+dft*dt/(1+dt*abs(dft))...
    +DiffInput*Winc/(1+dt*abs(DiffInput))...
    -0.5*Diff_H1*sqrtdt/(1+dt*abs(Diff_H1))...
    +0.5*Diff_H2*sqrtdt/(1+dt*abs(Diff_H2));
outcome=outcomeMid+AMPLTD(JumpSize,outcomeMid);


