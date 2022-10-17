function outcome=schemeTSRK(parms, StartValue, Winc, dtType)
% fixed-step Tamed Stochastic Runge-Kutta of order 1.0

if strcmp('fixed',dtType)
    dt=parms.dtuse;
else  % last step
    dt=parms.dtlast;
end

sqrtdt=sqrt(dt);
dft=DFT(parms.nonlinearity, StartValue);
DiffInput=DIFF(parms.sigma, StartValue, parms.ModelNoise);
H1=StartValue-DiffInput*Winc^2/(2*sqrtdt)...
    +DiffInput*sqrtdt/2;
H2=StartValue+DiffInput*Winc^2/(2*sqrtdt)...
    -DiffInput*sqrtdt/2;
Diff_H1=DIFF(parms.sigma, H1, parms.ModelNoise);
Diff_H2=DIFF(parms.sigma, H2, parms.ModelNoise);
outcome=StartValue+dft*dt/(1+dt*abs(dft))...
    +DiffInput*Winc/(1+dt*abs(DiffInput))...
    -0.5*Diff_H1*sqrtdt/(1+dt*abs(Diff_H1))...
    +0.5*Diff_H2*sqrtdt/(1+dt*abs(Diff_H2));


