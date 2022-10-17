function outcome=schemePMil(parms, StartValue, Winc, dt)
% fixed-step Projected Milstein method
input= min(1 , dt^(-parms.alpha)*abs(StartValue)^(-1)) *StartValue;
[dft,Diff,Ddiff]=coefficients(input, parms);
outcome=input+dft*dt+diffTerms(Winc,Diff,Ddiff,dt);
