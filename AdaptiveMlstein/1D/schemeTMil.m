function outcome=schemeTMil(parms, StartValue, Winc, dtType)
%  fixed-step Tamed Milstein method

if strcmp('fixed',dtType)
    dt=parms.dtuse;
else  % last step
    dt=parms.dtlast;
end

scale=1+dt*abs(StartValue)^4;
[dft,Diff,Ddiff]=coefficients(StartValue, parms);
outcome=StartValue+(dft*dt+diffTerms(Winc,Diff,Ddiff,dt))/scale;
