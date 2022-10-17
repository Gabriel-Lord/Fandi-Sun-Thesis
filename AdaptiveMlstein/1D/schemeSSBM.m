function outcome=schemeSSBM(parms, StartValue, Winc, dtType, funcCubicDrift)
%  fixed-step Split-Step Backward Milstein
if strcmp('fixed',dtType)  
    dt=parms.dtuse;
else  % last step
    dt=parms.dtlast;
end

% following is to solve the cubic equation for drift
p=(dt-1)/(dt*parms.nonlinearity);
q= StartValue/(dt*parms.nonlinearity);
s = sqrt(q^2/4 - p^3/27);
input =  funcCubicDrift(q/2 + s)+funcCubicDrift(q/2-s);
% substituting 'input' into the diffusion terms
Diff=DIFF(parms.sigma, input, parms.ModelNoise);
Ddiff=DDIFF(parms.sigma, input, parms.ModelNoise);
outcome=input+diffTerms(Winc,Diff,Ddiff,dt);

