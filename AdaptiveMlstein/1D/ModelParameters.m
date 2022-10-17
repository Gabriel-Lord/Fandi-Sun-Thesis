function parms=ModelParameters
parms.T=1;  % terminal time
parms.Nref=pow2(18); % number of reference steps
parms.rho=pow2(2);  % hmax/hmin
parms.dtref=parms.T/parms.Nref; % size of one reference step
parms.alpha=0.25; % parameter for PMil, for a cubic drift
parms.nonlinearity = 3;  % nonlinearity parameter for drift in 'DFT.m'
parms.sigma=0.2;  %  diffusion parameter in 'DIFF.m' and 'DDIFF.m'
% square root of reference stepsize, needed for reference Wiener increments
parms.sqrtdtref=sqrt(parms.dtref); 