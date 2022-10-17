function outcome=AMilCoreH(input,parms)
% The core function in the adaptive time-stepping 
%           strategy for adaptive Milstein method
% Outcome needs to be bounded by
%           hmin and hmax outside the function
outcome=parms.hmax/abs(input)^(1/parms.kappa);