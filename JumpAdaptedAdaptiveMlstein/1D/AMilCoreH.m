function outcome=AMilCoreH(input,parms)
% The core function in the adaptive time-stepping 
%           strategy for adaptive Milstein method
% Outcome needs to be bounded by
%           hmin and hmax outside the function
% % INPUTS:
% % 'input' -- (scaler) initial value 
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts
% % OUTPUTS:
% % 'outcome' -- (scalar)  of adaptive step h
outcome=parms.hmax/abs(input)^(1/parms.kappa);