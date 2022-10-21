function outcome=AMilCoreH(input,parms)
% The core function in the adaptive time-stepping 
%           strategy for adaptive Milstein method
% Outcome needs to be bounded by
%           hmin and hmax outside the function
% % INPUTS:
% % 'input' -- (2-by-1 vector) initial value vector  
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts
% % OUTPUTS:
% % 'outcome' -- (scalar)  of adaptive step h
outcome=parms.hmax/norm(input)^(1/parms.kappa);