function outcome=DFT(nonlinearity,u)
% drift function
% Notice: if change the drift here, 
%       the TSRK method in 'main_1D.m'
%               will need to be updated accordingly 
%       and the paramter parms.alpha 
%               might need to be updated
outcome= u-nonlinearity*u^3; 