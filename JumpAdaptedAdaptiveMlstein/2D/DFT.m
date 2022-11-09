function outcome=DFT(nonlinearity,u)
% drift function
% Notice: if change the drift here, 
%       the TSRK method in 'main2D.m'
%               will need to be updated accordingly 
%       and the paramter parms.alpha 
%               might need to be updated  
% % INPUTS: 
% % 'nonlinearity' -- (scalar) parameter in drift
% % 'u' -- (2-by-1 vector)  input value vector 
% % OUTPUTS:
% % 'outcome' -- (2-by-1 vector) of drift
outcome= [u(2); u(1)]-nonlinearity*[u(1); u(2)].^3;