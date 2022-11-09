function x=multinewton( x, t, c, parms)
% Performs multidimensional Newton s method for the function defined in f
% starting with x and running NumIters times. 
% % INPUTS:
% % 'x' -- (2-by-1 vector) initial vector (main variable) 
% %                  in the nonlinear implicit equation
% % 't' -- (scalar) time step size
% % 'c' -- (2-by-1 vector) constant term in nonlinear 
% %                  implicit equation 
% % 'parms' -- (object) parameter list from 'ModelParameters.m' 
% %           with updates from main scripts
% % OUTPUTS:
% % 'x' -- (2-by-1 vector) solved root

[y,dy]=funcDrift(x, t, c, parms.nonlinearity);
for j=1:parms.NumIters
    s=dy\y; 
    x=x-s;
    [y,dy]=funcDrift(x,t,c,parms.nonlinearity);
end
 