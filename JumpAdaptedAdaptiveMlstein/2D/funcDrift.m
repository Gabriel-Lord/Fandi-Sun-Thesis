function [y, dy]=funcDrift(x, t, c, nonlinearity)  
% function for solving multidimentional Newton method
% in 'multinewton.m'
% % INPUTS:
% % 'x' -- (2-by-1 vector) initial vector (main variable) 
% %                  in the nonlinear implicit equation
% % 't' -- (scalar) time step size
% % 'c' -- (2-by-1 vector) constant term in nonlinear 
% %                  implicit equation 
% % 'nonlinearity' -- (scalar) coefficient
% % OUTPUTS:
% % 'y' -- (2-by-1 vector) nonlinear implicit equation 
% %             for drift from SSBM method
% % 'dy' -- (2-by-2 matrix) Jacobian matrix of 'y'
% %             which is the derivative of vector 'y' 
% %             with respect of vector 'x'  
y=zeros(2,1);  
dy=zeros(2,2);    
% vector 'y'
y(1)=-nonlinearity*t*x(1)^3-x(1)+t*x(2)+c(1);
y(2)=-nonlinearity*t*x(2)^3-x(2)+t*x(1)+c(2);
% elements of Jacobian matrix of 'y'
dy(1,1)=-3*nonlinearity*t*x(1)^2-1;   dy(1,2)=t;
dy(2,1)=t;                                          dy(2,2)=-3*nonlinearity*t*x(2)^2-1;



