function [LA,Y2,Y4]=funcLevy(K,Y2,Y4,dW)
% main function in approximating Levy area
% % INPUTS
% % 'K' -- (scalar) number of mini steps 
% % 'Y2' 'Y4' -- (scalar) markers in Levy approximation
% % 'dW' -- (2-by-K vector) vector of all Wiener increments
% % OUTPUTS
% % 'LA' -- (scalar) Levy area approximation
% % 'Y2' 'Y4' -- (scalar) updated markers for later Levy approximation
Y1=0; Y3=0;
for a=1:K
    dWla=dW(:,a);
    Y1=Y1+Y2*dWla(1);
    Y2=Y2+dWla(2);
    Y3=Y3+Y4*dWla(2);
    Y4=Y4+dWla(1);
end
LA=Y1-Y3;