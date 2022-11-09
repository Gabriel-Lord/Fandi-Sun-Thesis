%BrownianBridge
%
% given increment dW on Dt split into
% interval size Dt1 and Dt2 with increments
% dW1 = Dt1/Dt dW + sqrt(Dt1*Dt2/Dt) zeta
% dW2 = dW-dW1
function [dW1,dW2]=bbridge(dW,Dt,Dt1)
%   if nargin<5
%     [J,~]=size(W);
%   end
  [J,~]=size(dW);
  Dt2=Dt-Dt1;
  tmp=(Dt1/Dt);
  z=randn(J,1);
  dW1 = tmp*dW + sqrt(tmp*Dt2)*z;
  dW2 = dW-dW1; 
return