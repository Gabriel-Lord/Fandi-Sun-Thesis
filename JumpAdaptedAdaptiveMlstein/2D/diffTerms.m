function outcome=diffTerms(W,Diff,Ddiff1,Ddiff2,dt,parms)
% to generate the terms from diffusion by Milstein method
% % INPUTS:
% % 'W' -- (2-by-1 vector)  Wiener increment for this step  
% % 'Diff' -- (2-by-2 matrix)  diffusion matrix
% % 'Ddiff1' 'Ddiff2' -- (2-by-2 matrix)  derivative of diffusion matrix 
% %           with respect to input vector
% % 'dt' -- (scalar)  step size  
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts  
% % OUTPUTS:
% % 'outcome' -- (2-by-1 vector) contains all terms from diffusion
if strcmp('diagonal',parms.ModelNoise)
    outcome=Diff*W+0.5*( Ddiff1*Diff(:,1)*(W(1)^2-dt)...
    +Ddiff2*Diff(:,2)*(W(2)^2-dt) );
elseif strcmp('commutative',parms.ModelNoise)
    outcome=Diff*W+0.5*( Ddiff1*Diff(:,1)*(W(1)^2-dt)...
    +Ddiff2*Diff(:,2)*(W(2)^2-dt) )...
    +0.5*( Ddiff1*Diff(:,2)+Ddiff2*Diff(:,1) )*W(1)*W(2);
else  % non-commutative
    outcome=Diff*W+0.5*( Ddiff1*Diff(:,1)*(W(1)^2-dt)...
    +Ddiff2*Diff(:,2)*(W(2)^2-dt) )...
    +0.5*( Ddiff1*Diff(:,2)+Ddiff2*Diff(:,1) )*W(1)*W(2)...
    +0.5*( Ddiff1*Diff(:,2)-Ddiff2*Diff(:,1) )*parms.LA;
end