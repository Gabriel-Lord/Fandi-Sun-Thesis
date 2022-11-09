function [dt, jumpIndct, JumpSize, JumpTimesTemp, zetaTemp]=...
    JumpSizeFinder(t, dt, JumpTimesTemp, zetaTemp)
% % To find the jump size if there is a jump in this step
% % INPUTS:
% % 't' -- (scalar) starting time of this step
% % 'dt' -- (scalar) step size
% % 'JumpTimesTemp' -- (1-by-n vector)  jump times in this MC realisation
% % 'zetaTemp' -- (1-by-n vector)  jump size in this MC realisation 
% % OUTPUTS:
% % 'dt' -- (scalar) the step size to use if there is a jump
% % 'jumpIndct' -- (scalar) indicator of whether there is a jump
% % 'JumpTimesTemp' -- (1-by-(n-1) vector)  
% %         updated jump times in this MC realisation  
% % 'zetaTemp' -- (m-by-(n-1) matrix)  
% %         updated jump size in this MC realisation 

% the following if statement works for both fixed step method
% and the adaptive method
% fixed-step methods only have 'JumpTimesTemp(1) == t+dt'
% adaptive method has 'JumpTimesTemp(1) <= t+dt'
if JumpTimesTemp(1) <= t+dt   
    % update 'dt' if 'JumpTimesTemp(1) <= t+dt'
    dt=JumpTimesTemp(1) - t;
    JumpSize=zetaTemp(:,1); % jump size
    % removing the first checked element 
    JumpTimesTemp=JumpTimesTemp(2:end);  
    zetaTemp=zetaTemp(:,2:end);  
    % indicator of whether this step has a jump
    jumpIndct=1;  
else  % if there is no jump in this step
    JumpSize=0;
    jumpIndct=0;
end
 