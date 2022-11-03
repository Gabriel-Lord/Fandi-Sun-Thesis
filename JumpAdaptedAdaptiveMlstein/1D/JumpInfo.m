function [JumpTimes, zeta]=JumpInfo(parms)
% % To generate the information of jumps including
% %           jump times and jump sizes
% % waiting time -- Exponential distribution
% % jump size -- Standard Normal distribution
% % INPUTS:
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main script
% % OUTPUTS:
% % 'JumpTimes' -- (M-by-n matrix) each row contains
% %           jump times in real line in one MC realisation
% % 'zeta' -- (M-by-n matrix) each row contains
% %           jump sizes occurred at their corresponding
% %           jump times in one MC realisation
 
zeta=zeros(parms.rank*parms.M,1); pi=zeros(parms.M,1);
% run through each MC realisation
for i=1:parms.M
    k=0; tau=0;
    while tau < parms.T
        k=k+1;
        % generate waiting time
        piTemp=exprnd(1/parms.lambda);
        if tau+piTemp > parms.T % true if waiting time steps over T
            if sum(pi(i,:)) == 0 % true if no jump has occurred
                piTempNew=exprnd(parms.lambda,1,5);
                piTemp=piTempNew(find(piTempNew<parms.T-tau,1));
            else  % true if some jumps have occurred
                break
            end
        end
        pi(i,k)=piTemp;
        % generate jump size using Standard Normal distribution
        zeta(parms.rank*i-(parms.rank-1):parms.rank*i, k)=randn(parms.rank,1);
        tau=tau+piTemp;
    end
end
pi(pi==0)=NaN;
% add up waiting times to generate jump times
JumpTimes=cumsum(pi,2);
JumpTimes(isnan(JumpTimes))=0;
if size(JumpTimes,1) ~= parms.M
    error('JumpTimes row number ~= M')
end
