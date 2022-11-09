function [AllTimes,  timesInc]=allTimePoints(parms, WTimes, jumpTimes) 
% % To generate the superposition of uniform time steps and jump times 
% % INPUTS:
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main script
% % 'WTimes' -- (1-by-Nref vector)  vector of reference uniform steps
% % 'jumpTimes' -- (M-by-n matrix)  matrix of jump times for each MC realisations
% % OUTPUTS:
% % 'AllTimes' -- (M-by-(Nref+n) matrix) 
% %           superposition of uniform time steps and jump times 
% % 'timesInc' -- (M-by-(Nref+n) matrix) each row contains 
% %           time difference of the superposition of uniform time steps and jump times 

% creating empty matrix for storing data later
AllTimes=zeros(parms.M, parms.maxLength);
timesInc=zeros(parms.M, parms.maxLength); 

% run each MC realisations
for i=1:size(jumpTimes,1)
    % combine uniform steps and jump times
    combinedTimes=sort([WTimes, jumpTimes(i,:)]);
    % remove all repeated times 
    combinedTimes(diff(combinedTimes) == 0)=[]; 
    % store all time points
    AllTimes(i, 1:length(combinedTimes))=combinedTimes;
    % store all time increments of all time points, mostly are dtref
    % some are smaller due to jumps
    timesInc(i, 1:length(combinedTimes)-1)=diff(combinedTimes);
end
 