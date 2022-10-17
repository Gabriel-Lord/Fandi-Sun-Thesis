function [XAMil, j, hfix]=schemeAMil(parms, XAMil, WKEEP, i, j, hfix)
% % To generate the terminal value of the proposed Adaptive Milstein method
% %           based on a reference solution
% % For fixed-step methods:
% % 'j' -- number of adaptive steps taken on all MC realisations
% % 'hfix' -- sum length of adaptive steps taken on all MC realisations

% Since each step 'h' is adaptive and recalculated after each mesh grid 
%           we do not know how many reference Wiener increaments
%           they will take for each step
% Therefore, 'Marker' below acts like a bookmark 
%           which will record the starting position of each adaptive step
%           on the corresponding reference line 
% 'Marker' is updated after each adaptive step is done
Marker=0;   

% we add each adaptive step to 't' 
% until it reaches terminal time 'parms.T'
t=0;
while abs(t-parms.T)>parms.Tol && t<=parms.T
    
    % core function of choosing 'h'
    HcoreFunc=AMilCoreH(XAMil,parms);
    % bounding the function value by hmin and hmax
    h=max(parms.hmin, min(parms.hmax, HcoreFunc));  
    
    % 'h' might not be a multiple of single reference step
    % so we take the closest floor
    Ninc_adpt=floor(h/parms.dtref);
    if Marker+Ninc_adpt > parms.Nref  % indicates last step
        Ninc_adpt= parms.Nref-Marker;  % take whatever is left
        huse=Ninc_adpt*parms.dtref;
        WincAMil = sum(WKEEP(i,Marker+1:end),2);
        Marker=parms.Nref; % update 'Marker' to the end
    else  % indicates adaptive step
        huse=Ninc_adpt*parms.dtref; % 'h' that we actually use
        WincAMil = sum(WKEEP(i,Marker+1:Marker+Ninc_adpt),2);
        Marker=Marker+Ninc_adpt; % update 'Marker'
        %  For all MC realisations except the last step, record
        hfix=hfix+huse;  % adaptive stepsizes 
        j=j+1; % number of adaptive steps
    end
    
    if huse <= parms.hmin % backstop case
        % we use PMil with timestep 'huse' for the backstop
        XAMil=schemePMil(parms, XAMil, WincAMil, huse);
    else
        % explicit Milstein
        [dft,Diff,Ddiff]=coefficients(XAMil, parms);
        XAMil=XAMil+dft*huse+diffTerms(WincAMil,Diff,Ddiff,huse);
    end
    t=t+huse; % push 't'
end