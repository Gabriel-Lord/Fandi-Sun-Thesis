function [outcome, j, hfix]=schemeAMil2DLevy(parms, StartValue, WKEEP, i, j, hfix)
% % To generate the terminal value of the proposed Adaptive Milstein method
% %           based on a reference solution
% % INPUTS:
% % 'parms' -- (object) parameter list from 'ModelParameters.m' with updates from
% %           main scripts
% % 'StartValue' -- (2-by-1 vector) initial vector of the process
% % 'WKEEP' -- (2-by-n vector) list of Wiener process from the reference solution  
% % 'j' -- (scalar) number of adaptive steps taken on all MC realisations
% % 'hfix' -- (scalar) sum length of adaptive steps taken on all MC realisations
% % OUTPUTS:
% % 'outcome' -- (2-by-1 vector) terminal value of the process 
% % 'j' -- (scalar) updated value
% % 'hfix' -- (scalar) updated value

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

% Y2 and Y4 are needed for Levy area computation
Y2=0; Y4=0;
while abs(t-parms.T)>parms.Tol && t<=parms.T

    % core function of choosing 'h'
    HcoreFunc=AMilCoreH(StartValue,parms);
    % bounding the function value by hmin and hmax
    h=max(parms.hmin, min(parms.hmax, HcoreFunc));
    
    % Levy area setting
    delta=h^2;  % size of ministep
    % number of reference ministeps within one ministep
    miniStepNo=floor(delta/parms.minilength); 
    % resize ministep to make it a multiple 
    % of reference ministeps
    delta=miniStepNo*parms.minilength;  
    % number of ministep within one h
    K=floor(h/delta);  
    % resize h to make it a multiple of ministeps
    huse=K*delta;   
    % number of reference ministep in one h
    miniStepNoTotal=miniStepNo*K; 


    if Marker+miniStepNoTotal > parms.Nref*parms.Kref  % last step
        miniStepNoTotal= parms.Nref*parms.Kref-Marker;
        huse=miniStepNoTotal*parms.minilength;
        WincAMil = sum(WKEEP((2*i-1):2*i,Marker+1:end),2);
        Marker=parms.Nref;
        % we remove Levy area approximation in the last step
        parms.LA=0; 
    else  %   before last step
        % mini Wiener increments 
        % for mini step within one h
        dWla=zeros(parms.rank,K);
        for a=1:K
            dWla(:,a)=sum(WKEEP((2*i-1):2*i, ...
                Marker+(a-1)*miniStepNo+1:Marker+a*miniStepNo) ,2);
        end
        % approximating Levy area
        [parms.LA,Y2,Y4]=funcLevy(K,Y2,Y4,dWla);
        % Wiener increment to use 
        WincAMil = sum(dWla,2);
        % updating indexes
        Marker=Marker+miniStepNoTotal; 
        hfix=hfix+huse;
        j=j+1;
    end 

    if huse <= parms.hmin % backstop case
        % we use PMil with timestep 'huse' for the backstop
        StartValue=schemePMil2D(parms, StartValue, WincAMil, huse);
    else
        % explicit Milstein
        [dft,Diff,Ddiff1,Ddiff2]=coefficients(StartValue, parms);
        StartValue=StartValue+dft*huse+diffTerms(WincAMil,Diff,Ddiff1,Ddiff2,huse,parms);
    end
    t=t+huse; % push 't'
end
outcome=StartValue;