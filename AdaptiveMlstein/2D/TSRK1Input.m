function output=TSRK1Input(i,j,input,diff1,diff2,W, sqrtdt)
% function embeded in the 'schemeTSRK2D.m'
% % INPUTS:
% % 'i' -- (scalar) index value of a 2-by-2 matrix
% % 'j' -- (scalar) index value of a 2-by-2 matrix
% % 'input' -- (2-by-1 vector) initial vector value of the process for this step
% % 'W' -- (2-by-1 vector)  Wiener increment from the reference solution for this step  
% % 'diff1' 'diff2' -- (2-by-1 vector) first and second column vector in the diffusion matrix
% % 'sqrtdt' -- (scalar) square root value of the step size
% % OUTPUTS:
% % 'outcome' -- (2-by-1 vector)   
if i==1
    if j==1
        output=input-(diff1*W(1)+diff1*W(2))*W(1)/(2*sqrtdt)...
            +diff1*sqrtdt/2;
    else
        output=input-(diff2*W(1)+diff2*W(2))*W(1)/(2*sqrtdt)...
            +diff2*sqrtdt/2;
    end
else
    if j==1
        output=input+(diff1*W(1)+diff1*W(2))*W(1)/(2*sqrtdt)...
            -diff1*sqrtdt/2;
    else
        output=input+(diff2*W(1)+diff2*W(2))*W(1)/(2*sqrtdt)...
            -diff2*sqrtdt/2;
    end
end
