function output=TSRK1Scale(dtuse,input)
% scaling function embeded in the 'schemeTSRK2D.m'
% % INPUTS:
% % 'dtuse' -- (scalar) step size
% % 'input' -- (2-by-1 vector)  
% % OUTPUTS:
% % 'outcome' -- (scalar)  
    output=1/(1+dtuse*norm(input));