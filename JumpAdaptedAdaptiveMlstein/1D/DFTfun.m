function output=DFTfun(nonlinearity,input)
output=input-nonlinearity*input.^3;
return