function [ Sx ] = RAA_sigmoid_scaling( ax, bx, Px )
%RAA_sigmoid_scaling calculates scaled parameter values of RAA output 
%   based on sigmoid function

Sx = 1 / (1+exp(-ax*(Px-bx)));

end

