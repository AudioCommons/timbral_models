function [ eta ] = adpt_thresholdmu( f )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    Nf = length(f);
    beta = [13 9 7 4 2.95 2.05 1.3 1 0.7 0.65 0.55 0.55 0.5 0.5 0.5 0.5 0.45 0.45 0.4 0.35 0.3 0.2 0.15 0.18 0.23 0.26 0.28 0.5 2.3 14];
    alpha = 2^(1/3);
    fmin = 20;
    I = (log(f/fmin)/log(alpha))+1;
    eta = zeros(size(f));
    for n = 1:Nf
        if (I(n) <= 1)
            eta(n) = beta(1) + (beta(1)-beta(2))*(1-I(n));
            continue;
        end
        if (I(n) >= 30)
            eta(n) = beta(30) + (beta(30)-beta(29))*(I(n)-30);
            continue;
        end
        N1 = floor(I(n));
        N2 = N1+1;
        eta(n) = beta(N1) + (I(n)-N1)*(beta(N2)-beta(N1));
    end


end

