function y=PeakRatio(ratio);

% function y=PeakRatio(ratio);
%
% determines relative lead/lag ratio from autorcorrelation function
% according to Eq. 2 and Appendix. Function cannot resolve which signal has
% the higher intensity lead or lag.
% 
% J. Braasch (2013) J. Acoust. Soc. Am. 134(1): 420-435.
% For section references see paper
%
% (c) 2013, Rensselaer Polytechnic Institute
% contact: jonasbraasch@gmail.com
%
% INPUT PARAMETERS:
% ratio = measured main/side peak ratio of auto-correlation function
% OUTPUT PARAMETERS:
% y     = relative lead/lag ratio between 0 and 1


if ratio>0.5
    ratio=0.5;
end % of if 
if ratio<0
    ratio=0;
end % of if 

m=1;
for n=0:0.01:1
   A=1;
   B(m)=n;
   C(m)=A^2+B(m)^2;
   D(m)=A*B(m);
   m=m+1;
end % of for 
index=find((D./C)>=ratio);
y=B(index(1));
