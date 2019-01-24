function y=de_conv(x,lag,amplitude,Fs);

% function y=de_conv(x,lag,amplitude,Fs);
%
% deconvolution filter 
% 
% J. Braasch (2013) J. Acoust. Soc. Am. 134(1): 420-435.
% For section references see paper
%
% (c) 2013, Rensselaer Polytechnic Institute
% contact: jonasbraasch@gmail.com
%
% INPUT PARAMETERS:
% x         = input signal
% lag       = delay in tabs
% amplitude = lead/lag amplitude ratio
% Fs        = sampling frequency (tested for 48 kHz)
% OUTPUT PARAMETERS:
% y         = deconvolved signal

maxlag=100; % max possilbe lag in ms
Fms=round(Fs./1000);
% slope to cut-off filter tail, see Eq. 9, p. 427
a=zeros(maxlag.*Fms,1);
a(1:6*Fms)=1;
a(6*Fms+1:20*Fms)=linspace(1,0,14*Fms).^4; 

% if lag is larger than maxlag 
if lag>maxlag.*Fms;
    lag=maxlag.*Fms; % set to maxlag
end % of if 

% Now, deconvolve signal, Eq. 3
att=amplitude; % redefinition to make equation below shorter
dur=length(x); % signal length/duration in tabs

y=zeros(dur+5*lag,1);
y(1:dur)=y(1:dur)+x;
for n=1:5 % filter order
    y(1+lag*n:lag*n+dur)=y(1+lag*n:lag*n+dur)+x.*(-att.*a(lag)).^n;
end 