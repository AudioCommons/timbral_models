function y=reconHW(x,Fs,f0)

% function y=reconHW(x,Fs,f0)
%
% reconstruction of full wave from halfwave rectified signal
% based on center frequency of auditory band
% see Eq. 6 and Fig. 8
% 
% J. Braasch (2013) J. Acoust. Soc. Am. 134(1): 420-435.
% For section references see paper
%
% (c) 2013, Rensselaer Polytechnic Institute
% contact: jonasbraasch@gmail.com
%
% INPUT PARAMETERS:
% x         = half-wave rectified auditory band-wide input signal
% Fs        = sampling frequency (tested for 48 kHz)
% f0        = center frequency of auditory band 
% OUTPUT PARAMETERS:
% y         = full-wave signal

y=x;

if f0>0
    T=round(0.5*Fs./f0);
    y(1:length(x)-T)=x(1:length(x)-T)-x(T+1:length(x));
end % of if 