function [B1,B2,A1,A2] = headphonefilter_Dorp2011(fs,order)
% function [B1,B2,A1,A2] = headphonefilter_Dorp2011(fs,order)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%       See also inline function in m20160714_fourth_report_MBBM_testing_vanDorp.m
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 11/07/2016
% Last update on: 11/07/2016 
% Last use on   : 11/07/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freqs=[1000 4000];

if nargin < 2
    order = 1; % butterworth: slope=20*order/oct
end

[B1,A1]=butter(order, freqs(1)/fs*2,'high');
[B2,A2]=butter(order, freqs(2)/fs*2,'low');

[H1,F]=freqz(B1,A1,[],fs);
[H2  ]=freqz(B2,A2,[],fs);

Ht_max = max(abs(H1.*H2));
B1 = sqrt(1/Ht_max)*B1;
B2 = sqrt(1/Ht_max)*B2;

if nargout == 0
    H1=freqz(B1,A1,[],fs);
    H2=freqz(B2,A2,[],fs);
    figure
    semilogx(F,20*log10([abs(H1) abs(H2)])); grid on
    title(sprintf('fs=%.0f [Hz]',fs))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
