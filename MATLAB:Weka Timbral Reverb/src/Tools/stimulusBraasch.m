function [out, signal1, signal2]=stimulusBraasch(Fs,mode, len,f,bw,itd1,itd2,isi,at,dc,nn,lag_level)
%Y=STIMULUS(Fs,mode,len,f,bw,itd1,itd2,isi,at,dc)
%uses the function SIGNAL to process mixed stimuli
%with the interstimulus-interval isi.
%Parameters:
%
%Fs:    sampling frequency
%mode:  integer to select a waveform:
%len:   Signal length in milliseconds
%bw:    Bandwidth of the FFT bandpass filters
%f:     For periodic waves: Frequency in Hz
%       For Bandpass Noise: Fc of the bandpass filter
%itd1,
%itd2:  ITDs in milliseconds
%       (positive for right channel first - left channel last,
%       negative for left channel first - right channel last)
%       Indices 1 and 2 correspond to the incidence of the
%       the two signals.
%isi:   ISI in milliseconds
%at:    attack time (ms)
%dc:    decay time (ms)
%nn:    0: specified operation, 1: switch off S1, 2: switch off S2

ISI=abs(isi);
signal1=signalBraasch(Fs,mode,len,f,bw,itd1,1,at,dc); % ILD is the difference in level between the ipsilteral and contralateral ears
signal2=signalBraasch(Fs,mode,len,f,bw,itd2,1,at,dc); % ILD should be a negative number so that the contralateral ear is attenuated


%% switch off either signal if nn = 1 or 2 - this is used for the reference condition so there is no reflection 
if (nn==1)
    signal1=zeros(size(signal1));
end
if (nn==2)
    signal2=zeros(size(signal2));
end

%%
isiSamples=fix(Fs*ISI/1000); % convert ISI to samples
spc=zeros(isiSamples,2);
signal2=[spc;signal2]; % this is the lag

%% pad the shorter signal with zeros to make them the same length

signal1=[signal1;zeros(length(signal2)-length(signal1),2)];
    
% make lag louder to test Haas effect. Level is already converted from dB
% to amplitude 
signal2 = signal2 .* lag_level; 

out=signal1+signal2;

norm_constant = max([max(abs(out(:,1))) max(abs(out(:,2)))]);
out=out./norm_constant; % normalise the output stimulus signal

% now, for the Farimaa model, normalize the left and right signals by the
% same value as well
signal1 = signal1 ./ norm_constant;
signal2 = signal2 ./ norm_constant;


%% swap the order (left to right) of the signals if itd is negative
if isi<0
    out=([out(:,2) out(:,1)]);
end % of if