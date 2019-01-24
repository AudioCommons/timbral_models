function [out, signal1, signal2]=stimulusCustomBraasch(Fs,signal1,signal2,isi,nn,lag_level)
%
%[out, signal1, signal2] = STIMULUSCUSTOMBRAASCH(Fs, signal1, siangl2, isi, nn, lag_level)
% creates a mixed binaural lead-lag signal 
% with the interstimulus-interval isi and with the lag level lag_level.
%
% Parameters:
%     Fs        : sampling frequency
%     signal1/2 : original signal to be synthesised as lead-lag pair
%     isi       : ISI in milliseconds
%     nn        : 0-specified operation, 1-switch off lead, 2-switch off lag
%     lag_level : Specific lag level

ISI=abs(isi);

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