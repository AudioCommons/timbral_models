function out=signalBraasch(Fs,mode,len,f,bw,id,LT,at,dc);

%Y=SIGNAL(Fs,mode,len,f,bw,itd,LT,at,dc,pw)
%generates a Signal for ITD/ILD-Experiments.
%Returns the Stereo-Signal-Vector Y.
%Parameters:
%
%Fs:    sampling frequency
%mode:  integer to select a waveform:
%           0 - Sine Wave
%           1 - Triangle Wave
%           2 - Bandpass Noise
%           3 - White Noise (Uniform Distribution)
%           4 - Two sine waves (1st & 3rd harmonics; f = f of 2nd harmonic)
%           5 - Three sine waves (1st, 2nd & 3rd harmonics; f = f of 2nd harmonic)
%           6 - peak train
%len:   Signal length in milliseconds

%   f:     For periodic waves: Frequency in Hz
%          For Bandpass Noise: Fc of the bandpass filter

%   bw:    Bandwidth of the FFT bandpass filters

%   id:    ITD in milliseconds  OR
%          ILD in dB
%       (positive for right channel first/loud - left channel last,
%       negative for left channel first/loud - right channel last)

%LT:   selects ITD or ILD :  0 for ILD, 1 for ITD
%at:    attack time (ms)
%dc:    decay time (ms)
%pw:    --- this has something to do with peak train

% ILD is for the difference in level between the ipsilateral and contralateral ears

at=at/1000;
dc=dc/1000;

%Calculating constants

%%  ITD
if(LT==1)
    lenSamples=fix(Fs*len/1000);
    itdSamples=round(Fs*id/1000);

    
    %% MODELS
    
    %% Sine Model
    if(mode==0)
        t=0:lenSamples-1;
       t=t/Fs;
       y=cos(2*pi*f*t);
       y=y';
    end

    %% Triangle Model
    if(mode==1)
    h1=(-1)*(1:(round(Fs/f/2)))./(1:(round(Fs/f/2))) + (2/round(Fs/f/2))*(1:round(Fs/f/2));
    h2=(1:(round(Fs/f/2)))./(1:(round(Fs/f/2))) - (2/round(Fs/f/2))*(1:round(Fs/f/2));
    h=[h1 h2];
    y(1:lenSamples)=0*(1:lenSamples)';
    m=1;
    while (m+length(h)-1)<=length(y)
     y(m:(m+length(h)-1))=h';
     m=m+length(h);
    end;
    if m<length(y) 
     data(m:length(y))=h(1:(1+length(y)-m))';
    end;
    y=y';
    end

    %% Band Pass Noise Model - this is the one used in precedence effect tests!
    if(mode==2)
        randn('seed',0);
        y=0*(1:lenSamples);
       u=f-round(bw/2); % find lower freq limit for bandwidth
      o=f+round(bw/2); % find upper freq limit for bandwidth
      u=round(lenSamples*u/Fs);
      o=round(lenSamples*o/Fs);
      y(u:o)=randn(1,(o-u+1))+i*randn(1,(o-u+1));
      y((lenSamples-o+1):(lenSamples-u+1))=conj(y(u:o));
      y=real(fft(y));   
      y=y';
    end

    %% Noise Model
    if(mode==3)
        rand('seed',0);
        y=rand(lenSamples,1);
    end

    
    %% 4 - Two sine waves (1st & 3rd harmonics; f = f of 2nd harmonic)
    if(mode==4)
       t=0:lenSamples-1;
       t=t/Fs;
       y=cos(pi*f*t)+cos(2*pi*f*1.5*t);
       y=y';
    end
    
    
    %% 5 - Three sine waves (1st, 2nd & 3rd harmonics; f = f of 2nd harmonic)
    if(mode==5)
       t=0:lenSamples-1;
       t=t/Fs;
       y=cos(pi*f*t)+cos(2*pi*f*t)+cos(2*pi*f*1.5*t);
       y=y';
    end

    
   %%  6 - peak train (click train) 5 ms noise separated by 100 ms
    if(mode==6)
        lenSamples=Fs; % 1 second signal length
        y=zeros(lenSamples,1);
        rand('seed',0);
        z=rand(1,fix(Fs*.005)); % 5 ms clicks 
        N=fix(100/1000 * Fs); % 100 ms gaps between clicks
       
        for n=1:N:lenSamples
            y(n:(n+length(z)-1))=z;
        end
    end

    
    %%
    %Scaling and removing Ey
    y=y-mean(y);
    y=y/max(abs(y));
    
    %%
    % Apply windowing function
    winlat=fix(Fs*at);
    winl=fix(Fs*dc);
    t=0:winl-1;
    t=t/((winl-1)/(.5*pi));
    decay=(cos(t).^2)';
    t2=0:winlat-1;
    t2=t2/((winlat-1)/(.5*pi));
    decay2=(cos(t2).^2)';
    attack=flipud(decay2); % flip upside down
    y(1:winlat)=y(1:winlat).*attack;
    y(length(y)-winl+1:length(y))=y(length(y)-winl+1:length(y)).*decay;
    
    %%
    % Generating the Stereo Signal Data
    dR=fix(abs(itdSamples/2))-fix(itdSamples/2); %  one of these will come out to be zero and the other will come out to be the ITD in samples
    dL=fix(abs(itdSamples/2))+fix(itdSamples/2); % this method just makes it work with pos or neg ITD

    yRight=zeros(abs(itdSamples)+length(y),1);
    yRight(dR+1:dR+length(y))=y; % depending on whether dR or dL is 0 or the ITD in samples, this puts  the zeros at the beginning of the signal
     % for the contra lateral ear or the zeros at the end of the signal for the ipsilateral ear

    yLeft=zeros(abs(itdSamples)+length(y),1);
    yLeft(dL+1:dL+length(y))=y;
    

        

    out=[yLeft yRight];
    %out=[yRight yLeft];
end





%% ILD  
if(LT==0)
    lenSamples=fix(Fs*len/1000);
    LeftLevel=10^(-(id)/20); % converts dB into amplitude
    RightLevel=10^((id)/20);
    
    %Sine Model
    if(mode==0)
        t=0:lenSamples-1;
       t=t/Fs;
       y=cos(2*pi*f*t);
       y=y';
    end

    %Triangle Model
    if(mode==1)
    h1=(-1)*(1:(round(Fs/f/2)))./(1:(round(Fs/f/2))) + (2/round(Fs/f/2))*(1:round(Fs/f/2));
    h2=(1:(round(Fs/f/2)))./(1:(round(Fs/f/2))) - (2/round(Fs/f/2))*(1:round(Fs/f/2));
    h=[h1 h2];
    y(1:lenSamples)=0*(1:lenSamples)';
    m=1;
    while (m+length(h)-1)<=length(y)
     y(m:(m+length(h)-1))=h';
     m=m+length(h);
    end;
    if m<length(y) 
     data(m:length(y))=h(1:(1+length(y)-m))';
    end;
    y=y';
    end

    %BP Noise Model
    if(mode==2)
        randn('seed',0);
        y=0*(1:lenSamples);
       u=f-round(bw/2);
      o=f+round(bw/2);
      u=round(lenSamples*u/Fs);
      o=round(lenSamples*o/Fs);
      y(u:o)=randn(1,(o-u+1))+i*randn(1,(o-u+1));
      y((lenSamples-o+1):(lenSamples-u+1))=conj(y(u:o));
      y=real(fft(y));   
      y=y';
    end

    %Noise Model
    if(mode==3)
        rand('seed',0);
        y=rand(lenSamples,1);
    end

    if(mode==4)
       t=0:lenSamples-1;
       t=t/Fs;
       y=cos(pi*f*t)+cos(2*pi*f*1.5*t);
       y=y';
    end
    
    if(mode==5)
       t=0:lenSamples-1;
       t=t/Fs;
       y=cos(pi*f*t)+cos(2*pi*f*t)+cos(2*pi*f*1.5*t);
       y=y';
    end

    if(mode==6)
        y=zeros(lenSamples,1);
        d=pw*Fs;
        for n=1:d:lenSamples
            y(n)=1;
        end
    end

    %Scaling and removing Ey
    y=y-mean(y);
    y=y/max(abs(y));
    y=y/10; %-20 dB

    %Apply windowing function
    winlat=fix(Fs*at);
    winl=fix(Fs*dc);
    t=0:winl-1;
    t=t/((winl-1)/(.5*pi));
    decay=(cos(t).^2)';
    t2=0:winlat-1;
    t2=t2/((winlat-1)/(.5*pi));
    decay2=(cos(t2).^2)';
    attack=flipud(decay2);
    y(1:winlat)=y(1:winlat).*attack;
    y(length(y)-winl+1:length(y))=y(length(y)-winl+1:length(y)).*decay;

    %Generating the Stereo Signal Data
    yRight=y*RightLevel;
    yLeft=y*LeftLevel;
    out=[yLeft yRight];
end