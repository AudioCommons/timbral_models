function [CC,ITD,ILD,lagL,lagR,BL,BR,LLARmode,cc]=prec_acmod(x1,x2,Fs,cfHz,maxLag,cc)

% function [ITD,ILD,lagL,lagR,BL,BR,LLARmode,info]=acmod(y,Fs);
%
% Implementation of a precedence effect model to simulate localization 
% dominance using an adaptive, stimulus parameter-based inhibition process
% 
% J. Braasch (2013) J. Acoust. Soc. Am. 134(1): 420-435.
% For section references see paper
%
% (c) 2013, Rensselaer Polytechnic Institute
% contact: jonasbraasch@gmail.com
%
% Modified by Ryan Chungeun Kim for Two!Ears software framework, 2015
% TODO: If found to be adequate, make the A1/A2 and B1/B2 options as
% additional input parameter
%
% INPUT PARAMETERS:
% x1, x2    : binaural input signal, time x frequency domain (after filterbank)
% Fs        : sampling frequency (tested for 48 kHz)
% cfHz      : vector of auditory filterbank centre frequencies (e.g., gammatone filterbank)
% maxLag    : maximum lag for correlation calculation
% cc        : structure to save various internal/intermediate calculation results

% OUTPUT PARAMETERS:
% ITD      = Interaural Time Difference collapsed over time and frequency
% ILD      = Interaural Level Difference collapsed over time and frequency
% lagL     = lag delay for left lag [ms]
% lagR     = lag delay for right lag [ms]
% BL       = left-lag amplitude, relative according to Eq. A5    
% BR       = right-lag amplitude, relative according to Eq. A5    
% LLARmode = Lag-to-Lead Amplitude Ratio mode, see Table 1

% DEPENDENT FUNCTIONS:
% prec_anaone.m    = analysis and lag removal for single channel
% calcXCorr.m      = Cross-correlation calculation

Fms=Fs./1000; % sampling frequency based on milliseconds
lagFms = ceil(Fms);     % lag to use for XCorr calculation (samples for 1ms)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monaural preprocessing (Section II.E) ^
% Lag removal (Section II.B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(cc.ac1) % new instance 
    % LAG FOR AUTOCORRELATION IS FIXED AT 10 ms! (Communication with Braasch)
    [yL_LgS,yL_LgL,lagL,BL,cc.ac1]=prec_anaone(x1, Fs, cfHz, round(0.01*Fs)); % left channel
    [yR_LgS,yR_LgL,lagR,BR,cc.ac2]=prec_anaone(x2, Fs, cfHz, round(0.01*Fs)); % right channel
else
    [yL_LgS,yL_LgL,lagL,BL,cc.ac1]=prec_anaone(x1, Fs, cfHz, round(0.01*Fs), cc.ac1); % left channel
    [yR_LgS,yR_LgL,lagR,BR,cc.ac2]=prec_anaone(x2, Fs, cfHz, round(0.01*Fs), cc.ac2); % right channel

end    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine correct LLAR mode
% Section II.D, Table 1 & Fig. 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(cc.cc1)  % ??
    cc.cc1=calcXCorr(sum(yL_LgS).',sum(yR_LgS).',lagFms); % calc ICC ofr LLAR mode 1
    cc.cc2=calcXCorr(sum(yL_LgL).',sum(yR_LgL).',lagFms); % calc ICC ofr LLAR mode 2
    cc.cc3=calcXCorr(sum(yL_LgS).',sum(yR_LgL).',lagFms); % calc ICC ofr LLAR mode 3
    cc.cc4=calcXCorr(sum(yL_LgL).',sum(yR_LgS).',lagFms); % calc ICC ofr LLAR mode 4  
else 
    cc.cc1=cc.cc1+calcXCorr(sum(yL_LgS).',sum(yR_LgS).',lagFms); % calc ICC ofr LLAR mode 1
    cc.cc2=cc.cc2+calcXCorr(sum(yL_LgL).',sum(yR_LgL).',lagFms); % calc ICC ofr LLAR mode 2
    cc.cc3=cc.cc3+calcXCorr(sum(yL_LgS).',sum(yR_LgL).',lagFms); % calc ICC ofr LLAR mode 3
    cc.cc4=cc.cc4+calcXCorr(sum(yL_LgL).',sum(yR_LgS).',lagFms); % calc ICC ofr LLAR mode 4
end

% cross-correlation to determine best combination
n1=max(cc.cc1); % calc ICC ofr LLAR mode 1
n2=max(cc.cc2); % calc ICC ofr LLAR mode 2
n3=max(cc.cc3); % calc ICC ofr LLAR mode 3
n4=max(cc.cc4); % calc ICC ofr LLAR mode 4
% Pick LLARmode based on highest coherence
[maxi,LLARmode]=max([n1 n2 n3 n4]);

% assign left and right signals x1/x2 from best LLAR mode
switch LLARmode
    case 1
        x1=yL_LgS;
        x2=yR_LgS;
    case 2
        x1=yL_LgL;
        x2=yR_LgL;
    case 3
        x1=yL_LgS;
        x2=yR_LgL;
    case 4
        x1=yL_LgL;
        x2=yR_LgS;
end % switch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determination of binaural cues in individual frequency bands
% Section III.A, see Box "CE?in Fig. 9
%
% Note: In this model version we do NOT select LLAR Mode 1 
% automatically if all four correlation values
% differ by less than 0.1 (p. 428)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ICCintT = zeros(length(cfHz), 2*maxLag+1);
ILDint = zeros(length(cfHz), 1);
Eint = zeros(length(cfHz), 1);

for n=1:length(cfHz) % loop over all frequency bands
    xL=x1(n,:)';
    xR=x2(n,:)';

%         ICC=calcXCorr(xL,xR,maxLag)'; % interaural cross-correlation
        ICC=calcXCorr(xL,xR,maxLag, 'coeff')'; % interaural cross-correlation, with normalisation option
        
%         eL=mean(sqrt(xL.^2)); % rms energy in left channel
%         eR=mean(sqrt(xR.^2)); % rms energy in right channel 
        % aren't the above (original acmod) wrong for rms? - corrected below
        eL=sqrt(mean(xL.^2)); % rms energy in left channel
        eR=sqrt(mean(xR.^2)); % rms energy in right channel
        if eL>0 && eR>0; % if signals in both channels exist
            % right/left order is switched from Braasch's original version
            % for AFE ILD comparability
            ILD_f = 20*log10(eR./eL); % ILD within frequency band
            En=eL.*eR; % Amplitude
        else 
            ILD_f=0;
            En=0;
        end % of if                    

    % Frequency weighting according to Stern et al. (1988)
    % additional comment from Braasch: "only important for conditions with
    % very narrow bandwidths"
    % frequency weighting filter coefficients: b1, b2, b3
    b1=-0.0938272;
    b2=0.000112586;
    b3=-0.0000000399154;
    f=cfHz(n);
    if f<1200
        w=10.^((-b1.*f-b2.*f.^2-b3.*f.^3)./10)./276;
    else
        w=10.^((-b1.*1200-b2.*1200.^2-b3.*1200.^3)./10)./276;
    end % of if 
    ICCintT(n,:)=ICC.*w; % Apply freq. weighting to ICC
   % whos
    ILDint(n)=ILD_f.*En./sum(En); % Energy weighted ILD
    Eint(n)=En; % integrated Energy
end % of for 

% The following corresponds to the "cumulative" cross-correlation and ILD
% calculation, but this is not necessary
% Comment/uncomment and use/test as desired
%%%%%%
% A1: from Braasch's original code - cumulative "sum"
%%%%%%
% if isempty(cc.ICCintT)      % when called for the very 1st time!
%     cc.ICCintT=ICCintT;     % Dimension: [nChannels x 2*MAXLAGS+1]
%     cc.ILDint=ILDint;
%     cc.Eint=Eint;
%     cc.lastFrame = 1;       % this will actually be 1 always
% else 
%     cc.ICCintT=cc.ICCintT+ICCintT;
%     cc.ILDint=(cc.ILDint.*cc.Eint+ILDint.*Eint)./(Eint+cc.Eint);
%     cc.Eint=Eint+cc.Eint;
%     cc.lastFrame = cc.lastFrame + 1;
% end

%%%%%%
% A2: with "cumulative moving average" version - this basically
% corresponds to LPF for fluctuating CC
%%%%%%
if isempty(cc.ICCintT)      % when called for the very 1st time!
    cc.ICCintT=ICCintT;     % Dimension: [nChannels x 2*MAXLAGS+1]
    cc.ILDint=ILDint;
    cc.Eint=Eint;
    cc.lastFrame = 1;       % this will actually be 1 always
else
    % Average CC over the cumulative number of frames so far
    cc.ICCintT=(cc.ICCintT+ICCintT)./(cc.lastFrame+1);
    cc.ILDint=(cc.ILDint.*cc.Eint+ILDint.*Eint)./(Eint+cc.Eint);
    cc.ILD=ILDint;
    cc.Eint=Eint+cc.Eint;
    cc.lastFrame = cc.lastFrame + 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determination of binaural cues across frequency bands
% See Section III.A.1. Decision device
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following corresponds to the derivation of outputs using or not using
%   the cumulative calculation made in option A part
% Comment/uncomment and use/test as desired
%%%%%%%
% B1: CC, ITD and ILD based on cumulative calculation (using cc structure) 
%%%%%%%
% CC = cc.ICCintT;                      % cumulative CC
% [maxCorr,indexITD]=max(CC, [], 2);
% ITD=(indexITD-maxLag-1)./Fs;          % return ITD in seconds
% ILD=cc.ILDint;

%%%%%%%
% B2: CC, ITD and ILD without cumulative calculation  
% (copying directly after calculation per frame instead of using cc structure) 
%%%%%%%
CC = ICCintT;                           % CC per frame
[maxCorr,indexITD]=max(CC, [], 2);
ITD=(indexITD-maxLag-1)./Fs;            % return ITD in seconds
ILD=ILDint;

% convert lag delays to ms:
lagL=lagL./Fms;
lagR=lagR./Fms;

