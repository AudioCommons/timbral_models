function [par, psi] = afeRAA(earSignals, fs, parConf)
%afeRAA Room Acoustic Analyzer (RAA) based on Auditory Front-End (AFE).
%   This is an implementation of the Room Acoustic Analyzer (RAA) 
%   developed by Jasper van Dorp Schuitman (2011), using the Two!Ears
%   Auditory Front-End (AFE) processors.
%   This function takes a binaural audio file as the input (.wav format), 
%   and calculates the perceptual predictors for 4 acoustical properties
%   of rooms: reverberance, clarity, apparent source width (ASW) and
%   listener envelopment (LEV). The additional input is a structure 
%   which contains a list of internal constants used to derive the 4
%   predictors. The parConf consists of the following:
% 
%         parConf.mu_psi = 0.34;            % (folowing a modification 
%         parConf.mu_psi_dip = -60.39*1e-3; % suggestedby Osses Vecchi)
% 
%         parConf.T_min = 63.1*1e-3;
% 
%         parConf.a_P_REV = 0.207; 
%         parConf.b_P_REV = 21.4;
% 
%         parConf.a_P_CLA = 2.76; 
%         parConf.b_P_CLA = 2.84;
% 
%         parConf.mu_ASW = 2e-2;      % or alpha1
%         parConf.nu_ASW = 5.63e2;    % or beta1
% 
%         parConf.a_P_ASW = 2.41;
%         parConf.b_P_ASW = 2.64;
% 
%         parConf.mu_LEV = 2.76e-2;   % or alpha2
%         parConf.nu_LEV = 6.80e2;    % or beta2
% 
%         parConf.a_P_LEV = 1.17; 
%         parConf.b_P_LEV = 2.20;
%       
%   Detailed descriptions of these parameters can be found in (van Dorp Schuitman, J., 2011)
%   and (van Dorp Schuitman, J., de Vries, D., & Lindau, A., 2013). 
%   The outputs consist of a number of intermediate outputs at various
%   stages, as well as the final outputs. 
%   Output par consists of the following:
%           FL: Level of foreground stream (MU) (*)
%           BL: Level of background stream (MU)
%           ITDf: ITD fluctuation in the foreground (sec)
%           ITDb: ITD fluctuation in the background (sec)
%           Llow: Level of the low-frequency part of the spectrum
%           pRev: Reverberance (unscaled)
%           pClar: Clarity (unscaled)
%           pASW: Apparent Source Width (unscaled)
%           pLEV: Listener Envelopment (unscaled)
%        	sRev: Reverberance (scaled)
%           sClar: Clarity (scaled)
%           sASW: Apparent Source Width (scaled)
%           sLEV: Listener Envelopment (scaled)
% (*) MU = Model Units
% 
%   Output psi consists of the following:
%           fs: Sampling frequency of the data (Hz)
%           numSamples: Number of time samples
%           itdNumFrames: Number of frames for which ITD is calculated
%           numBands: Number of frequency bands in the data
%           f: Vector with frequencies (Hz)
%           PsiL: Time-frequency representation for the left ear (Model Units)
%           PsiR: Time-frequency representation for the right ear (Model Units)
%           ITD: Inter-aural time difference as a function time and frequency
%           t_psi: Time axis for the Psi data
%           t_itd: Time axis for the ITD data
%
%   Coding by Ryan Chungeun Kim, Eindhoven University of Technology, 2016
% 
%   References:
%   van Dorp Schuitman, J. (2011). Auditory modeling for assessing 
%       room acoustics. Delft University of Technology.
%   van Dorp Schuitman, J., de Vries, D., & Lindau, A. (2013). 
%       Deriving content-specific measures of room acoustic perception 
%       using a binaural, nonlinear auditory model. The Journal of the 
%       Acoustical Society of America, 133(March), 1572–1585. 
%       http://doi.org/10.1121/1.4789357
%
%
% Note: This function uses the following script/functions for
%     internal operation (the files are in the Tools folder):
%   RAA_param_configuration.m: default parameter configuration script
%     (in case separate input parameter structure is not given)
%   RAA_group_indices.m: used to calculate the foreground/background streams
%   RAA_sigmoid_scaling.m: used to scale the model outputs


% If the second input is not given
if (nargin < 2)
    % Load default parameters as parConf structure
    run('RAA_param_configuration.m');
end

%% LOAD SIGNAL
% 
% Load a signal
psi.fs = fs;
[psi.numSamples,nChannels] = size(earSignals);

% Check if input is binaural
if nChannels ~= 2
    error('The afeRAA requires a binaural input signal.')
end

time = (1:psi.numSamples).'/(psi.fs);
psi.t_psi = time;

% Apply OME filtering as in vd Schuitman work
[b1, b2, a1, a2] = headphonefilter_Dorp2011(psi.fs);
outsig = filter(b1,a1,earSignals); % high-pass
earSignals_OME = filter(b2,a2,outsig); % low-pass

% Create a data object based on the ear signals
dObj = dataObject(earSignals_OME, psi.fs, 31, 2);

%% PLACE REQUEST AND CONTROL PARAMETERS
% 
% Request interaural time differences (ITDs)
requests = {'adaptation'};

% Parameters of preprocessing
pp_bLevelScaling = true;

% Parameters of the auditory filterbank processor
% Following van Dorp Schuitman's configuration (of frequency bands)
fb_type       = 'gammatone';
fb_lowFreqHz  = 167.7924;
fb_highFreqHz = 1836.4;
fb_nChannels  = 2;  
psi.numBands = fb_nChannels;

% Parameters of innerhaircell processor
ihc_method    = 'breebart';

% Parameters of adaptation loop processor
% ATH has not been applied (as in van Dorp Schuitman's model)
adpt_model = 'adt_vandorpschuitman';
adpt_lim = 0;       % No overshoot limitation

% Summary of parameters 
par_OME = genParStruct('pp_bLevelScaling', pp_bLevelScaling, ...
                   'fb_type',fb_type,'fb_lowFreqHz',fb_lowFreqHz,...
                   'fb_highFreqHz',fb_highFreqHz,'fb_nChannels',fb_nChannels,...
                   'ihc_method',ihc_method, ...
                   'adpt_model', adpt_model, 'adpt_lim', adpt_lim); 
               
%% PERFORM PROCESSING
% 
% Create a manager
mObj = manager(dObj, requests, par_OME);

% Request processing
mObj.processSignal();

% Grab CF for future use
cfHz = dObj.filterbank{1}.cfHz;
psi.f = cfHz.';

%% BINAURAL PROCESSING - ITD CALCULATION

% ITD frame step size (hSizeSec) is 25 ms (= numSamples/itdNumFrames/fs)
% ITD frame length (wSizeSec): 50 ms
w_wSizeSec = 0.05;
w_hSizeSec = 0.025;
wSize = 2*round(w_wSizeSec*psi.fs/2);
hSize = round(w_hSizeSec*psi.fs);
FsHzOut = 1/(w_hSizeSec);

% Specification of the double-sided exponential window
% w[n] = exp(-0.5*abs(n)/(tau_b*fs)) 
% This becomes "win" for the coming routines
tau_b = 0.03;       % time constant for the window function w[n]
n_win = (fix(-wSize/2):round(wSize/2)-1);
win = exp(-0.5*abs(n_win)/(tau_b*psi.fs)).';

% Range of tau for ITD calculation
maxDelaySec = 700*1e-6;         % 700 us

% Using adaptation loop output, without smoothing (see Figure 3.10)
in_l = dObj.adaptation{1}.Data(:);
in_r = dObj.adaptation{2}.Data(:);

[nSamples,nChannels] = size(in_l);

% How many frames are in the buffered input?
nFrames = floor((nSamples-(wSize-hSize))/hSize);
psi.itdNumFrames = nFrames;

% Determine maximum lag in samples
maxLag = ceil(maxDelaySec*psi.fs);

outITD = zeros(nFrames,nChannels);
psi.t_itd = (0:nFrames-1).'/FsHzOut;

% Loop on the time frame
for ii = 1:nFrames
    % Get start and end indexes for the current frame
    n_start = (ii-1)*hSize+1;
    n_end = (ii-1)*hSize+wSize;

    % Loop on the channel
    for jj = 1:nChannels

        % Extract frame for left and right input
        % NOTE: ITD calculation uses adaptation loop output 
        % BEFORE LPF SMOOTHING!!
        frame_l = win.*in_l(n_start:n_end,jj);
        frame_r = win.*in_r(n_start:n_end,jj);

        % Compute the frames in the Fourier domain
        X = fft(frame_l,2^nextpow2(2*wSize-1));
        Y = fft(frame_r,2^nextpow2(2*wSize-1));

        % Compute cross-power spectrum
        XY = X.*conj(Y);

        % Back to time domain
        c = real(ifft(XY));

        % Adjust to requested maximum lag and move negative
        % lags upfront
        if maxLag >= wSize
            % Then pad with zeros
            pad = zeros(maxLag-wSize+1,1);
            c = [pad;c(end-wSize+2:end);c(1:wSize);pad];
        else
            % Else keep lags lower than requested max
            c = [c(end-maxLag+1:end);c(1:maxLag+1)];
        end
        
        nLags = length(c);

        % Create a lag vector
        lags = (0:nLags-1).'-(nLags-1)/2;

        % Find the peak in the discretized crosscorrelation
        [c_peak,i] = max(c);

        % Lag of most salient peak
        lagInt = lags(i);

        if i>1 && i<nLags
            % Then interpolate using neighbor points
            c_l = c(i-1);    % Lower neighbor
            c_u = c(i+1);    % Upper neighbor

            % Estimate "true" peak deviation through parabolic
            % interpolation
            delta = 0.5*(c_l-c_u)/(c_l-2*c_peak+c_u);

            % Store estimate
            outITD(ii,jj) = (lagInt + delta)/FsHzOut;

        else
            % Do not interpolate if the peak is at a boundary
            outITD(ii,jj) = lagInt/FsHzOut;
        end

    end

end
psi.ITD = outITD;

%% MONAURAL PROCESSING - LPF smoothing

% LPF 1st order, 8 Hz cutoff frequency (Sec. 3.5.1)
% 20ms time constant (Eq. A.23), 1st order difference equation
% y[n] = (1-exp(-1/(tau_e*fs))*x[n] + exp(-1/(tau_e*fs))*y[n-1]

tau_e = 20*1e-3;

bEnvExtFilter = 1-exp(-1/(tau_e*psi.fs));
aEnvExtFilter = [1 -exp(-1/(tau_e*psi.fs))];

PsiL = filter(bEnvExtFilter, aEnvExtFilter, in_l);
PsiR = filter(bEnvExtFilter, aEnvExtFilter, in_r);
psi.PsiL = PsiL;
psi.PsiR = PsiR;

%% PEAK / DIP DETECTION
% Constants for peak/dip detection
% mu_psi = 7.49*1e-3;       % vd Schuitman's original values
% mu_psi_dip = -1.33*1e-3;    % vd Schuitman's
mu_psi = 0.34;              % Osses Vecchi's modified values
mu_psi_dip = -60.39*1e-3;    

T_min = 63.1*1e-3;

% use nSamples and nChannels
outDir_l  = zeros(nSamples, nChannels);
outDir_r  = zeros(nSamples, nChannels);

% Eq. 3.15 avg. absolute level
L_psi_l = 1/nSamples*sum(abs(in_l), 1); % dimension: (1 x length(cfHz))
L_psi_r = 1/nSamples*sum(abs(in_r), 1); 

% Eq. 3.14 threshold per frequency band
psi_min_l = parConf.mu_psi * L_psi_l;
psi_min_r = parConf.mu_psi * L_psi_r;

psi_min_dip_l = parConf.mu_psi_dip * L_psi_l;
psi_min_dip_r = parConf.mu_psi_dip * L_psi_r;

% Per frequency channel
for f=1:nChannels
    
    % Find time indices at which PsiL and PsiR are above threshold for peak
    idxAboveThreshold_l = find(PsiL(:, f) >= psi_min_l(f));
    idxAboveThreshold_r = find(PsiR(:, f) >= psi_min_r(f));
    
    % Group consecutive time indices of threshold-passing inputs
    consecIdxGroups_above_thr_l = RAA_group_indices(idxAboveThreshold_l);
    consecIdxGroups_above_thr_r = RAA_group_indices(idxAboveThreshold_r);
    % Groups of indices are returned as cells
    
    % Check if the duration of each index group exceeds parConf.T_min
    % Left and right separately
    % Left channel
    for gg=1:size(consecIdxGroups_above_thr_l, 1)
        dur = (consecIdxGroups_above_thr_l{gg}(end) - consecIdxGroups_above_thr_l{gg}(1))/psi.fs;
        % If this duration exceeds parConf.T_min then determine that group as a peak
        % and assign the corresponding input to the direct stream
        if(dur >= parConf.T_min)
            outDir_l(consecIdxGroups_above_thr_l{gg}, f) = ...
                PsiL(consecIdxGroups_above_thr_l{gg}, f);
        end
    end
    % Right channel
    for gg=1:size(consecIdxGroups_above_thr_r, 1)
        dur = (consecIdxGroups_above_thr_r{gg}(end) - consecIdxGroups_above_thr_r{gg}(1))/psi.fs;
        % If this duration exceeds parConf.T_min then determine that group as a peak
        % and assign the corresponding input to the direct stream
        if(dur >= parConf.T_min)
            outDir_r(consecIdxGroups_above_thr_r{gg}, f) = ...
                PsiR(consecIdxGroups_above_thr_r{gg}, f);
        end
    end
   
    % Find time indices at which PsiL and PsiR are below threshold for dip
    idxBelowThreshold_l = find(PsiL(:, f) <= psi_min_dip_l(f));
    idxBelowThreshold_r = find(PsiR(:, f) <= psi_min_dip_r(f));
    
    % Group consecutive time indices of threshold-passing inputs
    consecIdxGroups_below_thr_l = RAA_group_indices(idxBelowThreshold_l);
    consecIdxGroups_below_thr_r = RAA_group_indices(idxBelowThreshold_r);
    % Groups of indices are returned as cells
        
    % Check if the duration of each index group exceeds parConf.T_min
    % Left and right separately
    for gg=1:size(consecIdxGroups_below_thr_l, 1)
        dur = (consecIdxGroups_below_thr_l{gg}(end) - consecIdxGroups_below_thr_l{gg}(1))/psi.fs;
        % If this duration exceeds parConf.T_min then determine that group as a dip
        % and assign the corresponding input to the direct stream
        if(dur >= parConf.T_min)
            outDir_l(consecIdxGroups_below_thr_l{gg}, f) = ...
                PsiL(consecIdxGroups_below_thr_l{gg}, f);
        end
    end
    for gg=1:size(consecIdxGroups_below_thr_r, 1)
        dur = (consecIdxGroups_below_thr_r{gg}(end) - consecIdxGroups_below_thr_r{gg}(1))/psi.fs;
        % If this duration exceeds parConf.T_min then determine that group as a dip
        % and assign the corresponding input to the direct stream
        if(dur >= parConf.T_min)
            outDir_r(consecIdxGroups_below_thr_r{gg}, f) = ...
                PsiR(consecIdxGroups_below_thr_r{gg}, f);
        end
    end       
        
end

% The rest of the input signals become the reverberant streams
outRev_l = PsiL - outDir_l;
outRev_r = PsiR - outDir_r;
% outDir_l means Psi_L,direct (in the references)
% outDir_r means Psi_R,direct
% outRev_l means Psi_L,reverberant
% outRev_r means Psi_R,reverberant
psi.PsiLdir = outDir_l;
psi.PsiRdir = outDir_r;
psi.PsiLrev = outRev_l;
psi.PsiRrev = outRev_r;

% % Comment as necessary%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot streams for test comparison
% figure;
% subplot(3, 1, 1)
% plot(time, PsiL(:, 7));
% hold on; plot(time, psi_min_l(7)*ones(length(time), 1));
% hold on; plot(time, psi_min_dip_l(7)*ones(length(time), 1));
% subplot(3, 1, 2)
% plot(time, outDir_l(:, 7));
% subplot(3, 1, 3)
% plot(time, outRev_l(:, 7));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ITD SEGREGATION
% Compare ITD time indices to outDir and outRev indices

% 1s for nonzero time indices of L/R direct streams
idxL_dir = (outDir_l ~= 0);
idxR_dir = (outDir_r ~= 0);
zeta_DIR_t = idxL_dir+idxR_dir;
% zeta_DIR_t will be
% 0: for totally reverberant part
% 1: for parts with only left or right is in direct stream
% 2: for totally direct part
zeta_DIR_t = zeta_DIR_t/2;
% Now zeta_DIR_t will be 0, 0.5, or 1

% Allocate memory for zeta_DIR in ITD TIME FRAMES
zeta_DIR_itd = zeros(nFrames, nChannels);

% % ITD segregation method #1 (Ryan)
% % Run through zeta_DIR_t by ITD frames again
% for ii = 1:nFrames
%     % Get start and end indices for the current frame
%     % ITD frame length in the representation is hSize
%     n_start = (ii-1)*hSize+1;
%     n_end = ii*hSize;
%     
%     % Take one frame from zeta_DIR_t
%     frame_zeta_DIR = zeta_DIR_t(n_start:n_end, :);
%     % Equivalent to calculating the ratio of non-zero parts within frame
%     temp = frame_zeta_DIR .* ones(size(frame_zeta_DIR))/hSize;
%     % ratios of 0s, 0.5s, and 1s within the frame are multiplied by the
%     % actual zeta values
%     zeta_DIR_itd(ii, :) = sum(temp, 1);
% end

% ITD segregation method #2 (Original)
% Run through zeta_DIR_t by ITD frames again
for ii = 1:nFrames
    % Get start indices for the current frame
    % ITD frame length in the representation is hSize
    n_start = (ii-1)*hSize+1;
    
    % Check the centre of the frame - hSize/2? or wSize/2?
    % wSize/2 seems to make more sense, because
    % the centre of the frame was weighted the most in ITD calculation
    % (by the windowing process)
    % resample zeta_DIR_t at the time indices at centres of frames
    zeta_DIR_itd(ii, :) = ...
        zeta_DIR_t(n_start + round(wSize/2), :);
    
end

% ITD_DIR and ITD_REV calculation (eqs. 3.17 and 3.18)
outITD_DIR = zeta_DIR_itd.*outITD;
outITD_REV = outITD - outITD_DIR;
psi.ITDdir = outITD_DIR;
psi.ITDrev = outITD_REV;

%% Reverberance
% Using Psi_L,rev and Psi_R,rev
% (outRev_l and outRev_r)
% Eq. 3.20

L_REV = sum(sum(sqrt(outRev_l.^2 + outRev_r.^2))) / ...
    (nSamples*nChannels);
par.BL = L_REV;
P_REV = L_REV;

% Scaling (Eq. 6.4, Table 6.7, error corrected by vd Schuitman)
% * PREV: (a,b) = (0.207, 21.4)
% a_P_REV = 0.207; b_P_REV = 21.4;
S_REV = RAA_sigmoid_scaling(parConf.a_P_REV, parConf.b_P_REV, P_REV);
par.sREV = S_REV;

%% Clarity
% Using Psi_L,dir and Psi_R,dir
% (outDir_l and outDir_r)
% Eq. 3.21, 3.22

L_DIR = sum(sum(sqrt(outDir_l.^2 + outDir_r.^2))) / ...
    (nSamples*nChannels);
par.FL = L_DIR;
P_CLA = L_DIR/L_REV;
par.pClar = P_CLA;

% Scaling (Eq. 6.4, Table 6.7, error corrected by vd Schuitman)
% * PCLA: (a,b) = (1.17, 2.20)
% a_P_CLA = 1.17; b_P_CLA = 2.2;
S_CLA = RAA_sigmoid_scaling(parConf.a_P_CLA, parConf.b_P_CLA, P_CLA);
par.sClar = S_CLA;

%% ASW
% Needs two constants: mu_ASW (alpha1), nu_ASW (beta1)
% Eqs. 3.23 - 3.25
% mu_ASW = 2e-2;      % or alpha1
% nu_ASW = 5.63e2;    % or beta1

% Intermediate values - L_LOW, sigma_ITD,DIR
% L_LOW: average output stream level at low frequencies
% here the low frequency range (indices z0 and z1) must be specified
% According to the JASA paper, z0 is at 168Hz, and z1 at 387Hz
% Look back at cfHz and find the corresponding indices - in this case 
% cfHz(1) = 168, and cfHz(5) = 387
[min0, z0] = min(abs(cfHz - 168));
[min1, z1] = min(abs(cfHz - 387));
L_LOW = sum(sum(sqrt(PsiL(:, z0:z1).^2 + PsiR(:, z0:z1).^2))) / ...
    (nSamples*(z1-z0+1));
par.Llow = L_LOW;

% sigma_ITD,DIR (avg. standard deviation of direct ITD stream)
% Here the frequency range is specified with q0 and q1 (cfHz indices)
% According to the JASA paper, q0 is at 387Hz, and q1 at 1.84kHz
% cfHz(5) = 387, cfHz(16) = 1836.4
[min0, q0] = min(abs(cfHz - 387));
[min1, q1] = min(abs(cfHz - 1836.4));
% STD over time (frame), then average over frequency channels
sigma_ITD_DIR = mean(std(outITD_DIR(:, q0:q1)));
par.ITDf = sigma_ITD_DIR;

P_ASW = parConf.mu_ASW*L_LOW + log10(1 + parConf.nu_ASW*sigma_ITD_DIR*1e3);
par.pASW = P_ASW;

% Scaling (Eq. 6.4, Table 6.7, error corrected by vd Schuitman)
% a_P_ASW = 2.76; b_P_ASW = 2.84;
S_ASW = RAA_sigmoid_scaling(parConf.a_P_ASW, parConf.b_P_ASW, P_ASW);
par.sASW = S_ASW;

%% LEV (Eqs. 3.26, 3.27)
% Also needs two constants: mu_LEV (alpha2), nu_LEV(beta2)
% mu_LEV = 2.76e-2;   % or alpha2
% nu_LEV = 6.80e2;    % or beta2

% intermediate values - L_REV (already calculated for reverberance), and
% sigma_ITD_REV
% L_REV: mean reverberant stream level
% sigma_ITD_REV: avg. standard deviation of reverberant ITD stream
% Here the frequency range is specified with q0 and q1 (cfHz indices)
% According to the JASA paper, q0 is at 387Hz, and q1 at 1.84kHz
% cfHz(5) = 387, cfHz(16) = 1836.4
[min0, q0] = min(abs(cfHz - 387));
[min1, q1] = min(abs(cfHz - 1836.4));
% STD over time (frame), then average over frequency channels
sigma_ITD_REV = mean(std(outITD_REV(:, q0:q1)));
par.ITDb = sigma_ITD_REV;

P_LEV = parConf.mu_LEV*L_REV + log10(1 + parConf.nu_LEV*sigma_ITD_REV*1e3);
par.pLEV = P_LEV;

% Scaling (Eq. 6.4, Table 6.7, error corrected by vd Schuitman)
% * PLEV: (a,b) = (2.41, 2.64)
% a_P_LEV = 2.41; b_P_LEV = 2.64;
S_LEV = RAA_sigmoid_scaling(parConf.a_P_LEV, parConf.b_P_LEV, P_LEV);
par.sLEV = S_LEV;

end
