function [pitchHz, confidence, confThresPerc, pitchHzRaw] = estPitch_SACF(acf,lagsSec,pitchRangeHz,confThresPerc,orderMedFilt)

%   Developed with Matlab 8.2.0.701 (R2013b). Please send bug reports to:
%   
%   Authors :  Tobias May © 2014
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
% 
%   History :  
%   v.0.1   2014/03/05
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS 
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 5
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default parameter
if nargin < 3 || isempty(pitchRangeHz);  pitchRangeHz  = [80 400]; end
if nargin < 4 || isempty(confThresPerc); confThresPerc = 0.7;      end
if nargin < 5 || isempty(orderMedFilt);  orderMedFilt  = 3;        end


%% RESTRICT LAG RANGE
% 
% 
% Compute summary ACF
sacf = squeeze(mean(acf,2));

% Input size
[nFrames,nLags] = size(sacf);

% Restrict lags to plausible pitch range
rangeLagSec = 1./pitchRangeHz;

% Find corresponding lags
bValidLags = lagsSec >= min(rangeLagSec) & lagsSec <= min(max(rangeLagSec),nLags);

% Restrict lags to predefined pitch range
sacf = sacf(:,bValidLags);
lagsSec = lagsSec(bValidLags);


%% DETECT PITCH CANDIDATES
% 
% 
% Allocate memory
pitchHzRaw = zeros(nFrames,1);
confidence = zeros(nFrames,1);

% Loop over number of frames
for ii = 1 : nFrames
    
    % Detect local peaks
    [peakIdx,peakVal] = findpeaks_VB(sacf(ii,:));
        
    % Find maximum peak position and confidence value
    [maxVal,maxIdx] = max(peakVal);
    
    % Confidence value
    if isempty(maxVal)
        maxVal = 0;
    else
        confidence(ii) = maxVal;
    end
    
    % Only accept pitch estimate if confidence value is above 0
    if maxVal > 0
        % Pitch estimate in Hertz
        pitchHzRaw(ii) = 1/lagsSec(peakIdx(maxIdx));
    end
end


%% POST-PROCESSING
% 
% 
% Floor confidence value
confidence = max(confidence,0);

% Compute threshold
confThresPerc = max(confidence,[],1) * confThresPerc;

% Apply confidence threshold
bSetToZero = confidence < repmat(confThresPerc,[nFrames 1]);

% Set unreliable pitch estimates to zero
pitchHz = pitchHzRaw; pitchHz(bSetToZero) = 0;


%% POST-PROCESSING
% 
% 
% Apply median filtering to reduce octave errors
pitchHz = medfilt1(pitchHz,orderMedFilt);

% Replace all zeros with NANs
pitchHz(pitchHz==0) = NaN;

if nargout == 0 
    
   figure;
   ax(1) = subplot(211);
   plot(pitchHz)
   ax(2) = subplot(212);
   plot(confidence)
   hold on;
   plot([1 nFrames],[confThresPerc confThresPerc],'k--')
   
   linkaxes(ax,'x')
   axis tight;
end
