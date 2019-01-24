function [b,a,bwHz] = createFB_Mod(fsHz,cfModHz,Q,bDown2DC,bUp2Nyquist)
%createFB_Mod   Modulation filterbank based on second-order filters [1,2].
% 
%USAGE
%    [b,a,bwHz] = createFB_Mod(fsHz)
%    [b,a,bwHz] = createFB_Mod(fsHz,cfModHz,Q,bDown2DC,bUp2Nyquist)
%
%INPUT ARGUMENTS
%          fsHz : sampling frequency in Hz
%       cfModHz : center frequencies of modulation filters in Hertz
%                 (default, cfModHz = pow2(0:10))
%             Q : quality factor, either a scalar or a vector, specifying
%                 the selectivity of each modulation filter (default, Q = 1)
%      bDown2DC : implement the first modulation filter as a low-pass
%                 (default, bDown2DC = true) 
%   bUp2Nyquist : implement the last modulation filter as a high-pass
%                 (default, bUp2Nyquist = false) 
% 
%OUTPUT ARGUMENTS
%             b : cell array of numerator coefficients
%             a : cell array of denominator coefficients 
%          bwHz : -3 dB filter bandwidth of modulation filters in Hertz
% 
%   Use createFreqAxisLog to create a vector with log-scaped center
%   frequencies.
% 
%   createFB_Mod(...) plots the transfer function of the modulation
%   filterbank in a new figure.  
% 
%   See also createFreqAxisLog, createFB_MEL and createFB_Oct.
% 
%EXAMPLE
%   % Design and visualze a modulation filterbank at 16 kHz sampling
%   % frequency
%   createFB_Mod(16E3);
% 
%REFERENCES
%   [1] T. May and T. Dau, "Computational speech segregation based on an
%       auditory-inspired modulation analysis", The Journal of the
%       Acoustical Society of America, 136(6), pp.3350-3359, 2014.   
% 
%   [2] S. D. Ewert and T. Dau, "Characterizing frequency selectivity for
%       envelope fluctuations", The Journal of the Acoustical Society of
%       America, 108(3), pp.1181-1196, 2000.    


%   Developed with Matlab 8.3.0.532 (R2014a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2014
%              Technical University of Denmark (DTU)
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2013/07/06
%   v.0.2   2014/07/10 incorporated flag "bUp2Nyquist"
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS  
% 
% 
% Check for proper input arguments
if nargin < 1 || nargin > 5
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 5 || isempty(bUp2Nyquist); bUp2Nyquist = false;      end
if nargin < 4 || isempty(bDown2DC);    bDown2DC    = true;       end
if nargin < 3 || isempty(Q);           Q           = 1;          end
if nargin < 2 || isempty(cfModHz);     cfModHz     = pow2(0:10); end

% Filter order of bandpass filters
filterOrder = 2;

% Check filter order
if filterOrder ~= (round(filterOrder / 2) * 2)
    error('Filter order must be even-numbered!')    
else
    % Order of band-pass filters is 2 * N, thus N = filterOrder / 2 
    N = filterOrder / 2;
end

% Number of modulation filter
nFilter = numel(cfModHz);

% Replicate Q
if isscalar(Q)
    Q = repmat(Q,[nFilter 1]);
else
    if nFilter ~= numel(Q)
        error('Q factor and cfMod must be of equal size.')
    end
end


%% CREATE FILTERBANK
% 
% 
% Allocate memory
[w0,bwHz,f1,f2] = deal(zeros(nFilter,1));
[b,a]           = deal(cell(nFilter,1));
wn              = zeros(nFilter,2);

% Nyquist frequency
nyquistHz = fsHz/2;

% Loop over number of modulation filters
for ii = 1 : nFilter
    
    % Normalized center frequency (0,1), whereas 1 => fs/2
    w0(ii) = cfModHz(ii)/nyquistHz;
    
    % Check if center frequency is valid
    if w0(ii) >= 1
        error('Center frequencies must be smaller than fs/2')
    end
    
    % Filter bandwidth
    bwHz(ii) = w0(ii) * nyquistHz / Q(ii);

    % Passband frequencies 
    f1(ii) = cfModHz(ii) * (sqrt(1+(1/(4*Q(ii)^2))) - (1/(2*Q(ii))));
    f2(ii) = cfModHz(ii) * (sqrt(1+(1/(4*Q(ii)^2))) + (1/(2*Q(ii))));

    % Normalized passband frequencies (0,1), whereas 1 => fs/2
    wn(ii,:) = [f1(ii) f2(ii)] / nyquistHz;

    % Check if passband frequencies wn(ii,:) are valid
    if any(wn(ii,:) >= 1)
        error('Passband frequency must be smaller than fs/2')
    end
    
    % Derive filter coefficients
    if ii == 1 && bDown2DC
        [b{ii},a{ii}] = butter(N,w0(ii),'low');
    elseif ii == nFilter && bUp2Nyquist
        [b{ii},a{ii}] = butter(N,w0(ii),'high');
    else
        % Filter order is 2 * N
        [b{ii},a{ii}] = butter(N,wn(ii,:),'bandpass');
    end
    
    % Check if resulting filter coefficients are stable
    if ~isStable(b{ii},a{ii})
        error('IIR filter is not stable, reduce the filter order!')
    end    
end


%% PLOT FILTERBANK
% 
% 
% Plot filter transfer function
if nargout == 0 
    
    lineWidth = 2;
    strCol = {[0 0 0 ] [0.45 0.45 0.45]};
    figure; hold on;
    for ii = 1:nFilter
        [TF,freqs] = freqz(b{ii},a{ii},96E3,fsHz);
        plot(freqs,20*log10(abs(TF)),'linewidth',lineWidth,'color',strCol{1+mod(ii-1,2)});
    end
    
    set(gca,'xscale','log')
    title('Modulation filterbank')
    xlabel('Frequency (Hz)')
    ylabel('Filter attenuation (dB)')
    xlim([0 max(400,4*max(cfModHz))])
    ylim([-40 5])
end
