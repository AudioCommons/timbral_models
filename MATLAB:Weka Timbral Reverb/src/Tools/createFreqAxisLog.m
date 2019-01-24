function [cfHz,nFreq] = createFreqAxisLog(lowFreqHz,highFreqHz,nFreq)
%createFreqAxisLog   Create log-spaced frequency axis.
% 
%USAGE
%   [cfHz,nFreq] = createFreqAxisLog(lowFreqHz,highFreqHz)
%   [cfHz,nFreq] = createFreqAxisLog(lowFreqHz,highFreqHz,nFreq)
% 
%INPUT ARGUMENTS
%      lowFreqHz : low frequency limit in Hertz
%     highFreqHz : high frequency limit in Hertz
%          nFreq : number of frequencies (default, nFreq = [])
% 
%OUTPUT ARGUMENTS
%           cfHz : frequency axis in Hertz
%          nFreq : number of frequencies
% 
%   Starting at the lower frequency limit, frequencies are equally-spaced
%   on a log scale, until the high frequency limit is reached. Due to the
%   constant spacing, the highest frequency may not be included (e.g.
%   createFreqAxisLog(1,40)) 
% 
%   If a particular number of frequencies is requested, the spacing between
%   frequencies is adjusted to include both the low and the high frequency
%   limits (e.g. createFreqAxisLog(1,40,8)).

%   Developed with Matlab 8.3.0.532 (R2014a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2014
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2014/11/26
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 3
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 3; nFreq = []; end

% Check frequency range
if lowFreqHz <= 0
    error('Lower frequency limit must be larger than 0 Hz.')
end


%% CREATE LOG-SPACED FREQUENCY AXIS
% 
% 
if ~isempty(lowFreqHz) && ~isempty(highFreqHz) && ~isempty(nFreq)
    % 1. Frequency range and number of filters specified
    cfHz = pow2(linspace(log2(lowFreqHz),log2(highFreqHz),nFreq));
elseif ~isempty(lowFreqHz) && ~isempty(highFreqHz) && isempty(nFreq)
    % 2. Only frequency range is specified
    cfHz  = pow2(log2(lowFreqHz):log2(highFreqHz));
    nFreq = numel(cfHz);
else
    error('Not enough or incoherent input arguments.')
end