function bActivity = detectOnsetsOffsets(input,stepSizeSec,minStrengthdB,minSpread,fuseEventsWithinSec)


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 5
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default parameter
if nargin < 3 || isempty(minStrengthdB); minStrengthdB = 3; end
if nargin < 4 || isempty(minSpread);     minSpread     = 5; end
if nargin < 5 || isempty(fuseEventsWithinSec); fuseEventsWithinSec = 30E-3; end


%% **************************  DETECT ACTIVITY  ***************************
% 
% 
% Determine size of onset strength
[nFrames,nChannels] = size(input);

% Compute minimum distance of local peaks 
fuseNframes = round((fuseEventsWithinSec/stepSizeSec));

% Delete activity below "minStrengthdB"
input(input < minStrengthdB) = 0;

% Allocate memory 
bActivity = false(nFrames,nChannels);

% Loop over number of channels
for ii = 1 : nChannels

    % Detect local peaks
    [peakIdx,peakVal] = findpeaks_VB(input(:,ii));

    % Select peaks that are separated by at least "fuseNframes"
    peakIdx = selectPeaks(peakIdx,peakVal,fuseNframes);
    
    % Populate activity map
    bActivity(peakIdx,ii) = true;
end


%% ********************  EXTRACT CONNECTED FRAGMENTS  *********************
% 
% 
% Extract connected onset or offset fronts which spread across at least
% 'minOnsetSpread' adjacent channels. 
if ~license('test', 'image_toolbox')
    error(['Sorry, but the Image Processing Toolbox is required to ',...
        'extract coherent onsets and offsets.']);
else
    [L,nFragments] = bwlabel(bActivity,8);
end

% Loop over number of connected fragments
for ii = 1 : nFragments
    
    % Extract ii-th fragment
    idx = find(L==ii);
    
    % Check if activity spreads across at least 'minOnsetSpread' channels
    if numel(idx) < minSpread
       % Remove fragment 
       bActivity(idx) = false;
    end
end

