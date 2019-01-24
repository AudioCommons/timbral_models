function waveplot(data,xAxis,yAxis,zoom,bNorm)
%waveplot   Waterfall plot. 
%
%USAGE
%   waveplot(data)
%   waveplot(data,xAxis,yAxis,zoom,bNorm)
%   
%INPUT ARGUMENTS
%      data : data matrix arranged as  [nSamples  x nChannels]
%     xAxis : vector specifying x axis [nSamples  x 1]
%     yAxis : vector defining y axis   [nChannels x 1]
%      zoom : zoom factor              (default, zoom = 5)
%     bNorm : equalize the dynamic range of each channel prior to plotting 
%             (default, bNORM = true)
% 
%NOTE
%   This function is inspired by "waveplot", a MATLAB function
%   available at http://www.casabook.org/.
%
%EXAMPLE
%   nSamples = 500;
%   % Initialize gammatone parameter structure
%   GFB = createFB_GFB(20E3);
%   % Filter impulse with gammatone filtering 
%   fbOut = applyFB_GFB(GFB,[1; zeros(nSamples-1,1)]);
%   % Plot result
%   waveplot(fbOut,(1:nSamples)/20E3,GFB.cfHz);

%   Developed with Matlab 7.4.0.287 (R2007a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2008-2014
%              TUe Eindhoven and Philips Research  
%              t.may@tue.nl      tobias.may@philips.com
%
%   History :  
%   v.0.1   2008/05/11
%   v.0.2   2014/07/06 interpolate cf labels for y-axis
%           2014/10/22 RD: removed axis handle input argument for
%                      integration in the Two!Ears framework.
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
if nargin < 1 || nargin > 6
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default parameter
if nargin < 4 ||isempty(zoom);  zoom  = 5;    end
if nargin < 5 ||isempty(bNorm); bNorm = true; end

% Desired number of y-axis labels
nYLabels = 8;


%% ****************************  CREATE AXES  *****************************
% 
% 
% Determine size of input
[nSamples,nChannels] = size(data);

% We can never have more labels than channels
nYLabels = min(nYLabels,nChannels);

if exist('xAxis','var')
    % Ensure that x is a row vector
    xAxis = xAxis(:)';
    strX  = 'Time (s)';
else
    xAxis = 1:nSamples;
    strX  = 'Number of samples';
end

if exist('yAxis','var')
    % Ensure that x is a row vector
    yAxis = yAxis(:)';
    strY  = 'Center frequency (Hz)';
else
    yAxis = 1:nChannels;
    strY  = 'Number of channels';
end

% Check for proper dimensions
if any(nSamples ~= length(xAxis) || nChannels ~= length(yAxis))
    error(['The length of "xAxis" and "yAxis" must correspond to the ',...
           'first ("nSamples") and second ("nChannels") dimension of "data".'])
end


%% *****************************  PLOT DATA  ******************************
% 
% 
% Normalization
if bNorm
    % Normalize channels
    data = normalizeData(data,'meanvar');
end

% Get global maxima
m = max(max(data));
n = max(5 * zoom * data(:,nChannels)/m+nChannels * 10);

hold on;

% Loop over number of channels
for ii = nChannels : -1 : 1
   plot(xAxis,5 * zoom * data(:,ii)'/m+ii * 10,'k')
end

hold off;
 
% Find the spacing for the y-axis which evenly divides the y-axis
set(gca,'ytick',linspace(1,nChannels,nYLabels) * 10);
if nChannels > 1
    set(gca,'yticklabel',round(interp1(1:nChannels,yAxis,linspace(1,nChannels,nYLabels))));
else
    set(gca,'yticklabel',round(yAxis));
end

xlabel(strX)
ylabel(strY)

% Set axis limits
xlim([xAxis(1) xAxis(end)]);

if nChannels > 1
    ylim([0 n]);
end