function out = frameData(input,blockSize,stepSize,win,bZeroPad)
%frameData   Segment the input signal into overlapping frames.
% 
%USAGE
%      frames = frameData(input,blockSize,stepSize,win,bZeroPad);
% 
%INPUT ARGUMENTS
%       input : input signal [nSamples x 1]
%   blockSize : frame size in samples
%    stepSize : step size of frames in samples
%         win : string or vector defining analysis window
%    bZeroPad : Zero-pad input signal to fill last frame 
%               (default, bZeroPad = false)
% 
%OUTPUT ARGUMENTS
%      frames : frame-based signal [blockSize x nFrames]

%   Developed with Matlab 7.9.0.529 (R2009b). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2009 
%              TUe Eindhoven and Philips Research  
%              t.may@tue.nl      tobias.may@philips.com
%
%   History :
%   v.0.1   2009/11/23
%   v.0.2   2010/05/18 automatically determine the number of frames
%   v.0.3   2013/08/19 added zero-padding
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
if nargin < 4 || nargin > 5
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default parameter
if nargin < 5 || isempty(bZeroPad); bZeroPad = false; end
    
% Determine size of input signal
[nSamples,nChannels] = size(input);

% Check for proper size
if nChannels > 1     
    error('Monaural input is required!')
end


%% ****************************  ZERO-PADDING  ****************************
% 
% 
% Compute number of frames
nFrames = (nSamples-(blockSize-stepSize))/(stepSize);

% Append zeros
if bZeroPad
    % Number of frames (at least one)
    nFrames = ceil(max(1,nFrames));
    
    % Compute number of required zeros
    nZeros = (stepSize * nFrames + (blockSize-stepSize)) - nSamples;
    
    % Pad zeros
    input = [input; zeros(nZeros,1)];
else
    % No zero padding, exclude elementes that do not fit into last frame
    nFrames = max(0,floor(nFrames));
end
 
% Check if window is a window in samples or a string
if ischar(win)
    win = window(win,blockSize);
else
    if blockSize ~= length(win)    
       error('Mismatch between blocksize and window size.') 
    end
end
        
    
%% ***********************  FRAME-BASED PROCESSING  ***********************
% 
% 
% Framing (MEX processing)
out = frameDataMEX(input,blockSize,stepSize,win,nFrames);


%   ***********************************************************************
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   ***********************************************************************