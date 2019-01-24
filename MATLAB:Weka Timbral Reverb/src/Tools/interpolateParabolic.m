function [delta,rHat] = interpolateParabolic(r, maxIdx)
%interpolateParabolic   Multi-channel parabolic interpolation.
%
%USAGE
% [DELTA,MAXVAL] = interpolateParabolic(R,MAXIDX);
%
%INPUT ARGUMENTS
%        R : detection function [nSamples x nChan]
%   MAXIDX : initial maximum estimate (maximum argument(index) of R). In
%            case R is single channel data, interpolation can be performed
%            at multiple indices -> MAXIDX [nIdx x 1]
%            
%            Otherwise, if nChannels > 1, interpolation is performed at one
%            index per channel. So the dimension of MAXIDX must correspond
%            to the number of channels in R -> MAXIDX [nChan x 1]
%         
%OUTPUT ARGUMENTS
%    DELTA : fractional part of maximum index [nIdx x 1]  |  [nChan x 1]
%     RHAT : amplitude of estimated maximum   [nIdx x 1]  |  [nChan x 1]        
%
%   See also interpolateExponential and interpolateLagrange.


%   Developed with Matlab 7.4.0.287 (R2007a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2008-2009
%              TUe Eindhoven and Philips Research  
%              t.may@tue.nl      tobias.may@philips.com
%
%   History :   
%   v.0.1   2008/02/17
%   v.0.2   2009/01/18 allow multiple maxima for single channel input
%   v.0.3   2009/04/15 return integer value of rHat for invalid indices 
%   ***********************************************************************

% Check for proper input arguments
if nargin ~=2
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Check dimension
if any(size(r) == 1)
    % Ensure column vector
    r = r(:);
end

% Ensure that t is a column vector
t = maxIdx(:);

% Determine size of r
[blockSize,nFrames] = size(r);

% Single/multi-channel modus
if isequal(nFrames,1)
    % Process multiple maxima for single channel
    delta = zeros(size(t,1),1);
    % Set index shift to zero (only used for multi-channel input)
    idxShift = delta;
else
    % Check proper dimension
    if size(t,1) ~= nFrames
        error(['The second dimension of r must be identical to'...
               ' the dimension of t.'])
    end
    % Allocate memory
    delta = zeros(nFrames,1);
    % Offset to work vectorized with matrix indices 
    idxShift = (0:nFrames-1)' * blockSize;
end

% Check for valid indices and do not interpolate if the detected maximum is
% located at the edge of r
validIdx = t > 1 & t < blockSize;

% Three points around the maximum of r
rM1 = r(idxShift(validIdx) + t(validIdx) - 1);
rP0 = r(idxShift(validIdx) + t(validIdx) + 0);
rP1 = r(idxShift(validIdx) + t(validIdx) + 1);

% Parabolic interpolation to estimate fractional part of the maximum
delta(validIdx) = 0.5 * (rM1-rP1) ./ (rM1 - 2 * rP0 + rP1);

% Estimate amplitude
if nargout > 1
    % Allocate memory
    rHat = zeros(size(delta));
    % Estimate height of peak 
    rHat(validIdx) = rP0 - 0.25 * (rM1-rP1) .* delta(validIdx);
    % Use integer peak position for not valid indices (e.g. endpoints)
    rHat(~validIdx) = r(idxShift(~validIdx) + t(~validIdx));
end


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