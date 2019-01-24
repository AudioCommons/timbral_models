function bStable = isStable(b,a)
%isStable   Determine stability of IIR filter.
%   The IIR filter, defined by the coefficients b and a, is considered to 
%   be stable if all poles are within the unit circle.
% 
%USAGE 
%   bStable = isStable(b,a);
% 
%INPUT ARGUMENTS
%   b : filter coefficients of the numerator
%   a : filter coefficients of the denominator
% 
%OUTPUT ARGUMENTS
%   bStable : 1 if the IIR filter is stable, and 0 otherwise
%
%EXAMPLE
%   % Whitening filter coefficients
%   b = 1; a = [1 -0.97];
% 
%   % Check stability of IIR filter
%   isStable(b,a)


%% CHECK INPUT ARGUMENTS
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 2
    help(mfilename);
    error('Wrong number of input arguments!')
end


%% DETERMINE STABILITY OF IIR FILTERS
% 
% 
% Get zero (z) and pole (p) locations of the filter
[z,p] = tf2zp(b(:).',a(:).'); %#ok

% Determine if all poles are less than unit magnitude (within unit circle)
bStable = all(abs(p) < 1);


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