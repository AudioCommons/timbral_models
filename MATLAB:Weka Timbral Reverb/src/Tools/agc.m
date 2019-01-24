function out = agc(in,fsHz,tauSec,bMultiChannel,bPlot)
%agc   Automatic gain control.
%
%USAGE 
%             out = agc(in,fsHz,tauSec,bMultiChannel)
%
%INPUT ARGUMENTS
%              in : input signal arranged as [nSamples x nChannels]
%        stepSize : sampling frequency in Hertz 
%          tauSec : time constant of RMS estimation 
%   bMultiChannel : If true, multi-channel differences are preserved
%                   (default, bMultiChannel = true)
%           bPlot : plot AGC processing (default, bPlot = false)
% 
%OUTPUT ARGUMENTS
%             out : normalized signal [nSamples x nChannels]


%   Developed with Matlab 8.2.0.701 (R2013b). Please send bug reports to:
%   
%   Author  :  Tobias May © 2014
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2014/10/23
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS 
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 5
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 3 || isempty(tauSec);        tauSec        = 2;     end
if nargin < 4 || isempty(bMultiChannel); bMultiChannel = true;  end
if nargin < 5 || isempty(bPlot);         bPlot         = false; end

% Prevent division overflow
epsilon = 1E-8;


%% APPLY AUTOMATIC GAIN CONTROL
%
%
% Time constant
tau = tauSec * fsHz;
ax2 = [1 -exp(-1/tau)];
bx2 = sum(ax2);

% Initialize filter using an average of 1 tau
if bMultiChannel
    sm = repmat(max(mean(in(1:min(size(in,1),round(tau)),:).^2)),[1 size(in,2)]);
else
    sm = mean(in(1:min(size(in,1),round(tau)),:).^2);
end

% Estimate normalization constant
normFactor = sqrt(filter(bx2,ax2,in.^2,-ax2(2)*sm)) + epsilon;

% Preserve multi-channel differences
if bMultiChannel
    % Use maximum normalization constant across all channels
    normFactor = repmat(max(normFactor,[],2),[1 size(in,2)]);
end

% Check normalization constant
if any(normFactor(:)==0) || any(~isfinite(normFactor(:)))
    error('Normalization constant is zero/not finite')
else
    % Normalize signal to have an estiamted RMS of one
    out = in./normFactor;
end


%% DISPLAY EFFECT OF AGC
%
%
% Plot result
if nargout == 0 || bPlot == true
    figure;
    ax(1) = subplot(2,1,1);
    plot(in)
    title('Before AGC')
    ax(2) = subplot(2,1,2);
    plot(out)
    title('After AGC')
    
    axis tight;
    linkaxes(ax,'x');
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