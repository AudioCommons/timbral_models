function [xcorr,lags] = calcXCorr(sig1,sig2,maxLags,scale)
%calcXCorr   FFT-based cross-correlation function.
%
%USAGE
%	        LAGS = calcXCorr(SIG1,SIG2)
%	[XCORR,LAGS] = calcXCorr(SIG1,SIG2,MAXLAGS,SCALE)
%
%INPUT ARGUMENTS
%      SIG1 : first input signal, either vector or matrix of dimension
%             [nSamples x nChannels]
%      SIG2 : second input signal either vector or matrix of dimension
%             [nSamples x nChannels]. 
%             If input are matrices, the cross-correlation function is
%             computed over the columns. 
%   MAXLAGS : the correlation sequenze is computed over the lag range
%             -MAXLAGS:MAXLAGS 
%             (default, MAXLAGS = max(size(SIG1,1),size(SIG2,1))-1)
%     SCALE : string specifying normalization 
%             'biased'   - scale XCORR by 1/nSamples
%             'unbiased' - scale XCORR by 1/(nSamples-abs(lags))
%             'coeff'    - scale XCORR by the autocorrelation at lag zero
%             'none'     - no scaling
%             (default, SCALE = 'none')
%
%OUTPUT ARGUMENTS
%     XCORR : cross-correlation function between SIG1 and SIG2
%             [2*MAXLAGS+1 x nChannels]
%      LAGS : time lags of cross-correlation function [-MAXLAG:MAXLAG x 1]


%   Developed with Matlab 8.2.0.701 (R2013b). Please send bug reports to:
%   
%   Author  :  Tobias May © 2014
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
% 
%   History :   
%   v.0.1   2014/03/04
%   v.0.2   2014/03/05 fixed zero-padding
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS 
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 4
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 4 || isempty(scale); scale = 'none'; end

% FFT approach for matrices...
[blockSize1,nChannels1,check1] = size(sig1);
[blockSize2,nChannels2,check2] = size(sig2);

% Check if number of channels are identical
if ~isequal(nChannels1,nChannels2)
    error('Number of channels must be identical for both input matrices!');
end

% Check dimensionaly of input signals
if any([check1 check2] > 1)
    error('Require either vector or a matrix as input!');
end

% Work along non-singleton dimension
if (blockSize1 == 1 && blockSize2 == 1)
    % Transpose data
    sig1 = transpose(sig1);
    sig2 = transpose(sig2);
    
    % Swap channel to blockSize
    blockSize1 = nChannels1;
    blockSize2 = nChannels2;  
end

% Determine maximum number of samples ...
M = max(blockSize1,blockSize2);

% Set default number of lags if necessary...
if nargin < 3 || isempty(maxLags);  maxLags = M - 1; end


%% COMPUTE SPECTRA
% 
% 
% Transform both input signals
X = fft(sig1,2^nextpow2(2*M-1));
Y = fft(sig2,2^nextpow2(2*M-1));

% Compute cross-power spectrum
XY = X.*conj(Y);


%% BACK TO TIME DOMAIN
% 
% 
% Inverse FFT
c = real(ifft(XY));   

% Define lags 
lags = (-maxLags:maxLags).';

% Keep only the lags we want and move negative lags in front 
if maxLags >= M,
    % Pad with zeros
    padZeros = zeros(maxLags-M+1,nChannels1);
	c        = [padZeros; c(end-M+2:end,:);c(1:M,:); padZeros];
else
	c = [c(end-maxLags+1:end,:);c(1:maxLags+1,:)];
end


%% NORMALIZATION
% 
% 
% Normalization
switch lower(scale)
    case 'none'
        xcorr = c;
    case 'biased'
        xcorr = c / M;
    case 'unbiased'
        scale = M-abs(lags'); scale(scale<=0)=1;
        xcorr = c./repmat(scale(:),[1 nChannels1]);
    case 'coeff'
        % Compute autocorrelation at lag zero
        powL = sum(sig1.^2,1);
        powR = sum(sig2.^2,1);
        
        % Normalization
        xcorr = c ./ repmat(eps + sqrt(powL .* powR),[length(lags) 1]);
    otherwise
        error(['Normalization method ''',lower(scale),...
               ''' is not supported.'])
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