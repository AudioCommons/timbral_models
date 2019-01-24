function [out,states] = leakyIntegrator(in,fs,decaySec,states)
%leakyIntegrator   Perform leaky integration.
%
%USAGE
%            out = leakyIntegrator(in,fs)
%   [out,states] = leakyIntegrator(in,fs,decaySec,states)
%
%INPUT ARGUMENTS
%             in : input data [nSampels x nChannels x ... ] The integration
%                  is performed across the first dimension. 
%             fs : sampling frequency of input data in Hertz
%       decaySec : time constant of leaky integrator in seconds 
%                  (default, decaySec = 8E-3 )
%         states : filter states [1 x nChannels x ... ] 
%                  (default, states = zeros(1,nChannels,... ))
%
%OUTPUT ARGUMENTS
%            out : output data [nSamples x nChannels x ... ]
%         states : integrator filter states [1 x nChannels x ... ]

%   Developed with Matlab 8.2.0.701 (R2013b). Please send bug reports to:
%   
%   Authors :  Tobias May © 2014
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
% 
%   History :  
%   v.0.1   2014/03/08
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
if nargin < 3 || isempty(decaySec); decaySec = 8E-3; end


%% CHECK FILTER STATES
% 
% 
% Determine dimension of input
dim = size(in);

% Check filter states
if exist('states','var') && ~isequal(size(states),dim(2:end))
    error('Filter states dimension mismatch!')
else
    % Initialize states with zeros
    states = zeros([1 dim(2:end)]); 
end


%% APPLY FILTER
% 
% 
% Filter decay
intDecay = exp(-(1/(fs * decaySec)));

% Integration gain
intGain = 1 - intDecay;

% Apply integration filter
[out, states] = filter(intGain, [1 -intDecay], in, states);