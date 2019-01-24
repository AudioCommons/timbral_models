function out = applyCenterClipping(frames,method,alpha)
% 
%USAGE
%    frames = applyCenterClipping(frames,method,alpha)
%
%INPUT PARAMETERS
% 
%OUTPUT PARAMETERS
% 

%   Developed with Matlab 8.2.0.701 (R2013b). Please send bug reports to:
%   
%   Authors :  Tobias May © 2014
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
% 
%   History :  
%   v.0.1   2014/03/06
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS 
% 
% 
% Check for proper input arguments
if nargin ~= 3
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Check valid range
if ~isfinite(alpha) || alpha > 1 || alpha < 0
    error('Clipping level ''alpha'' must be within [0, 1].')
end

% Determine input size
[nSamples,nFrames] = size(frames);

 % Allocate memory
out = zeros(nSamples,nFrames);


%% COMPUTE FRAME-BASED MAXIMA AND MINIMA
% 
% 
% Find frame-based threshold values
maxVal = repmat(abs(max(frames,[],1)) * alpha,[nSamples 1]);
minVal = repmat(abs(min(frames,[],1)) * alpha,[nSamples 1]);

% Find all samples greater and smaller than the predefined threshold
bAbove = frames >=  maxVal;
bBelow = frames <= -minVal;


%% PERFORM CENTER CLIPPING
% 
% 
% Apply center clipping method
switch lower(method)
    case 'clc'
        % Clip and compress
        out(bAbove) = frames(bAbove) - maxVal(bAbove);
        out(bBelow) = frames(bBelow) + minVal(bBelow);
    case 'clp'
        % Clip
        out(bAbove) = frames(bAbove);
        out(bBelow) = frames(bBelow);        
    case  'sgn'
        % Sign
        out(bAbove) =  1;
        out(bBelow) = -1;                
    otherwise
        error('Center clipping method ''%s'' is not supported.',lower(method));
end


