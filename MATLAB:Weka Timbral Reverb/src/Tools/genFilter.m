function [b,a] = genFilter(type,fsHz)
%genFilter   Create various types of filters. 
%
%USAGE
%   [b,a] = genFilter(type,fs)
%
%INPUT ARGUMENTS
%   type : filter type
%   fsHz : sampling frequency in Hertz
% 
%OUTPUT ARGUMENTS
%   b : filter coefficients (numerator)
%   a : filter coefficients (denominator)
% 
%REFERENCES
%   [1] 
% 

%   Developed with Matlab 8.2.0.701 (R2013b). Please send bug reports to:
%   
%   Authors :  Tobias May © 2014
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
% 
%   History :  
%   v.0.1   2014/03/07
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS 
% 
% 
% Check for proper input arguments
if nargin ~= 2
    help(mfilename);
    error('Wrong number of input arguments!')
end


%% CREATE FILTER COEFFICIENTS
% 
%
% Select filter type
switch lower(type)
    case 'middleear_meddisomard_1997'
        
        % Design second-order bandpass
        [b, a] = butter(2,[450 8500] / (fsHz / 2));
        
    case 'removedc'
        
        % Remove DC component
        [b, a] = butter(4, 50 / (fsHz / 2), 'high');
        
    otherwise
        error('Filter type ''%s'' is not supported',lower(type));
end