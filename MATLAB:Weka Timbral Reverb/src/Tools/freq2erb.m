% freqtoerb.m - converts to frequencyscale erbscale.
%               Uses a scaling based on the equivalent rectangular bandwidth
%               (ERB) of an auditory filter at centre frequency fc:
%               ERB(fc) = 24.7 + fc[Hz]/9.265 (Glasberg and Moore, JASA 1990). 
%
% Usage: erb = freqtoerb(freq)
%
% freq = input vector in Hz
%
% erb  = erbscaled output vector

% Author: Stephan Ewert, Universitaet Oldenburg.
% $Revision: 1 $  $Date: 1999/11/18 12:01:37 $

function erb = freq2erb(freq)

erb = 9.265*log(1+freq./(24.7*9.265));

