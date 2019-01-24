% erbtofreq.m - converts erbscale to frequencyscale
%               Uses a scaling based on the equivalent rectangular bandwidth
%               (ERB) of an auditory filter at centre frequency fc:
%               ERB(fc) = 24.7 + fc[Hz]/9.265 (Glasberg and Moore, JASA 1990). 
%
% Usage: freq = erbtofreq(erb)
%
% erb  = erbscaled input vector  
%
% freq = output vector in Hz

% Author: Stephan Ewert, Universitaet Oldenburg.
% $Revision: 1 $  $Date: 1999/11/18 12:01:37 $

function freq = erb2freq(erb)

freq = 24.7*9.265*(exp(erb/9.265)-1);

