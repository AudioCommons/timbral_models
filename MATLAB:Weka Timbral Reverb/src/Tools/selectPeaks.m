function [sepPeakIdx,sepPeakVal] = selectPeaks(peakIdx,peakVal,minDist)

%   Developed with Matlab 7.7.0.471 (R2008b). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2009
%              TUe Eindhoven and Philips Research  
%              t.may@tue.nl      tobias.may@philips.com
%
%   History :   
%   v.1.0   2009/02/17
%   ***********************************************************************


%% ***********************  CHECK INPUT ARGUMENTS  ************************
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 3
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 3 || isempty(minDist); minDist = 1; end


%% **************************  PEAK SELECTION  ****************************
% 
% 
% Scan all peaks from highest to lowest and retain only those separated by
% at least minDist samples.  
if minDist > 2 && ~isempty(peakIdx)
    % Sort peaks
    [peakVal,sortIdx] = sort(peakVal,'descend');

    % Initialize first peak
    sepPeakVal = peakVal(1); 
    sepPeakIdx = peakIdx(sortIdx(1));
    
    % select highest peak first
    for ii = 2 : length(peakVal) 
        % Analyze remaining peaks 
        [np, junk] = size(sepPeakVal); 
        % Reset counter
        count = 0; 
        % Loop over number of remaining peaks
        for jj = 1 : np
            % Check distance
            if(abs(peakIdx(sortIdx(ii)) - sepPeakIdx(jj)) < minDist) 
                % Increase counter
                count = count + 1;
                break;
            end
        end
        if (count == 0) 
            sepPeakVal    = [sepPeakVal;    peakVal(ii)         ]; 
            sepPeakIdx = [sepPeakIdx; peakIdx(sortIdx(ii))];
        end
    end
else
    sepPeakIdx = peakIdx;
    sepPeakVal = peakVal;
end

% Re-arrange
[sepPeakIdx,sortIdx] = sort(sepPeakIdx,'ascend');
sepPeakVal           = sepPeakVal(sortIdx);

