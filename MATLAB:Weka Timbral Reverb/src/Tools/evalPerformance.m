function [pcCorrect,rmse] = evalPerformance(azRef,azEst,thresDeg)
%evalPerformance   Evaluate sound source localization performance.
%
%USAGE 
%   [pcCorrect,rmse] = evalPerformance(azRef,azEst,thresDeg)
%
%INPUT ARGUMENTS
%       azRef : vector with true sound source positions      [K x 1]
%       azEst : vector with estimated sound source positions [N x 1]
%    thresDeg : absolute error boundary in degree (default, thresDeg = 10)
%
%OUTPUT ARGUMENTS
%   pcCorrect : Percentage of correctly localized sound source positions
%        rmse : Root mean square error of correctly localized positions

%   Developed with Matlab 8.2.0.701 (R2013b). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2013
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
% 
%   History :  
%   v.0.1   2014/01/26
%   v.0.2   2014/01/31 check if azRef and azEst are vectors
%   ***********************************************************************


%% 1. CHECK INPUT ARGUMENTS  
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 3
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default parameter
if nargin < 3 || isempty(thresDeg); thresDeg = 10; end


%% 2. EVALUATE LOCALIZATION PERFORMANCE
% 
% 
if any(size(azRef)<2)
    % Ensure column vectors
    azRef = azRef(:);
    azEst = azEst(:);
else
   error('''azRef'' and ''azEst'' must be vectors.') 
end

% Calculate localization error in degree
azErrorDeg = calcLocalizationError(azRef,azEst);

% Percentage of correctly localized sound sources wihtin "thresDeg"
bCorrect   = azErrorDeg <= thresDeg;
pcCorrect  = 100 * sum(bCorrect)/numel(azRef); 

% RMSE of correctly localized sound source positions
rmse       = sqrt(mean(azErrorDeg(bCorrect).^2));


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