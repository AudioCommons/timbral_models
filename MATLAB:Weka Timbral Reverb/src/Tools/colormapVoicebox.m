function cmap = colormapVoicebox(nColors,bInvert)
%colormapVoicebox   Black to green color map with good grayscale linearity.
%
%USAGE 
%      cmap = colormapVoicebox
%      cmap = colormapVoicebox(nColors,bInvert)
%
%INPUT ARGUMENTS
%   nColors : number of distinct colors (default, nColors = 64)
%   bInvert : invert color map (default, bInvert = true)
% 
%OUTPUT ARGUMENTS
%      cmap : color map [nCols x 3]
% 
%EXAMPLE
%   % Plot a 2D image
%   figure;imagesc(peaks(30));
%
%   % Apply new colormap
%   colormap(colormapVoicebox)
% 
%   See also colormapParula.


%   Developed with Matlab 8.3.0.532 (R2014a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2014
%              Technical University of Denmark (DTU)
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2014/09/26
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS
% 
% 
% Check for proper input arguments
if nargin > 2 || nargout == 0
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values 
if nargin < 1 || isempty(nColors); nColors = 64;    end
if nargin < 2 || isempty(bInvert); bInvert = true;  end

% Color map taken from spgrambw (Voicebox toolbox)
cmap = [0 0 0; 7 0 17; 14 0 33; 21 0 50; 29 0 67; 36 0 84; 43 0 100; 
    50 0 117; 57 0 134; 64 0 150; 72 0 167; 80 3 164; 89 7 156; 97 11 149; 
    106 15 142; 114 19 134; 123 23 127; 131 27 119; 140 31 112; 149 35 105; 
    157 39 97; 166 43 90; 174 47 82; 183 51 75; 192 55 68; 200 59 60; 
    209 63 53; 217 67 45; 226 71 38; 234 75 31; 243 79 23; 252 83 16; 
    255 88 12; 255 95 12; 255 102 11; 255 109 11; 255 116 10; 255 123 10; 
    255 130 9; 255 137 9; 255 144 8; 255 151 8; 255 158 7; 255 165 7; 
    255 172 6; 255 179 6; 255 186 5; 255 193 4; 255 200 4; 255 207 3; 
    255 214 3; 255 221 2; 255 228 2; 255 235 1; 255 242 1; 255 249 0; 
    255 252 22; 255 252 55; 255 253 88; 255 253 122; 255 254 155; 
    255 254 188; 255 255 222; 255 255 255]/255;

% Use linear interpolation
cmap = interp1((1:size(cmap,1))/size(cmap,1),cmap,(1:nColors)/nColors);

% Invert scale and start with the darkest color
if bInvert
    cmap = flipud(1-cmap);
end