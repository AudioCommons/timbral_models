function arrange(fig1,fig2)
%arrange   Arrange two figures on the screen. 
%   

%   Developed with Matlab 8.2.0.701 (R2013b). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2013
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
% 
%   History :  
%   v.0.1   2014/01/22
%   ***********************************************************************

% Check for proper input arguments
if nargin ~= 2
    help(mfilename);
    error('Wrong number of input arguments!')
end

set(0,'Units','pixels');
scnsize = get(0,'ScreenSize');
    
position = get(fig1,'Position');
outerpos = get(fig1,'OuterPosition');
borders  = outerpos - position;

edge = -borders(1)/2;
pos1 = [edge,...
        scnsize(4) * (1/2),...
        scnsize(3)/2 - edge,...
        scnsize(4)/2];

pos2 = [scnsize(3)/2 + edge,...
        pos1(2),...
        pos1(3),...
        pos1(4)];

set(fig1,'OuterPosition',pos1)
set(fig2,'OuterPosition',pos2)
    

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