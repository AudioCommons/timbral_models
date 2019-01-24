function varargout = frameDataMEX(varargin)
%frameDataMEX   Private method, which is normally shadowed by the
%   corresponding MEX routine. It will ONLY be executed, if the underlying
%   MEX routine is not available.  

%   Developed with Matlab 7.9.0.529 (R2009b). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2009 
%              TUe Eindhoven and Philips Research  
%              t.may@tue.nl      tobias.may@philips.com
%
%   History :
%   v.0.1   2009/11/23
%   ***********************************************************************

% Report error message
error(['MEX function "',mfilename,'" is not available for your ',...
       'operating system (',computer,'). Run "mex" to built ',...
       'required MEX binaries.'])

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