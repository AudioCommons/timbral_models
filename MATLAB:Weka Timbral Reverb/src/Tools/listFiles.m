%LISTFILES  List all files of directory and its sub-directories.
%   LISTFILES('a_directory') lists the files in a directory and its
%   sub-directories up to a depth of four sub-directories. Pathnames
%   and wildcards may be used.
%   For example, LISTFILES('a_directory', '*.m') lists all the M-files
%   in a directory and its sub-directories up to a depth of four
%   sub-directories.
%   LISTFILES('a_directory', '*.m', 6) lists the files of a directory
%   and its sub-directories up to a depth of six sub-directories. Whereas
%   LISTFILES('a_directory', '*.m',-1) lists the files of all existing
%   sub-directories with infinite recursion. The default recursion depth
%   is to list files of up to four sub-directories.
%
%   D = LISTFILES('a_directory') returns the results in an M-by-1
%   structure with the fields: 
%       name  -- filename (incl. the 'a_directory' path)
%       date  -- modification date
%       bytes -- number of bytes allocated to the file
%       isdir -- 1 if name is a directory and 0 if not
%
%   Hint:
%   To convert the struct array D into a cell-array C, which contains only
%   filenames (one per row), you may do so by typing: C = {D.name}'
%
%   See also DIR, LISTDIRS.

% Auth: Sven Fischer
% Vers: v0.821
function [ stFileList ] = listFiles(szCurDir, szFileMask, iRecursionDepth)

%-------------------------------------------------------------------------%
% Check input arguments.
narginchk(0,3);
% Check output arguments.
nargoutchk(0,1);
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Check function arguments and set default values, if necessary.
%-------------------------------------------------------------------------%
if( (nargin<1) || (isempty(szCurDir       )) ), szCurDir        =  ''; end
if( (nargin<2) || (isempty(szFileMask     )) ), szFileMask      = '*'; end
if( (nargin<3) || (isempty(iRecursionDepth)) ), iRecursionDepth =   4; end
%-------------------------------------------------------------------------%

stFileList = dir( fullfile( szCurDir, szFileMask ) );
stFileList = stFileList(  ~[stFileList.isdir]  );
for( k = [ 1 : length(stFileList) ] )
  stFileList(k).name = fullfile( szCurDir, stFileList(k).name );
end

% If we have to process sub-directories recursively...
if( (iRecursionDepth > 0) || (iRecursionDepth == -1) )
    % Decrease recursion counter by one, if not set to infinite...
    if( iRecursionDepth > 0 ), iRecursionDepth = iRecursionDepth - 1; end

    % Get a list of all existing sub-directories (exclusive '.' and '..').
    stSubDirs = dir( szCurDir );
    stSubDirs = stSubDirs(  [stSubDirs.isdir]  );
    if( ~isempty(stSubDirs) )
        if( strcmp(stSubDirs(1).name,  '.') ), stSubDirs(1) = []; end
        if( strcmp(stSubDirs(1).name, '..') ), stSubDirs(1) = []; end

        % Process all subdirectories and append all each file list to the
        % list created above.
        for( k = [ 1 : length(stSubDirs) ] )
            szSubDir = fullfile( szCurDir, stSubDirs(k).name);
            stFileList = [ stFileList; ...
                listFiles( szSubDir, szFileMask, iRecursionDepth) ];
        end
    end
end

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