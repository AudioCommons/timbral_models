function data = data_goode1994(varargin)
%DATA_GOODE1994 Returns data points from the Goode et al. (1994) paper
%   Usage: data = data_goode1994_data(flag)
%
%   DATA_GOODE1994(flag) returns various data from Goode et al. (1994). 
%   (See the particular flags for more details on the nature of the data).
%
%   The flag may be one of:
%
%     'noplot'    Don't plot, only return data. This is the default.
%
%     'plot'      Plot the data.
%
%     'fig1_104'  (default) Returns the stapes footplate diplacement at an SPL of 104 dB 
%                 as shown in Fig. 1. (Currently the only option).
%
%   Examples:
%   ---------
% 
%   To plot the Fig. 1 from Goode et al. (1994) use :
%
%     data_goode1994('fig1_104','plot');
%
%   References:
%     R. Goode, M. Killion, K. Nakamura, and S. Nishihara. New knowledge
%     about the function of the human middle ear: development of an improved
%     analog model. The American journal of otology, 15(2):145-154, 1994.
%     
%
%   Url: http://amtoolbox.sourceforge.net/doc/humandata/data_goode1994.php

% Copyright (C) 2009-2014 Peter L. Søndergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.6
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%   AUTHOR: Peter L. Søndergaard


%% ------ Check input options --------------------------------------------

% Define input flags
definput.flags.type={'fig1_104'};
definput.flags.plot = {'noplot','plot'};

% Parse input options
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

%% ------ Data points from the paper ------------------------------------
%
% Data for the given figure

if flags.do_fig1_104
  % data derived from Goode et al. 1994, Figure 1
  data = [...
    400	0.19953; ...
    600	0.22909; ...
    800	0.21878; ...
    1000 0.15136; ...
    1200 0.10000; ...
    1400 0.07943; ...
    1600 0.05754; ...
    1800 0.04365; ...
    2000 0.03311; ...
    2200 0.02754; ...
    2400 0.02188; ...
    2600 0.01820; ...
    2800 0.01445; ...
    3000 0.01259; ...
    3500 0.00900; ...
    4000 0.00700; ...
    4500 0.00457; ...
    5000 0.00500; ...
    5500 0.00400; ...
    6000 0.00300; ...
    6500 0.00275];

  if flags.do_plot
    figure;
    loglog(data(:,1)/1000,data(:,2),'+');
    xlabel('Frequency (kHz)');
    ylabel('Peak-to-peak displacement ({\mu}m)');
    %axis([0.1,10,min(data(:,2))*1e6,max(data(:,2))*1e6]);
    axis([0.1,10,0.001,10]);
  end;
end;


