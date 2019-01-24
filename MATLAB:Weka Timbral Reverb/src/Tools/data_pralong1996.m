function data = data_pralong1996(varargin)
%DATA_PRALONG1996 Head phone data from Pralong & Carlile (1996)
%   Usage: data = data_pralong1996(flag)
%
%   DATA_PRALONG1996(flag) returns data points from the Pralong & Carlile
%   (1996) paper.
%
%   The flag may be one of:
%
%     'noplot'  Don't plot, only return data. This is the default.
%
%     'plot'    Plot the data.
%  
%     'fig1e'   Data from Fig. 1(e), Gain of Sennheiser 250 Linear
%               circumaural headphones. This is the default
%
%   Examples:
%   ---------
%
%   Figure 1e can be displayed using:
%
%     data_pralong1996('plot');
%
%   References:
%     D. Pralong and S. Carlile. The role of individualized headphone
%     calibration for the generation of high fidelity virtual auditory space.
%     J. Acoust. Soc. Am., 100:3785-3793, 1996.
%     
%
%   Url: http://amtoolbox.sourceforge.net/doc/humandata/data_pralong1996.php

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
  
% Define input flags
definput.flags.type={'fig1e'};
definput.flags.plot = {'noplot','plot'};

[flags,keyvals]  = ltfatarghelper({},definput,varargin);

if flags.do_fig1e
  data = [...
      125,	1; ...
      250,	1; ...
      500,	1; ...
      1000,	0.994850557; ...
      1237.384651,	0.994850557; ...
      1531.120775,	0.994850557; ...
      1894.585346,	1.114513162; ...
      2002.467159,	1.235743262; ...
      2344.330828,	1.867671314; ...
      2721.273584,	2.822751493; ...
      3001.403462,	2.180544843; ...
      3589.453635,	1.442755787; ...
      4001.342781,	1.173563859; ...
      4441.534834,	1.37016005; ...
      5004.212211,	1.599690164; ...
      5495.887031,	1.37016005; ...
      5997.423738,	1.114513162; ...
      6800.526258,	0.648125625; ...
      6946.931144,	0.631609176; ...
      7995.508928,	0.276505667; ...
      8414.866811,	0.084335217; ...
      9008.422743,	0.084335217; ...
                 ];
  
  
  if flags.do_plot
    figure;
    x=data(:,1);
    freqresp=20*log10(data(:,2));
    semilogx(x,freqresp);
    xlabel('Frequency (Hz)');
    ylabel('Gain (dB)');
  end;
  
end;

