function data = data_lopezpoveda2001(varargin)
%DATA_LOPEZPOVEDA2001  Data from Lopez-Poveda & Meddis (2001)
%   Usage: data = data_lopezpoveda2001(flag)
%
%   DATA_LOPEZPOVEDA2001(flag) returns data points from the paper by 
%   Lopez-Poveda and Meddis (2001).
%
%   The flag may be one of:
%
%     'noplot'       Don't plot, only return data. This is the default.
%
%     'plot'         Plot the data.
%  
%     'fig2a'        Data from Fig. 2(a), outer ear filter.
%
%     'fig2b'        Data from Fig. 2(b), middle ear filter.
%
%   For Fig. 2b you can choose between:
%
%     'goode'        Return the data points derived from Goode et al. (1994)
%                    This is the default.
%
%     'lopezpoveda'  Return the data points just read from Fig. 2b
%                    of Lopez-Poveda and Meddis (2001)
%
%   Examples:
%   ---------
%
%   To display Figure 2a, use:
%
%     data_lopezpoveda2001('fig2a','plot');
%
%   To display Figure 2b, use:
%
%     data_lopezpoveda2001('fig2b','plot');
%
%   References:
%     E. Lopez-Poveda and R. Meddis. A human nonlinear cochlear filterbank.
%     J. Acoust. Soc. Am., 110:3107-3118, 2001.
%     
%     
%
%   Url: http://amtoolbox.sourceforge.net/doc/humandata/data_lopezpoveda2001.php

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

%   AUTHOR: Peter L. Søndergaard, Katharina Egger


%% ------ Check input options --------------------------------------------

% Define input flags
definput.flags.plot = {'noplot','plot'};
definput.flags.type = {'missingflag','fig2a','fig2b'};
definput.flags.datapoints = {'goode','lopezpoveda'};

% Parse input options
[flags,keyvals]  = ltfatarghelper({},definput,varargin);

if flags.do_missingflag
  flagnames=[sprintf('%s, ',definput.flags.type{2:end-2}),...
             sprintf('%s or %s',definput.flags.type{end-1},definput.flags.type{end})];
  error('%s: You must specify one of the following flags: %s.',upper(mfilename),flagnames);
end;

%% ------ Data points from the paper ------------------------------------
%
% Data for the given figure

if flags.do_fig2a

  data=data_pralong1996;
  
  if flags.do_plot

    fs=22050;
    bout=headphonefilter(fs);

    % Manually calculate the frequency response.
    fout = 20*log10(abs(fftreal(bout)));

    % Half the filter length.
    n2=length(fout);

    figure;
    hold on;
    % Plot the measured data
    x=data(:,1);
    freqresp=20*log10(data(:,2));
    semilogx(x,freqresp,'ro');
    
    % Plot the filter
    x_filter=linspace(0,fs/2,n2);
    semilogx(x_filter,fout);
    hold off;
  end;
  
end;

if flags.do_fig2b
  
  if flags.do_goode  
      % get peak-to-peak displacement data
      data = data_goode1994;

      % get velocity (proportional voltage) acc. to formula in Goode et al. (1994), page 147
      data(:,2) = data(:,2) * 1e-6 * 2 * pi .* data(:,1);  

      % to get data at 0dB SPL (assumed that stapes velocity is linearly related to pressure
      data(:,2) = data(:,2) * 10^(-104/20);

      % to get stapes PEAK velocity, multiply amplitudes by sqrt(2)
      data(:,2) = data(:,2).*sqrt(2);   


      % extrapolated data points, directly read from figure 2b) of Lopez-Poveda
      % and Meddis (2001)
      extrp = [100 1.181E-09; ...
               200	2.363E-09; ...
               7000 8.705E-10; ...
               7500 8.000E-10; ...
               8000 7.577E-10; ...
               8500 7.168E-10; ...
               9000 6.781E-10; ...
               9500 6.240E-10; ...
               10000 6.000E-10];

      if flags.do_plot
        figure;
        loglog(data(:,1)/1000,data(:,2),'ok', 'MarkerFaceColor', 'k');
        hold on
        loglog(extrp(:,1)/1000,extrp(:,2),'ok');
        xlabel('Frequency (kHz)');
        ylabel('Stapes velocity (m/s) at 0dB SPL)');
        axis([0.1,10,1e-10,1e-7]);
      end;

      data = [extrp(extrp(:,1) < data(1,1),:); data; extrp(extrp(:,1) > data(end,1),:)];

  elseif flags.do_lopezpoveda
    
    % data read from Poveda Fig.2, excl. extrapolated points
    data = [400	4.728E-09; ...
    600	7.577E-09; ...
    800	1.000E-08; ...
    1000 8.235E-09; ...
    1200 6.240E-09; ...
    1400 5.585E-09; ...
    1600 5.000E-09; ...
    1800 4.232E-09; ...
    2000 3.787E-09; ...
    2200 3.000E-09; ...
    2400 2.715E-09; ...
    2600 2.498E-09; ...
    2800 2.174E-09; ...
    3000 1.893E-09; ...
    3500 1.742E-09; ...
    4000 1.516E-09; ...
    4500 1.117E-09; ...
    5000 1.320E-09; ...
    5500 1.214E-09; ...
    6000 9.726E-10; ...
    6500 9.460E-10];

    % data read from Poveda Fig.2, extrapolated points
    extrp = [100 1.181E-09; ...
    200	2.363E-09; ...
    7000 8.705E-10; ...
    7500 8.000E-10; ...
    8000 7.577E-10; ...
    8500 7.168E-10; ...
    9000 6.781E-10; ...
    9500 6.240E-10; ...
    10000 6.000E-10];

    data = [extrp(extrp(:,1) < data(1,1),:); data; extrp(extrp(:,1) > data(end,1),:)];

    if flags.do_plot
      figure;
      loglog(data(:,1)/1000,data(:,2),'ok', 'MarkerFaceColor', 'k');
      hold on
      loglog(extrp(:,1)/1000,extrp(:,2),'ok');
      xlabel('Frequency (kHz)');
      ylabel('Stapes velocity (m/s) at 0dB SPL)');
      axis([0.1,10,1e-10,1e-7]);
    end;

  end
  
end;

