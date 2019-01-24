function par = init_rt_estimate_e( fs )
% par = init_rt_estimate_e( fs )
% executes initialization for the function 
% rt_estimate_frame.m to perform a blind estimation of the reverberation time
% (RT) bu frame-wise processing in the time-domain.
%
% INPUT 
% fs : sampling frequency (default = 24kHz)
%
% OUTPUT
% par: struct containing all parameters and buffer for executing the
%      function rt_estimate_frame.m 
%
% author: Heiner Loellmann, IND, RWTH Aachen University
%
% created: August 2011

% general paraemters
if nargin == 0
    par.fs = 24e3;
else
   par.fs = fs;
end

no = par.fs / 24e3 ;  % correction factor to account for different sampling frequency

% pararmeters for pre-selection of suitable segments
if par.fs>8e3
    par.down = 2;                           % rate for downsampling applied before RT estimation to reduce computational complexity
else
    par.down = 1;
end
par.N_sub = round( no*700/par.down);  % sub-frame length (after downsampling)
par.N_shift = round( no*200/par.down);  % frame shift (before downsampling)
par.nos_min = 3;                        % minimal number of subframes to detect a sound decay
par.nos_max = 7;                        % maximal number of subframes to detect a sound decay
par.N = par.nos_max*par.N_sub;          % maximal frame length (after downsampling)

% parameters for ML-estimation
Tmax = 1.1;                  % max RT being considered
Tmin = 0.2;                  % min RT being considered
par.bin = 0.1;                             % step-size for RT estimation
par.Tquant = ( Tmin : par.bin : Tmax );    % set of qunatized RTs considered for maximum search
par.a = exp( - 3*log(10) ./ ( (par.Tquant) .* (par.fs/par.down)));   % corresponding decay rate factors
par.La = length( par.a );                  % no. of considered decay rate factors (= no of. RTs)

% paramters for histogram-based approach to reduce outliers (order statistics)
par.buffer_size = round( no*800/par.down);% buffer size
par.buffer = zeros( 1, par.buffer_size ); % buffer with previous indices to update histogram
par.no_bins  = par.La;                    % no. of histogram bins
par.hist_limits = Tmin - par.bin/2 : par.bin :  Tmax + par.bin/2 ; % limits of histogram bins
par.hist_rt = zeros(1,par.no_bins);       % histogram with ML estimates
par.hist_counter = 0;                     % counter increased if histogram is updated

% paramters for recursive smoothing of final RT estimate
par.alpha = 0.995;              % smoothing factor
par.RT_initial = 0.3;          % initial RT estimate
par.RT_last = par.RT_initial;  % last RT estimate
par.RT_raw  =  par.RT_initial; % raw RT estimate obtained by histogram-approach