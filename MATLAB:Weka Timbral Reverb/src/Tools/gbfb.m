function features = gbfb(log_mel_spec, omega_max, size_max, nu, distance)
% usage: features = gbfb(log_mel_spec)
%  log_mel_spec    log Mel-spectrogram or similar spectro-temporal representation
%
% usage: features = gbfb(log_mel_spec, omega_max, size_max, nu, distance)
%  omega_max       upper center modulation frequencies in [radian radian]
%  size_max        maximum filter size (lower center modulation frequency) in [bands, frames]
%  nu              half-waves under the envelope in [spectral temporal] dimension
%  distance        [spectral temporal] spacing of filters (<1)
%
% - Gabor filter bank features v3.0 -
%
% Autor    : Marc René Schädler
% Email    : marc.r.schaedler@uni-oldenburg.de
% Institute: Medical Physics / Carl-von-Ossietzky University Oldenburg, Germany
%
%-----------------------------------------------------------------------------
%
% Release Notes:
% v2.0 - Inital release
% v3.0 - Fixed hann_win bug
% v3.0 - Removed log Mel-spectrogram calculation
% v3.0 - Added caching and simplified feature vector composition
%

 
%% Default settings and checks

% Default upper center modulation frequencies
if nargin < 2 || isempty(omega_max)
  omega_max = [pi/2 pi/2];
end

% Get number of bands from input
num_bands = size(log_mel_spec,1);

% Default number of half-waves under the envelope
if nargin < 4 || isempty(nu)
  nu = [3.5 3.5];
end

% Default maximum filter size (lower center modulation frequencies)
if nargin < 3 || isempty(size_max)
  size_max = [3*num_bands 40];
end

% Default spacing of filters
if nargin < 5 || isempty(distance)
  distance = [0.3 0.2];
end

% Maximum context
context = floor(size_max(2)/2);


%% Calculate Gabor filter bank features

% Temporally pad log Mel-spectrogram by repeating first and last frames
log_mel_spec = [repmat(log_mel_spec(:,1),1,context) log_mel_spec repmat(log_mel_spec(:,end),1,context)];

% Calculate center modulation frequencies.
[omega_n, omega_k] = gfbank_calcaxis(omega_max, size_max, nu, distance);
omega_n_num = length(omega_n);
omega_k_num = length(omega_k);

% Generate filters for all pairs of spectral and temporal modulation
% frequencies except for the redundant ones.
gfilters = cell(omega_k_num,omega_n_num);
for i=1:omega_k_num
  for j=1:omega_n_num
    if ~(omega_k(i) < 0 && omega_n(j) == 0)
      gfilters{i,j} = gfilter_gen(omega_k(i), omega_n(j), nu(1), nu(2), size_max(1), size_max(2));
    end
  end
end

% Filter mel spectrogram with filter bank filters and select representative channels.
gfilters_output = cell(numel(gfilters),1);
for i=1:length(gfilters_output)
  gfilter = gfilters{i};
  if ~isempty(gfilter)
    % Filter mel spectrogram with Gabor filter.
    log_mel_spec_filtered = gfilter_filter(gfilter, log_mel_spec);
    % Select representative channels from filtered Mel-spectrogram.
    gfilters_output{i} = gfilter_rep(gfilter, log_mel_spec_filtered);
  end
end
features = cell2mat(gfilters_output);

% Use the real part of the filter output
features = real(features);

% Remove padded context
features = features(:,(1+context):(end-context));
end


function [omega_n, omega_k] = gfbank_calcaxis(omega_max, size_max, nu, distance)
% Calculates the modulation center frequencies iteratively.
% Termination condition for iteration is reaching omega_min, which is
% derived from size_max.
omega_min = (pi .* nu) ./ size_max;

% Eq. (2b)
c = distance .* 8 ./ nu;
% Second factor of Eq. (2a)
space_n = (1 + c(2)./2) ./ (1 - c(2)./2);
count_n = 0;
omega_n(1) = omega_max(2);
% Iterate starting with omega_max in spectral dimension
while omega_n/space_n > omega_min(2)
  omega_n(1+count_n) = omega_max(2)/space_n.^count_n;
  count_n = count_n + 1;
end
omega_n = fliplr(omega_n);
% Add DC
omega_n = [0,omega_n];
% Second factor of Eq. (2a)
space_k = (1 + c(1)./2) ./ (1 - c(1)./2);
count_k = 0;
omega_k(1) = omega_max(1);
% Iterate starting with omega_max in temporal dimension
while omega_k/space_k > omega_min(1)
  omega_k(1+count_k) = omega_max(1)/space_k.^count_k;
  count_k = count_k + 1;
end
% Add DC and negative MFs for spectro-temporal opposite 
% filters (upward/downward)
omega_k = [-omega_k,0,fliplr(omega_k)];
end


function gfilter = gfilter_gen(omega_k, omega_n, nu_k, nu_n, size_max_k, size_max_n)
% Generates a gabor filter function with:
%  omega_k       spectral mod. freq. in rad
%  omega_n       temporal mod. freq. in rad
%  nu_k          number of half waves unter the envelope in spectral dim.
%  nu_n          number of half waves unter the envelope in temporal dim.
%  size_max_k    max. allowed extension in spectral dimension
%  size_max_n    max. allowed extension in temporal dimension

% Build a config id string
config = strrep(sprintf('c%.0f', [omega_k omega_n nu_k nu_n size_max_k size_max_n]*10000),'-','_');

% Load cache
persistent cache;

% Only generate filters which are not cached
if isempty(cache) || ~isfield(cache, config)
  % Calculate windows width
  w_n = 2*pi / abs(omega_n) * nu_n / 2;
  w_k = 2*pi / abs(omega_k) * nu_k / 2;

  % If the size exceeds the max. allowed extension in a dimension set the
  % corresponding mod. freq. to zero.
  if w_n > size_max_n
    w_n = size_max_n;
    omega_n = 0;
  end
  if w_k > size_max_k
    w_k = size_max_k;
    omega_k = 0;
  end

  % Separable hanning envelope, cf. Eq. (1c).
  env_n = hann_win(w_n);
  env_k = hann_win(w_k);
  envelope = env_k * env_n.';
  [win_size_k, win_size_n] = size(envelope);

  % Sinusoid carrier, cf. Eq. (1c).
  n_0 = (win_size_n+1) / 2;
  k_0 = (win_size_k+1) / 2;
  [n,k] = meshgrid (1:win_size_n, 1:win_size_k);
  sinusoid = exp(1i*omega_n*(n - n_0) + 1i*omega_k*(k - k_0));

  % Eq. 1c
  gfilter  = envelope .* sinusoid;

  % Compensate the DC part by subtracting an appropiate part
  % of the envelope if filter is not the DC filter.
  envelope_mean = mean(mean(envelope));
  gfilter_mean = mean(mean(gfilter));
  if (omega_n ~= 0) || (omega_k ~= 0)
    gfilter = gfilter - envelope./envelope_mean .* gfilter_mean;
  else
    % Add an imaginary part to DC filter for a fair real/imag comparison.
    gfilter = gfilter + 1i*gfilter;
  end
  
  % Normalize filter to have gains <= 1.
  gfilter = gfilter ./ max(max(abs(fft2(gfilter))));
  
  % Save to cache
  cache.(config) = gfilter;
else
  % Load from cache
  gfilter = cache.(config);
end
end


function log_mel_spec_filt = gfilter_filter(gfilter, log_mel_spec)
% Applies the filtering with a 2D Gabor filter to log_mel_spec
% This includes the special treatment of filters that do not lie fully
% inside the spectrogram
if any(any(gfilter < 0))
  % Compare this code to the compensation for the DC part in the
  % 'gfilter_gen' function. This is an online version of it removing the
  % DC part of the filters by subtracting an appropriate part of the
  % filters' envelope.
  gfilter_abs_norm = abs(gfilter) ./ sum(sum(abs(gfilter)));
  gfilter_dc_map = fftconv2(ones(size(log_mel_spec)), gfilter,'same');
  env_dc_map = fftconv2(ones(size(log_mel_spec)), gfilter_abs_norm,'same');
  dc_map = fftconv2(log_mel_spec, gfilter_abs_norm,'same') ./ env_dc_map .* gfilter_dc_map;
else
  % Dont' remove the DC part if it is the DC filter.
  dc_map = 0;
end
% Filter log_mel_spec with the 2d Gabor filter and remove the DC parts.
log_mel_spec_filt = fftconv2(log_mel_spec, gfilter,'same') - dc_map;
end


function mel_spec_rep = gfilter_rep(gfilter, mel_spec)
% Selects the center channel by choosing k_offset and those with k_factor
% channels distance to it in spectral dimension where k_factor is approx.
% 1/4 of the filters extension in the spectral dimension.
k_factor = floor(1/4 * size(gfilter,1));
if k_factor < 1
    k_factor = 1;
end
k_offset = mod(floor(size(mel_spec,1)/2),k_factor);
k_idx = (1+k_offset):k_factor:size(mel_spec,1);
% Apply subsampling.
mel_spec_rep = mel_spec(k_idx,:);
end


function window_function = hann_win(width)
% A hanning window of "width" with the maximum centered on the center sample
x_center = 0.5;
step = 1/width;
right = x_center:step:1;
left = x_center:-step:0;
x_values = [left(end:-1:1) right(2:end)].';
valid_values_mask = (x_values > 0) & (x_values < 1);
window_function = 0.5 * (1 - ( cos(2*pi*x_values(valid_values_mask))));
end


function out = fftconv2(in1, in2, shape)
% 2D convolution in terms of the 2D FFT that substitutes conv2(in1, in2, shape).
size_y = size(in1,1)+size(in2,1)-1;
size_x = size(in1,2)+size(in2,2)-1;
fft_size_x = 2.^ceil(log2(size_x));
fft_size_y = 2.^ceil(log2(size_y));
in1_fft = fft2(in1,fft_size_y,fft_size_x);
in2_fft = fft2(in2,fft_size_y,fft_size_x);
out_fft = in1_fft .* in2_fft;
out_padd = ifft2(out_fft);
out_padd = out_padd(1:size_y,1:size_x);
switch shape
  case 'same'
    y_offset = floor(size(in2,1)/2);
    x_offset = floor(size(in2,2)/2);
    out = out_padd(1+y_offset:size(in1,1)+y_offset,1+x_offset:size(in1,2)+x_offset);
  case 'full'
    out = out_padd;
end
end

