from __future__ import division
import numpy as np
import soundfile as sf
from scipy.signal import spectrogram
from . import timbral_util
import scipy.stats
from sklearn import linear_model


def timbral_warmth(fname, dev_output=False, phase_correction=False, clip_output=False, max_FFT_frame_size=8192,
                   rolloff_ratio=0.85, max_WR = 12000):
    """
     This function estimates the perceptual Warmth of an audio file.

     Version 0.2

     This has been calibrated on a small data set and is yet to be subjectively validated.

     Required parameter
    :param fname:                   Audio filename to be analysed, including full file path and extension.

    Optional parameters
    :param dev_output:              Bool, when False return the warmth, when True return all extracted features.
    :param phase_correction:        If the inter-channel phase should be estimated when performing a mono sum.
                                    Defaults to False.
    :param max_FFT_frame_size:      Frame size for calculating spectrogram, default to 8192.
    :param rolloff_ratio:           Ratio for calculating the spectral rolloff, detault so 0.85 meaning the frequency
                                    where 85% of the energy sits below.
    :param max_WR:                  Maximun allowable warmth region frequency, defaults to 12000.

    :return:                        Estimated warmth of audio file.

    Copyright 2018 Andy Pearce

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

    """
    # use pysoundfile to read audio
    audio_samples, fs = sf.read(fname, always_2d=False)

    # channel reduction
    audio_samples = timbral_util.channel_reduction(audio_samples, phase_correction=phase_correction)

    # normalise audio
    audio_samples /= max(abs(audio_samples))

    '''
      Estimate a fundamental frequency
    '''
    # get FFT of signal
    audio_length = len(audio_samples)
    if audio_length < max_FFT_frame_size:
        freq, time, spec = spectrogram(audio_samples, fs, nperseg=audio_length, nfft=max_FFT_frame_size)
    else:
        freq, time, spec = spectrogram(audio_samples, fs, nperseg=max_FFT_frame_size, nfft=max_FFT_frame_size)
        if spec.shape[1] > 1:
            spec = np.sum(spec, axis=1)
            spec = spec.flatten()

    # flatten and normalise spectrogram
    spec = np.array(list(spec)).flatten()
    this_shape = spec.shape
    spec /= max(abs(spec))

    # peak picking algorithm
    peak_idx, peak_value, peak_x = timbral_util.detect_peaks(spec, freq=freq, fs=fs)
    # find lowest peak
    fundamental = np.min(peak_x)
    fundamental_idx = np.min(peak_idx)

    '''
      Warmth region calculation
    '''
    # estimate the Warmth region
    WR_upper_f_limit = fundamental * 3.5
    if WR_upper_f_limit > max_WR:
        WR_upper_f_limit = 12000
    tpower = np.sum(spec)

    WR_upper_f_limit_idx = int(np.where(freq>WR_upper_f_limit)[0][0])
    WR_energy = np.sum(spec[:WR_upper_f_limit_idx])
    WR_ratio = WR_energy / float(tpower)

    '''
      Spectral rolloff
    '''
    tpower_above_fundamental = np.sum(spec[fundamental_idx:])
    thresh = rolloff_ratio * tpower_above_fundamental
    cumsum = np.cumsum(spec[fundamental_idx:]) / tpower_above_fundamental
    rolloff_idx = np.where(cumsum > rolloff_ratio)[0][0]
    rolloff = freq[rolloff_idx+fundamental_idx]

    # estimate of fundamental-to-260Hz :  overall energy
    # need to check if the fundamental is greater than the
    top_level_idx = int(np.where(freq>260)[0][0])
    # sum energy up to this bin
    low_energy = np.sum(spec[:top_level_idx])
    # sum all energy
    tpower = np.sum(spec)
    # take ratio
    ratio = low_energy / float(tpower)

    '''
      Spectral centroid
    '''
    # spectral centroid
    SC = np.sum(freq * spec) / float(np.sum(spec))

    '''
      Difference metrics
    '''
    difference = rolloff - WR_upper_f_limit
    WR_to_rolloff = np.log10(abs(difference)) * np.sign(difference)  # cannot evaluate negative numbers in log

    '''
      Continuous/onsetoffset richness
    '''
    hp_ratio = timbral_util.get_percussive_audio(audio_samples, return_ratio=True)

    '''
      HF decay
      - linear regression of the values above WR
    '''
    above_WR_spec = np.log10(spec[WR_upper_f_limit_idx:])
    above_WR_freq = np.log10(freq[WR_upper_f_limit_idx:])
    np.ones_like(above_WR_freq)
    metrics = np.array([above_WR_freq, np.ones_like(above_WR_freq)])

    model = linear_model.LinearRegression(fit_intercept=False)
    model.fit(metrics.transpose(), above_WR_spec)
    decay_score = model.score(metrics.transpose(), above_WR_spec)


    if dev_output:
        return WR_to_rolloff, np.log10(rolloff), np.log10(SC), WR_ratio, hp_ratio, decay_score
    else:
        coefficients = np.array([5.9564072511755555, -10.217111190953734, -18.6868387124319, 13.490838413321452,
                                 -18.3376612134757, 16.122840820804377, 105.6071349465704])

        metric1 = WR_to_rolloff
        metric2 = np.log10(rolloff)
        metric3 = np.log10(SC)
        metric4 = WR_ratio
        metric5 = hp_ratio
        metric6 = decay_score
        metric_ones = np.ones_like(metric1)
        metrics = np.array([metric1, metric2, metric3, metric4, metric5, metric6, metric_ones])

        warmth = np.sum(coefficients * metrics)

        if clip_output:
            warmth = timbral_util.output_clip(warmth)

        return warmth
