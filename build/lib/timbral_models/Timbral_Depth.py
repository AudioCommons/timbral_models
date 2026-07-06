from __future__ import division
import numpy as np
import soundfile as sf
from scipy.signal import spectrogram
import scipy.stats
from . import timbral_util


def timbral_depth(fname, fs=0, dev_output=False, phase_correction=False, clip_output=False, threshold_db=-60,
                  low_frequency_limit=20, centroid_crossover_frequency=2000, ratio_crossover_frequency=500,
                  db_decay_threshold=-40):
    """
     This function calculates the apparent Depth of an audio file.
     This version of timbral_depth contains self loudness normalising methods and can accept arrays as an input
     instead of a string filename.

     Version 0.4

     Required parameter
      :param fname:                        string or numpy array
                                           string, audio filename to be analysed, including full file path and extension.
                                           numpy array, array of audio samples, requires fs to be set to the sample rate.

     Optional parameters
      :param fs:                           int/float, when fname is a numpy array, this is a required to be the sample rate.
                                           Defaults to 0.
      :param phase_correction:             bool, perform phase checking before summing to mono.  Defaults to False.
      :param dev_output:                   bool, when False return the depth, when True return all extracted
                                           features.  Default to False.
      :param threshold_db:                 float/int (negative), threshold, in dB, for calculating centroids.
                                           Should be negative.  Defaults to -60.
      :param low_frequency_limit:          float/int, low frequency limit at which to highpass filter the audio, in Hz.
                                           Defaults to 20.
      :param centroid_crossover_frequency: float/int, crossover frequency for calculating the spectral centroid, in Hz.
                                           Defaults to 2000
      :param ratio_crossover_frequency:    float/int, crossover frequency for calculating the ratio, in Hz.
                                           Defaults to 500.

      :param db_decay_threshold:           float/int (negative), threshold, in dB, for estimating duration.  Should be
                                           negative.  Defaults to -40.

      :return:                             float, aparent depth of audio file, float.

     Copyright 2018 Andy Pearce, Institute of Sound Recording, University of Surrey, UK.

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
    '''
      Read input
    '''
    audio_samples, fs = timbral_util.file_read(fname, fs, phase_correction=phase_correction)

    '''
      Filter audio
    '''
    # highpass audio - run 3 times to get -18dB per octave - unstable filters produced when using a 6th order
    audio_samples = timbral_util.filter_audio_highpass(audio_samples, crossover=low_frequency_limit, fs=fs)
    audio_samples = timbral_util.filter_audio_highpass(audio_samples, crossover=low_frequency_limit, fs=fs)
    audio_samples = timbral_util.filter_audio_highpass(audio_samples, crossover=low_frequency_limit, fs=fs)

    # running 3 times to get -18dB per octave rolloff, greater than second order filters are unstable in python
    lowpass_centroid_audio_samples = timbral_util.filter_audio_lowpass(audio_samples,crossover=centroid_crossover_frequency,fs=fs)
    lowpass_centroid_audio_samples = timbral_util.filter_audio_lowpass(lowpass_centroid_audio_samples,crossover=centroid_crossover_frequency,fs=fs)
    lowpass_centroid_audio_samples = timbral_util.filter_audio_lowpass(lowpass_centroid_audio_samples,crossover=centroid_crossover_frequency,fs=fs)

    lowpass_ratio_audio_samples = timbral_util.filter_audio_lowpass(audio_samples,crossover=ratio_crossover_frequency,fs=fs)
    lowpass_ratio_audio_samples = timbral_util.filter_audio_lowpass(lowpass_ratio_audio_samples,crossover=ratio_crossover_frequency,fs=fs)
    lowpass_ratio_audio_samples = timbral_util.filter_audio_lowpass(lowpass_ratio_audio_samples,crossover=ratio_crossover_frequency,fs=fs)

    '''
      Get spectrograms and normalise
    '''
    # normalise audio
    lowpass_ratio_audio_samples *= (1.0 / max(abs(audio_samples)))
    lowpass_centroid_audio_samples *= (1.0 / max(abs(audio_samples)))
    audio_samples *= (1.0 / max(abs(audio_samples)))

    # set FFT parameters
    nfft = 4096
    hop_size = int(3 * nfft / 4)
    # get spectrogram
    if len(audio_samples) > nfft:
        freq, time, spec = spectrogram(audio_samples, fs, 'hamming', nfft, hop_size,
                                       nfft, 'constant', True, 'spectrum')
        lp_centroid_freq, lp_centroid_time, lp_centroid_spec = spectrogram(lowpass_centroid_audio_samples, fs,
                                                                           'hamming', nfft, hop_size, nfft,
                                                                           'constant', True, 'spectrum')
        lp_ratio_freq, lp_ratio_time, lp_ratio_spec = spectrogram(lowpass_ratio_audio_samples, fs, 'hamming', nfft,
                                                                  hop_size, nfft, 'constant', True, 'spectrum')

    else:
        # file is shorter than 4096, just take the fft
        freq, time, spec = spectrogram(audio_samples, fs, 'hamming', len(audio_samples), len(audio_samples)-1,
                                       nfft, 'constant', True, 'spectrum')
        lp_centroid_freq, lp_centroid_time, lp_centroid_spec = spectrogram(lowpass_centroid_audio_samples, fs,
                                                                           'hamming',
                                                                           len(lowpass_centroid_audio_samples),
                                                                           len(lowpass_centroid_audio_samples)-1,
                                                                           nfft, 'constant', True, 'spectrum')
        lp_ratio_freq, lp_ratio_time, lp_ratio_spec = spectrogram(lowpass_ratio_audio_samples, fs, 'hamming',
                                                                  len(lowpass_ratio_audio_samples),
                                                                  len(lowpass_ratio_audio_samples)-1,
                                                                  nfft, 'constant', True, 'spectrum')



    threshold = timbral_util.db2mag(threshold_db)


    '''
      METRIC 1 - limited weighted mean normalised lower centroid
    '''
    # define arrays for storing metrics
    all_normalised_lower_centroid = []
    all_normalised_centroid_tpower = []

    # get metrics for each time segment of the spectrogram
    for idx in range(len(time)):
        # get overall spectrum of time frame
        current_spectrum = spec[:, idx]
        # calculate time window power
        tpower = np.sum(current_spectrum)
        all_normalised_centroid_tpower.append(tpower)

        # estimate if time segment contains audio energy or just noise
        if tpower > threshold:
            # get the spectrum
            lower_spectrum = lp_centroid_spec[:, idx]
            lower_power = np.sum(lower_spectrum)

            # get lower centroid
            lower_centroid = np.sum(lower_spectrum * lp_centroid_freq) / float(lower_power)

            # append to list
            all_normalised_lower_centroid.append(lower_centroid)
        else:
            all_normalised_lower_centroid.append(0)

    # calculate the weighted mean of lower centroids
    weighted_mean_normalised_lower_centroid = np.average(all_normalised_lower_centroid,
                                                         weights=all_normalised_centroid_tpower)
    # limit to the centroid crossover frequency
    if weighted_mean_normalised_lower_centroid > centroid_crossover_frequency:
        limited_weighted_mean_normalised_lower_centroid = np.float64(centroid_crossover_frequency)
    else:
        limited_weighted_mean_normalised_lower_centroid = weighted_mean_normalised_lower_centroid



    '''
     METRIC 2 - weighted mean normalised lower ratio
    '''
    # define arrays for storing metrics
    all_normalised_lower_ratio = []
    all_normalised_ratio_tpower = []

    # get metrics for each time segment of the spectrogram
    for idx in range(len(time)):
        # get time frame of broadband spectrum
        current_spectrum = spec[:, idx]
        tpower = np.sum(current_spectrum)
        all_normalised_ratio_tpower.append(tpower)

        # estimate if time segment contains audio energy or just noise
        if tpower > threshold:
            # get the lowpass spectrum
            lower_spectrum = lp_ratio_spec[:, idx]
            # get the power of this
            lower_power = np.sum(lower_spectrum)
            # get the ratio of LF to all energy
            lower_ratio = lower_power / float(tpower)
            # append to array
            all_normalised_lower_ratio.append(lower_ratio)
        else:
            all_normalised_lower_ratio.append(0)

    # calculate
    weighted_mean_normalised_lower_ratio = np.average(all_normalised_lower_ratio, weights=all_normalised_ratio_tpower)

    '''
      METRIC 3 - Approximate duration/decay-time of sample 
    '''
    all_my_duration = []

    # get envelpe of signal
    envelope = timbral_util.sample_and_hold_envelope_calculation(audio_samples, fs)
    # estimate onsets
    onsets = timbral_util.calculate_onsets(audio_samples, envelope, fs)

    # get RMS envelope - better follows decays than the sample-and-hold
    rms_step_size = 256
    rms_envelope = timbral_util.calculate_rms_enveope(audio_samples, step_size=rms_step_size)

    # convert decay threshold to magnitude
    decay_threshold = timbral_util.db2mag(db_decay_threshold)
    # rescale onsets to rms stepsize - casting to int
    time_convert = fs / float(rms_step_size)
    onsets = (np.array(onsets) / float(rms_step_size)).astype('int')

    for idx, onset in enumerate(onsets):
        if onset == onsets[-1]:
            segment = rms_envelope[onset:]
        else:
            segment = rms_envelope[onset:onsets[idx + 1]]

        # get location of max RMS frame
        max_idx = np.argmax(segment)
        # get the segment from this max until the next onset
        post_max_segment = segment[max_idx:]

        # estimate duration based on decay or until next onset
        if min(post_max_segment) >= decay_threshold:
            my_duration = len(post_max_segment) / time_convert
        else:
            my_duration = np.where(post_max_segment < decay_threshold)[0][0] / time_convert

        # append to array
        all_my_duration.append(my_duration)

    # calculate the lof of mean duration
    mean_my_duration = np.log10(np.mean(all_my_duration))


    '''
      METRIC 4 - f0 estimation with peak picking
    '''
    # get the overall spectrum
    all_spectrum = np.sum(spec, axis=1)
    # normalise this
    norm_spec = (all_spectrum - np.min(all_spectrum)) / (np.max(all_spectrum) - np.min(all_spectrum))
    # set limit for peak picking
    cthr = 0.01
    # detect peaks
    peak_idx, peak_value, peak_freq = timbral_util.detect_peaks(norm_spec, cthr=cthr, unprocessed_array=norm_spec,
                                                                freq=freq)
    # estimate peak
    pitch_estimate = np.log10(min(peak_freq)) if peak_freq[0] > 0 else 0


    # get outputs
    if dev_output:
        return limited_weighted_mean_normalised_lower_centroid, weighted_mean_normalised_lower_ratio, mean_my_duration, \
               pitch_estimate, weighted_mean_normalised_lower_ratio * mean_my_duration, \
               timbral_util.sigmoid(weighted_mean_normalised_lower_ratio) * mean_my_duration
    else:
        '''
         Perform linear regression to obtain depth
        '''
        # coefficients from linear regression
        coefficients = np.array([-0.0043703565847874465, 32.83743202462131, 4.750862716905235, -14.217438690256062,
                                 3.8782339862813924, -0.8544826091735516, 66.69534393444391])

        # what are the best metrics
        metric1 = limited_weighted_mean_normalised_lower_centroid
        metric2 = weighted_mean_normalised_lower_ratio
        metric3 = mean_my_duration
        metric4 = pitch_estimate
        metric5 = metric2 * metric3
        metric6 = timbral_util.sigmoid(metric2) * metric3

        # pack metrics into a matrix
        all_metrics = np.zeros(7)

        all_metrics[0] = metric1
        all_metrics[1] = metric2
        all_metrics[2] = metric3
        all_metrics[3] = metric4
        all_metrics[4] = metric5
        all_metrics[5] = metric6
        all_metrics[6] = 1.0

        # perform linear regression
        depth = np.sum(all_metrics * coefficients)

        if clip_output:
            depth = timbral_util.output_clip(depth)

        return depth

