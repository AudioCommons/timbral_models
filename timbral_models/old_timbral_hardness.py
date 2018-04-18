import numpy as np
import librosa
import soundfile as sf
from scipy.signal import butter, lfilter, spectrogram
import matplotlib.pyplot as plt
import essentia.standard as es
import pyfilterbank
from essentia import Pool, array


def db2mag(dB):
    mag = 10 ** (dB / 20.0)
    return mag


def get_percussive_audio(audio_samples, show_ratio=False):
    """
      Gets the percussive comonent of the audio file
    :param audio_samples:
    :return:
    """
    ''' Edit the percussive margin to affect extent of extraction '''
    D = librosa.core.stft(audio_samples)
    H, P = librosa.decompose.hpss(D)
    # H, P = librosa.decompose.hpss(D, margin=(1.0, 5.0))
    percussive_audio = librosa.core.istft(P)
    harmonic_audio = librosa.core.istft(H)

    if show_ratio:
        # frame by frame RMS energy
        percussive_energy = []
        harmoinic_energy = []
        for frame in es.FrameGenerator(percussive_audio, frameSize=1024, hopSize=512):
            percussive_energy.append(es.RMS()(frame))
        for frame in es.FrameGenerator(harmonic_audio, frameSize=1024, hopSize=512):
            harmoinic_energy.append(es.RMS()(frame))
        ratio = []
        t_power = []
        for i in range(len(percussive_energy)):
            if percussive_energy[i] != 0 and harmoinic_energy != 0:
                ratio.append(percussive_energy[i] / (percussive_energy[i] + harmoinic_energy[i]))
                t_power.append((percussive_energy[i] + harmoinic_energy[i]))

        ratio = np.average(ratio, weights=t_power)

        # hard limit max
        # if ratio > 0.5:
        #     ratio = 0.5

        return ratio
    else:

        # plt.subplot(2, 1, 1)
        # plt.plot(audio_samples)
        # plt.subplot(2, 1, 2)
        # plt.plot(percussive_audio)
        # plt.show()

        return percussive_audio


def filter_audio_highpass(audio_samples, crossover, fs, order=5):
    """ Calculate and apply a high-pass filter, with a -3dB point of crossover.

    :param audio_samples: the audio file as an array
    :param crossover: the crossover frequency of the filter
    :param fs: the sampling frequency of the audio file
    :param order: order of the filter, defaults to 5
    :return: audio file filtered

    """
    nyq = 0.5 * fs
    xfreq = crossover / nyq
    b, a = butter(order, xfreq, 'high')
    y = lfilter(b, a, audio_samples)
    return y


def butter_bandpass(lowcut, highcut, fs, order=2):
    """ Design a butterworth bandpass filter """
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def filter_audio_bandpass(audio_samples, f0, noct, fs, order=2):
    """ Calculate and apply an n/octave butterworth bandpass filter, centred at f0 Hz.

    :param audio_samples: the audio file as an array
    :param fs: the sampling frequency of the audio file
    :param f0: the centre frequency of the bandpass filter
    :param bandwidth: the bandwidth of the filter
    :param order: order of the filter, defaults to 2
    :return: audio file filtered

    """
    fd = 2 ** (1.0 / (noct * 2))
    lowcut = f0 / fd
    highcut = f0 * fd

    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, audio_samples)
    return y


def return_loop(onset_loc, envelope, function_time_thresh, hist_threshold, hist_time_samples, nperseg=512):
    """ This function is used by the timbral_hardness method.
     This looks backwards in time from the attack time and attempts to find the exact onset point by
     identifying the point backwards in time where the envelope no longer falls.
     This function includes a hyteresis to account for small deviations in the attack due to the
     envelope calculation.

     Function looks 10ms (function_time_thresh) backwards from the onset time (onset_loc), looking for any sample
     lower than the current sample.  This repeats, starting at the minimum value until no smaller value is found.
     Then the function looks backwards over 200ms, checking if the increase is greater than 10% of the full envelope's
     dynamic range.

        onset_loc:              The onset location estimated by librosa (converted to time domain index)
        envelope:               Envelope of the audio file
        function_time_thresh:   Time threshold for looking backwards in time.  Set in the timbral_hardness code
                                to be the number of samples that equates to 10ms
        hist_threshold:         Level threshold to check over 200ms if the peak is small enough to continue looking
                                backwards in time.
        hist_time_samples:      Number of samples to look back after finding the minimum value over 10ms, set to 200ms.
    """

    # define flag for exiting while loop
    found_start = False

    while not found_start:
        # get the current sample value
        current_sample = envelope[int(onset_loc)]
        # get the previous 10ms worth of samples
        if onset_loc - function_time_thresh >= 0:
            evaluation_array = envelope[onset_loc - function_time_thresh - 1:onset_loc]
        else:
            evaluation_array = envelope[:onset_loc - 1]

        if min(evaluation_array) - current_sample <= 0:
            '''
             If the minimum value within previous 10ms is less than current sample,
             move to the start position to the minimum value and look again.
            '''
            min_idx = np.argmin(evaluation_array)
            new_onset_loc = min_idx + onset_loc - function_time_thresh - 1

            if new_onset_loc > nperseg:
                onset_loc = new_onset_loc
            else:
                ''' Current index is close to start of the envelope, so exit with the idx as 512 '''
                return 0

        else:
            '''
             If the minimum value within previous 10ms is greater than current sample,
             introduce the time and level hysteresis to check again.
            '''
            # get the array of 200ms previous to the current onset idx
            if (onset_loc - hist_time_samples - 1) > 0:
                hyst_evaluation_array = envelope[onset_loc - hist_time_samples - 1:onset_loc]
            else:
                hyst_evaluation_array = envelope[:onset_loc]

            # values less than current sample
            all_match = np.where(hyst_evaluation_array < envelope[onset_loc])

            # if no minimum was found within the extended time, exit with current onset idx
            if len(all_match[0]) == 0:
                return onset_loc

            # get the idx of the closest value which is lower than the current onset idx
            last_min = all_match[0][-1]
            last_idx = int(onset_loc - len(hyst_evaluation_array) + last_min)

            # get the dynamic range of this segment
            segment_dynamic_range = max(hyst_evaluation_array[last_min:]) - min(hyst_evaluation_array[last_min:])

            # compare this dynamic range against the hyteresis threshold
            if segment_dynamic_range >= hist_threshold:
                '''
                 The dynamic range is greater than the threshold, therefore this is a separate audio event.
                 Return the current onset idx.
                '''
                return onset_loc
            else:
                '''
                 The dynamic range is less than the threshold, therefore this is not a separate audio event.
                 Set current onset idx to minimum value and repeat.
                '''
                if last_idx >= nperseg:
                    onset_loc = last_idx
                else:
                    '''
                     The hysteresis check puts the new threshold too close to the start
                    '''
                    return 0


def sample_and_hold_envelope_calculation(audio_samples, fs, decay_time=0.2, hold_time=0.01):
    """
     Calculates the envelope of audio_samples with a 'sample and hold' style function.
     This ensures that the minimum attack time is not limited by low-pass filtering,
     a common method of obtaining the envelope.

    :param audio_samples:
    :param fs: sampling frequency
    :param decay_time: decay time after peak hold
    :param hold_time:
    :return:
    """
    # rectify the audio signal
    abs_samples = abs(audio_samples)
    envelope = []

    # set parameters for envelope function
    decay = max(abs_samples) / (decay_time * fs)  # decay rate relative to peak level of audio signal
    hold_samples = hold_time * fs  # number of samples to hold before decay
    hold_counter = 0
    previous_sample = 0.0

    # perform the sample, hold, and decay function to obtain envelope
    for sample in abs_samples:
        if sample >= previous_sample:
            envelope.append(sample)
            previous_sample = sample
            hold_counter = 0
        else:
            # check hold length
            if hold_counter < hold_samples:
                hold_counter += 1
                envelope.append(previous_sample)
            else:
                out = previous_sample - decay
                if out > sample:
                    envelope.append(out)
                    previous_sample = out
                else:
                    envelope.append(sample)
                    previous_sample = sample

    # convert to numpy array
    return np.array(envelope)


def get_spectral_features(audio, fs, lf_limit=20, scale='hz', cref=27.5, power=2, window_type='none',
                          rollon_thresh=0.05):
    """
     This function calculates the spectral centroid and spectral spread of an audio array.

     :param audio:      Audio array
     :param fs:         Sample rate of audio file
     :param lf_limit:   Low frequency limit, in Hz, to be analysed.  Defaults to 20Hz.
     :param scale:      The frequency scale that calculations should be made over.  if no argument is given, this
                        defaults to 'hz', representing a linear frequency scale.  Options are 'hz', 'mel', 'erb',
                        or 'cents'.
     :param cref:       The reference frequency for calculating cents.  Defaults to 27.5Hz.
     :param power:      The power to raise devaition from specteal centroid, defaults to 2.
     :return:           Returns the spectral centroid and spectral spread.
    """
    # use a hanning window
    if window_type == 'hann':
        window = np.hanning(len(audio))
    elif window_type == 'none':
        window = np.ones(len(audio))
    else:
        raise ValueError('Window type must be set to either \'hann\' or \'none\'')

    next_pow_2 = int(pow(2, np.ceil(np.log2(len(window)))))
    # get frequency domain representation
    spectrum = np.fft.fft((window * audio), next_pow_2)
    spectrum = np.absolute(spectrum[0:int(len(spectrum) / 2) + 1])

    tpower = np.sum(spectrum)

    if tpower > 0:
        freq = np.arange(0, len(spectrum), 1) * (fs / (2.0 * (len(spectrum) - 1)))

        # find lowest frequency index, zeros used to unpack result
        lf_limit_idx = np.where(freq >= lf_limit)[0][0]
        spectrum = spectrum[lf_limit_idx:]
        freq = freq[lf_limit_idx:]

        # convert frequency to desired frequency scale
        if scale == 'hz':
            freq = freq
        elif scale == 'mel':
            freq = 1127.0 * np.log(1 + (freq / 700.0))
        elif scale == 'erb':
            freq = 21.4 * np.log10(1 + (0.00437 * freq))
        elif freq == 'cents':
            # for cents, a
            freq = 1200.0 * np.log2((freq / cref) + 1.0)
        else:
            raise ValueError('Frequency scale type not recognised.  Please use \'hz\', \'mel\', \'erb\', or \'cents\'.')

        # calculate centroid and spread
        centroid = sum(spectrum * freq) / float(sum(spectrum))

        # old calculation of spread
        deviation = np.abs(freq - centroid)
        spread = np.sqrt(np.sum((deviation ** 2) * spectrum) / np.sum(spectrum))

        # new calculation of spread according to librosa
        # spread = np.sqrt(np.sum(spectrum * (deviation ** power))) #** (1. / power))

        cumulative_spectral_power = spectrum[0]
        counter = 0
        rollon_threshold = np.sum(spectrum) * rollon_thresh
        while cumulative_spectral_power < rollon_threshold:
            counter += 1
            cumulative_spectral_power = np.sum(spectrum[:counter])

        if counter == 0:
            counter = 1

        rollon_frequency = freq[counter]
        unitless_centroid = centroid / rollon_frequency

        return centroid, spread, unitless_centroid
    else:
        return 0


def calculate_attack_time(envelope_samples, fs, calculate_attack_segment=True, thresh_no=8, normalise=True, m=3,
                          calculation_type='min_effort', gradient_calulation_type='all', return_descriptive_data=False,
                          max_attack_time=-1):
    """

    :param envelope_samples:
    :param fs:
    :param calculate_attack_segment:
    :param thresh_no:
    :param normalise:
    max_attack_time:    sets the maximum allowable attack time.  Defaults to -1, indicating that there is no maximum
                        attack time.  This value should be set in seconds.
    :param m:                           value used for computation of effort thresholds, set as 3 in CUIDADO project
    :return:
    """
    if normalise:
        # normalise the segments
        normalise_factor = float(max(envelope_samples))
        envelope_samples /= normalise_factor

    if calculate_attack_segment:
        # identify pre-attack segment
        peak_idx = np.argmax(envelope_samples)
        if peak_idx == 0:
            # exit on error
            return 0
        # min_pre_peak_idx = np.argmin(envelope_samples[:peak_idx])
        min_pre_peak_idx = np.where(envelope_samples[:peak_idx] == min(envelope_samples[:peak_idx]))[-1][-1]

        # redefine the envelope samples as just the min to the peak
        envelope_samples = envelope_samples[min_pre_peak_idx:peak_idx + 1]
    else:
        min_pre_peak_idx = 0

    # calculate the appropriate start and end of the attack using the selected method
    if calculation_type == 'min_effort':
        # get threshold time array
        threshold_step = 1.0 / (thresh_no + 2)  # +2 is to ignore the 0 and 100% levels.
        dyn_range = max(envelope_samples) - min(envelope_samples)
        thresh_level = np.linspace(threshold_step, (1 - threshold_step), thresh_no + 1)
        thresh_level = (thresh_level * dyn_range) + min(envelope_samples)

        # predefine an array for when each threshold is crossed
        threshold_idxs = np.zeros(thresh_no + 1)

        # get indexes for when threshold is crossed
        for j in range(len(thresh_level)):
            threshold_hold = np.argmax(envelope_samples >= thresh_level[j])
            # threshold_idxs[j] = threshold_hold + min_pre_peak_idx
            threshold_idxs[j] = threshold_hold

        # calculate effort values (distances between thresholds)
        effort = np.diff(threshold_idxs)

        # get the mean effort value
        effort_mean = np.mean(effort)
        effort_threshold = effort_mean * m

        # find start and stop times foxr the attack
        th_start = np.argmax(effort <= effort_threshold)

        # need to use remaining effort values
        effort_hold = effort[th_start:]
        th_end = np.argmax(effort_hold >= effort_threshold)  # this returns a 0 if value not found
        if th_end == 0:
            th_end = len(effort_hold) - 1  # make equal to the last value

        # apply correction for holding the values
        th_end = th_end + th_start

        # get the actual start and stop index
        th_start_idx = threshold_idxs[th_start]
        th_end_idx = threshold_idxs[th_end]

        if th_start_idx == th_end_idx:
            th_start_idx = threshold_idxs[0]
            th_end_idx = threshold_idxs[-1]

        if th_start_idx == th_end_idx:
            attack_time = 1.0 / fs
        else:
            attack_time = (th_end_idx - th_start_idx + 1.0) / fs

        if max_attack_time > 0:
            if attack_time > max_attack_time:
                # how many samples is equivalent to the maximum?
                max_attack_time_sample = int(fs * max_attack_time)  # convert to integer
                th_end_idx = th_start_idx + max_attack_time_sample
                attack_time = (th_end_idx - th_start_idx + 1.0) / fs

        start_level = envelope_samples[int(th_start_idx)]
        end_level = envelope_samples[int(th_end_idx)]

        # specify exceptions for a step functions crossing both thresholds
        if start_level == end_level:
            if th_start_idx > 0:
                # if a previous sample is avaiable, take the previous starting sample
                start_level = envelope_samples[int(th_start_idx) - 1]
            else:
                # set start level to zero if onset is at the first sample (indicating a step function at time zero)
                start_level = 0.0

        # is there enough data to calculate the mean
        if gradient_calulation_type == 'mean':
            if (end_level - start_level) < 0.2 or (th_end_idx - th_start_idx) < 2:
                # force calculation type to all
                gradient_calulation_type = 'all'
                print 'unable to calculate attack gradient with the \'mean\' method, reverting to \'all\' method.'

        if gradient_calulation_type == 'mean':
            # calculate the gradient based on the weighted mean of each attack
            threshold_step = dyn_range / (thresh_no + 2)

            gradient_thresh_array = np.arange(start_level, end_level + (threshold_step * dyn_range),
                                              (threshold_step * dyn_range))
            cross_threshold_times = np.zeros(len(gradient_thresh_array))
            cross_threshold_values = np.zeros(len(gradient_thresh_array))
            gradient_envelope_segment = envelope_samples[th_start_idx:th_end_idx + 1]

            for i in range(len(cross_threshold_values)):
                hold = np.argmax(gradient_envelope_segment >= gradient_thresh_array[i])
                cross_threshold_times[i] = hold[0] / float(fs)
                cross_threshold_values[i] = gradient_envelope_segment[hold[0]]

            pente_v = np.diff(cross_threshold_values) / np.diff(cross_threshold_times)

            # calculate weighted average of all gradients with a gausian dsitribution
            m_threshold = 0.5 * (gradient_thresh_array[:-1] + gradient_thresh_array[1:])
            weight_v = np.exp(-(m_threshold - 0.5) ** 2 / (0.5 ** 2))

            attack_gradient = np.sum(pente_v * weight_v) / np.sum(weight_v)

        elif gradient_calulation_type == 'all':
            # calculate the attack gradient from th_start_idx to th_end_idx
            attack_gradient = (end_level - start_level) / attack_time

        '''
          More stuff to return if we want extra information to be displayed
        '''
        thresholds_to_return = [calculation_type, th_start_idx + min_pre_peak_idx, th_end_idx + min_pre_peak_idx,
                                threshold_idxs + min_pre_peak_idx]

    elif calculation_type == 'fixed_threshold':
        # set threshold values for fixed threshold method
        fixed_threshold_start = 20
        fixed_threshold_end = 90

        # get dynamic range
        dyn_range = max(envelope_samples) - min(envelope_samples)

        # get thresholds relative to envelope level
        lower_threshold = (fixed_threshold_start * dyn_range * 0.01) + min(envelope_samples)
        upper_threshold = (fixed_threshold_end * dyn_range * 0.01) + min(envelope_samples)

        # calculate start index
        th_start_idx = np.argmax(envelope_samples >= lower_threshold)
        # th_start_idx = th_start_idx[0]

        # find the end idx after the start idx
        th_end_idx = np.argmax(envelope_samples[th_start_idx:] >= upper_threshold)
        th_end_idx = th_end_idx + th_start_idx

        if th_start_idx == th_end_idx:
            attack_time = 1.0 / fs
        else:
            attack_time = (th_end_idx - th_start_idx + 1.0) / fs

        # compare attack time to maximum permissible attack time
        if max_attack_time > 0:
            if attack_time > max_attack_time:
                # how many samples is equivalent to the maximum?
                max_attack_time_sample = int(fs * max_attack_time)  # convert to integer
                th_end_idx = th_start_idx + max_attack_time_sample
                attack_time = (th_end_idx - th_start_idx + 1.0) / fs

        # calculate the gradient

        # find the level of the first sample used
        start_level = envelope_samples[int(th_start_idx)]
        # find the level of the last sample used
        end_level = envelope_samples[int(th_end_idx)]

        # specify exceptions for a step functions crossing both thresholds
        if start_level == end_level:
            if th_start_idx > 0:
                # if a previous sample is avaiable, take the previous starting sample
                start_level = envelope_samples[int(th_start_idx) - 1]
            else:
                # set start level to zero if onset is at the first sample (indicating a step function at time zero)
                start_level = 0.0

        attack_gradient = (end_level - start_level) / attack_time

        '''
          More details to be returned if desired
        '''
        thresholds_to_return = [calculation_type, th_start_idx + min_pre_peak_idx, th_end_idx + min_pre_peak_idx]

    else:
        raise ValueError('calculation_type must be set to either \'fixed_threshold\' or \'min_effort\'.')

    # convert attack time to logarithmic scale
    attack_time = np.log10(attack_time)

    # revert attack gradient metric if envelope has been normalised
    if normalise:
        attack_gradient *= normalise_factor

    '''
      Calculate the temporal centroid
    '''
    hold_env = envelope_samples[int(th_start_idx):int(th_end_idx) + 1]
    t = np.arange(0, len(hold_env)) / float(fs)
    temp_centroid = np.sum(t * hold_env) / np.sum(hold_env)
    temp_centroid /= float(len(hold_env))

    if return_descriptive_data:
        return attack_time, attack_gradient, int(th_start_idx + min_pre_peak_idx), temp_centroid, thresholds_to_return
    else:
        return attack_time, attack_gradient, int(th_start_idx + min_pre_peak_idx), temp_centroid


def calculate_onsets(audio_samples, envelope_samples, fs, look_back_time=20, hysteresis_time=300, hysteresis_percent=10,
                     onset_in_noise_threshold=10, threshold_correction='onset_strength',
                     minimum_onset_time_separation=100, method='librosa', nperseg=512):
    '''
     Calculate onset idx using librosa

     I'm still adding new features to this so remember to update the definitions of all the input features
    '''
    if method == 'librosa':
        onsets = librosa.onset.onset_detect(audio_samples, fs, backtrack=True, units='samples')
    elif method == 'essentia_hfc':
        od1 = es.OnsetDetection(method='hfc')
        w = es.Windowing(type='hann')
        fft = es.FFT()  # this gives us a complex FFT
        c2p = es.CartesianToPolar()  # and this turns it into a pair (magnitude, phase)
        pool = Pool()
        spectrum = es.Spectrum()
        # let's get down to business
        for frame in es.FrameGenerator(audio_samples, frameSize=1024, hopSize=512):
            mag, phase, = c2p(fft(w(frame)))
            pool.add('features.hfc', od1(mag, phase))

        # Phase 2: compute the actual onsets locations
        hold_onsets = es.Onsets()
        onsets = hold_onsets(array([pool['features.hfc']]), [1])
        onsets *= fs

    elif method == 'essentia_complex':
        od2 = es.OnsetDetection(method='complex')
        # let's also get the other algorithms we will need, and a pool to store the results
        w = es.Windowing(type='hann')
        fft = es.FFT()  # this gives us a complex FFT
        c2p = es.CartesianToPolar()  # and this turns it into a pair (magnitude, phase)
        pool = Pool()
        spectrum = es.Spectrum()
        # let's get down to business
        for frame in es.FrameGenerator(audio_samples, frameSize=1024, hopSize=512):
            mag, phase, = c2p(fft(w(frame)))
            pool.add('features.complex', od2(mag, phase))

        # Phase 2: compute the actual onsets locations
        hold_onsets = es.Onsets()
        onsets = hold_onsets(array([pool['features.complex']]), [1])
        onsets *= fs
    else:
        raise ValueError('method for calculating onsets musy be \'librosa\', \'essentia_hfc\', or \'essentia_complex\'')

    # set values for return_loop method
    time_thresh = int(look_back_time * 0.001 * fs)  # 10 ms default look-back time, in samples
    hysteresis_samples = int(hysteresis_time * fs * 0.001)  # hysteresis time, in samples
    envelope_dyn_range = max(envelope_samples) - min(envelope_samples)
    hysteresis_thresh = envelope_dyn_range * hysteresis_percent * 0.01

    # only conduct analysis if there are onsets detected
    if np.size(onsets):
        # empty array for storing exact onset idxs
        corrected_onsets = []

        for onset_idx in onsets:
            # if the onset is 1 or 0, it's too close to the start to be corrected (1 is here due to zero padding)
            if onset_idx > 0:
                # actual onset location in samples (librosa uses 512 window size by default)
                onset_loc = onset_idx  # * 512
                onset_loc = np.array(onset_loc).astype('int')

                # only calculate if the onset is NOT at the end of the file, whilst oter onsets exist.
                # If the only onset is at the end, calculate anyway.
                if not corrected_onsets:
                    onset_hold = return_loop(onset_loc, envelope_samples, time_thresh, hysteresis_thresh,
                                             hysteresis_samples, nperseg=nperseg)
                    corrected_onsets.append(onset_hold)
                else:
                    if (onset_loc + 511) < len(envelope_samples):
                        onset_hold = return_loop(onset_loc, envelope_samples, time_thresh, hysteresis_thresh,
                                                 hysteresis_samples, nperseg=nperseg)
                        corrected_onsets.append(onset_hold)
            else:
                corrected_onsets.append(0)

        # zero is returned from return_loop if no valid onset identified
        # remove zeros (except the first)
        zero_loc = np.where(np.array(corrected_onsets) == 0)[0]
        # ignore if the first value is zero
        if list(zero_loc):
            if zero_loc[0] == 0:
                zero_loc = zero_loc[1:]
        corrected_onsets = np.delete(corrected_onsets, zero_loc)

        # remove duplicates
        hold_onsets = []
        for i in corrected_onsets:
            if i not in hold_onsets:
                hold_onsets.append(i)
        corrected_onsets = hold_onsets

        '''
         Methods of removing erronious onsets.  
         Method strings can be 'envelope_level_threshold', 'onset_strength', or 'none'
        '''
        # if threshold_correction == 'envelope_level_threshold':
        '''
         Remove repeated onsets and compare onset segments against the dynamic range
         to remove erroneous onsets in noise.  If the onset segment (samples between
         adjacent onsets) has a dynamic range less than 10% of total dynamic range,
         remove this onset.
        '''
        if len(corrected_onsets) > 1:
            thd_corrected_onsets = []
            last_value = corrected_onsets[-1]
            threshold = onset_in_noise_threshold * envelope_dyn_range * 0.01

            for i in reversed(range(len(corrected_onsets))):
                if corrected_onsets[i] == corrected_onsets[-1]:
                    segment = envelope_samples[corrected_onsets[i]:]
                else:
                    segment = envelope_samples[corrected_onsets[i]:corrected_onsets[i + 1]]

                # get the onset segment
                # segment = envelope_samples[corrected_onsets[i]:next_onset_idx]

                # only conduct if the segment if greater than 1 sample long
                if len(segment) > 1:
                    # find attack portion SNR
                    peak_idx = np.argmax(segment)
                    if peak_idx > 0:
                        # get the dynamic range of the attack portion
                        seg_dyn_range = max(segment) - min(segment[:peak_idx])
                        if seg_dyn_range >= threshold:
                            pass
                        else:
                            corrected_onsets = np.delete(corrected_onsets, i)
                    else:
                        corrected_onsets = np.delete(corrected_onsets, i)
                else:
                    corrected_onsets = np.delete(corrected_onsets, i)

        # remove onsets that are too close together, favouring the earlier onset
        if len(corrected_onsets) > 1:
            minimum_onset_time_separation_samples = fs * 0.001 * minimum_onset_time_separation
            time_separation = np.diff(corrected_onsets)
            # while loop for potential multiple itterations
            while len(corrected_onsets) > 1 and min(time_separation) < minimum_onset_time_separation_samples:
                onsets_to_remove = []
                # some onsets are closer together than the minimum value
                for i in range(len(corrected_onsets) - 1):
                    # are the last two onsets too close?
                    if abs(corrected_onsets[i + 1] - corrected_onsets[i]) < minimum_onset_time_separation_samples:
                        onsets_to_remove.append(i + 1)

                # remove onsets too close together
                corrected_onsets = np.delete(corrected_onsets, onsets_to_remove)
                time_separation = np.diff(corrected_onsets)

        # elif threshold_correction == 'onset_strength':
        '''
          Correct onsets by comparing to the onset strength.

          If there in an onset strength of 3 or greater between two onsets, then the onset if valid.  
          Otherwise, discard the onset.
        '''
        thd_corrected_onsets = []

        # get the onset strength
        onset_strength = librosa.onset.onset_strength(audio_samples, fs)

        strength_onset_times = np.array(corrected_onsets) / 512
        strength_onset_times.clip(min=0)

        corrected_original_onsets = []
        corrected_strength_onsets = []
        for onset_idx in reversed(range(len(corrected_onsets))):
            current_strength_onset = strength_onset_times[onset_idx]
            if current_strength_onset == strength_onset_times[-1]:
                onset_strength_seg = onset_strength[current_strength_onset:]
            else:
                onset_strength_seg = onset_strength[current_strength_onset:strength_onset_times[onset_idx + 1]]

            if max(onset_strength_seg) < 3:
                strength_onset_times = np.delete(strength_onset_times, onset_idx)
            else:
                thd_corrected_onsets.append(corrected_onsets[onset_idx])

        # elif threshold_correction == 'none':
        #     thd_corrected_onsets = corrected_onsets

        # insert method here for correction

    else:
        return []

    thd_corrected_onsets.sort()
    return thd_corrected_onsets


def get_bandwidth_array(audio_samples, fs, nperseg=512, overlap_step=32, rolloff_thresh=0.01,
                        rollon_thresh_percent=0.05, log_bandwidth=False, return_centroid=False,
                        low_bandwidth_method='Percentile', normalisation_method='RMS_Time_Window'):
    noverlap = nperseg - overlap_step
    # get spectrogram
    f, t, spec = spectrogram(audio_samples, fs, window='boxcar', nperseg=nperseg, noverlap=noverlap, scaling='density',
                             mode='magnitude')

    # normalise the spectrogram
    if normalisation_method == 'Single_TF_Bin':
        spec /= np.max(spec)
    elif normalisation_method == 'RMS_Time_Window':
        spec /= np.max(np.sqrt(np.sum(spec * spec, axis=0)))
    else:
        raise ValueError('Bandwidth normalisation method must be \'Single_TF_Bin\' or \'RMS_Time_Window\'')

    # initialise lists for storage
    rollon = []
    rolloff = []
    bandwidth = []
    centroid = []
    centroid_power = []

    # calculate the bandwidth curve
    for time_count in range(len(t)):
        seg = spec[:, time_count]
        tpower = np.sum(seg)
        if tpower > 0.0:
            if low_bandwidth_method == 'Percentile':
                # get the spectral rollon
                rollon_counter = 1
                cumulative_power = np.sum(seg[:rollon_counter])
                rollon_thresh = tpower * rollon_thresh_percent

                while cumulative_power < rollon_thresh:
                    rollon_counter += 1
                    cumulative_power = np.sum(seg[:rollon_counter])
                rollon.append(f[rollon_counter - 1])
            elif low_bandwidth_method == 'Cutoff':
                rollon_idx = np.where(seg >= rolloff_thresh)[0]
                if len(rollon_idx):
                    rollon_idx = rollon_idx[0]
                    rollon.append(f[rollon_idx])
            else:
                raise ValueError('low_bandwidth_method must be \'Percentile\' or \'Cutoff\'')

            # get the spectral rolloff
            rolloff_idx = np.where(seg >= rolloff_thresh)[0]
            if len(rolloff_idx):
                rolloff_idx = rolloff_idx[-1]
                rolloff.append(f[rolloff_idx])
                if log_bandwidth:
                    bandwidth.append(np.log(f[rolloff_idx] / float(f[rollon_counter - 1])))
                else:
                    bandwidth.append(f[rolloff_idx] - f[rollon_counter - 1])
            else:
                bandwidth.append(0)
        else:
            bandwidth.append(0)

        if tpower > 0.05:
            centroid.append(np.sum(seg * f) / np.sum(seg))
            centroid_power.append(tpower)
    if return_centroid:
        return bandwidth, t, f, np.average(centroid, weights=centroid_power)
    else:
        return bandwidth, t, f


def calculate_bandwidth_gradient(bandwidth_segment, t):
    if bandwidth_segment:
        max_idx = np.argmax(bandwidth_segment)
        if max_idx > 0:
            min_idx = np.where(np.array(bandwidth_segment[:max_idx]) == min(bandwidth_segment[:max_idx]))[0][-1]

            bandwidth_change = bandwidth_segment[max_idx] - bandwidth_segment[min_idx]
            time_to_change = (max_idx - min_idx) * (t[1] - t[0])

            bandwidth_gradient = bandwidth_change / time_to_change
        else:
            bandwidth_gradient = False
    else:
        bandwidth_gradient = False
    return bandwidth_gradient


def timbral_hardness(fname, dev_output=False, max_attack_time=0.3, bandwidth_thresh_db=-40, weighting='None',
                     phase_correction=False):
    """
     This function calculates the apparent hardness of an audio file.
     Currently the function only works on the left channel of multi-channel audio files.

      :param fname: Audio filename to be analysed, including full file path and extension.
      :return: returns the hardness as a float representing the relative hardness.
    """
    import scipy.stats
    # max_attack_time = 0.3  # limit the maximum attack time of all calculations to 0.3 seconds, set to -1 to ignore
    # nperseg = 4096
    nperseg = 4096
    # use pysoundfile to read audio
    audio_samples, fs = sf.read(fname, always_2d=False)

    # get sum all channels to mono
    num_channels = np.shape(audio_samples)
    if len(num_channels) > 1:
        # crudely check for out of phase signals
        if phase_correction:
            r, pval = scipy.stats.pearsonr(audio_samples[:, 0], audio_samples[:, 1])
            if r < -0.5:
                audio_samples = audio_samples[:, 0]  # [:,1] *= -1.0
            else:
                audio_samples = np.sum(audio_samples, axis=1)
        else:
            audio_samples = np.sum(audio_samples, axis=1)

    if weighting != 'None':
        audio_samples = pyfilterbank.splweighting.weight_signal(audio_samples, fs, weighting)

    # normalise audio
    audio_samples /= max(abs(audio_samples))

    # get percussive ratio
    HP_ratio = get_percussive_audio(audio_samples, show_ratio=True)

    # zero pad the signal
    audio_samples = np.lib.pad(audio_samples, (nperseg + 1, 0), 'constant', constant_values=(0.0, 0.0))

    envelope = sample_and_hold_envelope_calculation(audio_samples, fs, decay_time=0.1)
    envelope_time = np.arange(len(envelope)) / float(fs)
    # calculate the onsets
    original_onsets = calculate_onsets(audio_samples, envelope, fs, nperseg=nperseg)
    onset_strength = librosa.onset.onset_strength(audio_samples, fs)

    '''
      Calculate the spectrogram so that the bandwidth can be created
    '''

    mag = db2mag(bandwidth_thresh_db)
    bandwidth, t, f, centroid = get_bandwidth_array(audio_samples, fs, nperseg=nperseg, overlap_step=32,
                                                    rolloff_thresh=mag, return_centroid=True)
    # centroid = np.mean(centroid)
    all_bandwidth_attack_time = []
    all_bandwidth_attack_gradient = []
    all_bandwidth_attack_temp_centroid = []
    all_bandwidth_max = []
    all_attack_time = []
    all_attack_gradient = []
    all_attack_temp_centroid = []
    all_const_thresh_attack_time = []
    all_const_thresh_attack_gradient = []
    all_const_thresh_attack_temp_centroid = []
    all_max_strength = []
    all_attack_centroid = []
    all_attack_unitless_centroid = []
    all_max_strength_bandwidth = []

    '''
      Get bandwidth onset times and max bandwidth
    '''
    # If onsets don't exist, set it to time zero
    if not original_onsets:
        original_onsets = [0]

    # If onsets exist, convert the onsets to time domain
    if original_onsets:
        onsets = np.array(original_onsets) - nperseg  # 512 # nperseg=512
        onsets[onsets < 0] = 0
        bandwidth_onset = np.array(onsets) / 32  # overlap_step=32

        for onset_count in range(len(bandwidth_onset)):
            onset = bandwidth_onset[onset_count]
            if onset == bandwidth_onset[-1]:
                bandwidth_seg = np.array(bandwidth[onset:])
            else:
                next_onset = bandwidth_onset[onset_count + 1]
                bandwidth_seg = np.array(bandwidth[onset:next_onset])

            '''
              bandwidth maximum needs to be dependent on the attack time, however, the attack time calculation 
              function is causing the values to be screwed, due to normalisation of the array I think.  
              I'll look into this now.
            '''
            # making a copy of the bandqwidth segment
            hold_bandwidth_seg = list(bandwidth_seg)

            '''
              I need a different sample rate for this, because the bandwidth goes at the same rate as the 
              spectrogram, not the audio sample rate.

              So the bandwidth_rate will be the sample rate devided by the overlap step.  Then the bandwidth 
              sample rate will be in seconds.
            '''

            bandwidth_fs = fs / float(32)  # overlap_step=32
            # bandwidth_attack = calculate_attack_time(bandwidth_seg, bandwidth_fs, return_descriptive_data=True)
            if max(bandwidth_seg) > 0:
                bandwidth_attack = calculate_attack_time(bandwidth_seg, bandwidth_fs,
                                                         calculation_type='fixed_threshold',
                                                         max_attack_time=max_attack_time)
            else:
                bandwidth_attack = []

            # attack_time, attack_gradient, th_start_idx, temp_centroid
            if bandwidth_attack:
                all_bandwidth_attack_time.append(bandwidth_attack[0])
                all_bandwidth_attack_gradient.append(bandwidth_attack[1])
                all_bandwidth_attack_temp_centroid.append(bandwidth_attack[3])
                start_idx = bandwidth_attack[2]
                if max_attack_time > 0:
                    max_attack_time_samples = int(max_attack_time * bandwidth_fs)
                    if len(hold_bandwidth_seg[start_idx:]) > start_idx + max_attack_time_samples:
                        all_bandwidth_max.append(max(hold_bandwidth_seg[start_idx:start_idx + max_attack_time_samples]))
                    else:
                        all_bandwidth_max.append(max(hold_bandwidth_seg[start_idx:]))
                else:
                    all_bandwidth_max.append(max(hold_bandwidth_seg[start_idx:]))

            '''
              Get the attack time stuff
            '''
            onset = original_onsets[onset_count]
            if onset == original_onsets[-1]:
                attack_seg = np.array(envelope[onset:])
                strength_seg = np.array(onset_strength[(onset / 512):])
                audio_seg = np.array(audio_samples[onset:])
            else:
                attack_seg = np.array(envelope[onset:original_onsets[onset_count + 1]])
                strength_seg = np.array(onset_strength[(onset / 512):(original_onsets[onset_count + 1] / 512)])
                audio_seg = np.array(audio_samples[onset:original_onsets[onset_count + 1]])

            attack_time = calculate_attack_time(attack_seg, fs, max_attack_time=max_attack_time)

            all_attack_time.append(attack_time[0])
            all_attack_gradient.append(attack_time[1])
            all_attack_temp_centroid.append(attack_time[3])
            th_start_idx = attack_time[2]

            attack_time = calculate_attack_time(attack_seg, fs, calculation_type='fixed_threshold',
                                                max_attack_time=max_attack_time)
            all_const_thresh_attack_time.append(attack_time[0])
            all_const_thresh_attack_gradient.append(attack_time[1])
            all_const_thresh_attack_temp_centroid.append(attack_time[3])

            all_max_strength.append(max(strength_seg))
            if bandwidth_attack:
                all_max_strength_bandwidth.append(max(strength_seg))

            '''
              Get the spectral centroid of the attack
            '''
            # identify the start of the attack
            centroid_int_time = 0.125  # 125ms maximum integration time for spectral centroid
            centroid_int_samples = int(centroid_int_time * fs)  # number of samples for attack time integration

            # find start of attack section
            max_idx = np.argmax(attack_seg)
            min_idx = np.where(attack_seg[:max_idx] == min(attack_seg[:max_idx]))[-1][-1]

            if th_start_idx + centroid_int_samples >= len(audio_seg):
                audio_seg = audio_seg[th_start_idx:]
            else:
                audio_seg = audio_seg[th_start_idx:th_start_idx + centroid_int_samples]

            spectral_features_hold = get_spectral_features(audio_seg, fs)

            if spectral_features_hold:
                all_attack_centroid.append(spectral_features_hold[0])
                all_attack_unitless_centroid.append(spectral_features_hold[2])

        mean_bandwidth_attack_time = np.mean(all_bandwidth_attack_time)
        mean_bandwidth_attack_gradient = np.mean(all_bandwidth_attack_gradient)
        mean_bandwidth_attack_temp_centroid = np.mean(all_bandwidth_attack_temp_centroid)
        mean_bandwidth_max = np.mean(all_bandwidth_max)
        mean_attack_time = np.mean(all_attack_time)
        mean_attack_gradient = np.mean(all_attack_gradient)
        mean_attack_temp_centroid = np.mean(all_attack_temp_centroid)
        mean_const_thresh_attack_time = np.mean(all_const_thresh_attack_time)
        mean_const_thresh_attack_gradient = np.mean(all_const_thresh_attack_gradient)
        mean_const_thresh_attack_temp_centroid = np.mean(all_const_thresh_attack_temp_centroid)
        mean_max_strength = np.mean(all_max_strength)
        mean_attack_centroid = np.mean(all_attack_centroid)
        mean_attack_unitless_centroid = np.mean(all_attack_unitless_centroid)

        # weight all metrics prior to taking the mean
        mean_weighted_bandwidth_attack_time = np.average(all_bandwidth_attack_time, weights=all_max_strength_bandwidth)
        mean_weighted_bandwidth_attack_gradient = np.average(all_bandwidth_attack_gradient,
                                                             weights=all_max_strength_bandwidth)
        mean_weighted_bandwidth_attack_temp_centroid = np.average(all_bandwidth_attack_temp_centroid,
                                                                  weights=all_max_strength_bandwidth)
        mean_weighted_bandwidth_max = np.average(all_bandwidth_max, weights=all_max_strength_bandwidth)
        mean_weighted_attack_time = np.average(all_attack_time, weights=all_max_strength)
        mean_weighted_attack_gradient = np.average(all_attack_gradient, weights=all_max_strength)
        mean_weighted_attack_temp_centroid = np.average(all_attack_temp_centroid, weights=all_max_strength)
        mean_weighted_const_thresh_attack_time = np.average(all_const_thresh_attack_time, weights=all_max_strength)
        mean_weighted_const_thresh_attack_gradient = np.average(all_const_thresh_attack_gradient,
                                                                weights=all_max_strength)
        mean_weighted_const_thresh_attack_temp_centroid = np.average(all_const_thresh_attack_temp_centroid,
                                                                     weights=all_max_strength)
        mean_weighted_attack_centroid = np.average(all_attack_centroid, weights=all_max_strength)
        mean_weighted_attack_unitless_centroid = np.average(all_attack_unitless_centroid, weights=all_max_strength)

    else:

        mean_bandwidth_attack_time = 0
        mean_bandwidth_attack_gradient = 0
        mean_bandwidth_attack_temp_centroid = 0
        mean_bandwidth_max = 0
        mean_attack_time = 0
        mean_attack_gradient = 0
        mean_attack_temp_centroid = 0
        mean_const_thresh_attack_time = 0
        mean_const_thresh_attack_gradient = 0
        mean_const_thresh_attack_temp_centroid = 0
        mean_max_strength = 0
        mean_attack_centroid = 0
        mean_attack_unitless_centroid = 0

        # weight all metrics prior to taking the mean
        mean_weighted_bandwidth_attack_time = 0
        mean_weighted_bandwidth_attack_gradient = 0
        mean_weighted_bandwidth_attack_temp_centroid = 0
        mean_weighted_bandwidth_max = 0
        mean_weighted_attack_time = 0
        mean_weighted_attack_gradient = 0
        mean_weighted_attack_temp_centroid = 0
        mean_weighted_const_thresh_attack_time = 0
        mean_weighted_const_thresh_attack_gradient = 0
        mean_weighted_const_thresh_attack_temp_centroid = 0
        mean_weighted_attack_centroid = 0
        mean_weighted_attack_unitless_centroid = 0

    # return the values
    if dev_output:
        return mean_bandwidth_attack_time, mean_bandwidth_attack_gradient, mean_bandwidth_attack_temp_centroid, \
               mean_bandwidth_max, mean_attack_time, mean_attack_gradient, mean_attack_temp_centroid, \
               mean_const_thresh_attack_time, mean_const_thresh_attack_gradient, mean_const_thresh_attack_temp_centroid, \
               mean_max_strength, mean_attack_centroid, mean_attack_unitless_centroid, \
               mean_weighted_bandwidth_attack_time, mean_weighted_bandwidth_attack_gradient, \
               mean_weighted_bandwidth_attack_temp_centroid, mean_weighted_bandwidth_max, \
               mean_weighted_attack_time, mean_weighted_attack_gradient, mean_weighted_attack_temp_centroid, \
               mean_weighted_const_thresh_attack_time, mean_weighted_const_thresh_attack_gradient, \
               mean_weighted_const_thresh_attack_temp_centroid, mean_weighted_attack_centroid, \
               mean_weighted_attack_unitless_centroid, centroid, HP_ratio
    else:
        '''
         Apply regression model
        '''

        all_metrics = np.ones(8)

        all_metrics[0] = mean_attack_time
        all_metrics[1] = mean_attack_gradient
        all_metrics[2] = mean_attack_centroid
        all_metrics[3] = np.array(mean_attack_time) * np.array(mean_attack_gradient)
        all_metrics[4] = np.array(mean_attack_time) * np.array(mean_attack_centroid)
        all_metrics[5] = np.array(mean_attack_gradient) * np.array(mean_attack_centroid)
        all_metrics[6] = np.array(mean_attack_time) * np.array(mean_attack_gradient) * np.array(mean_attack_centroid)

        coefficients = np.array([-365.232930949, -3276.19905108, 426.807386961, -624.176789287, 104.143569144,
                                 988.622652051, 189.821027679, -1407.5701818])

        hardness = np.sum(all_metrics * coefficients)
        return hardness

    # percussive_audio = get_percussive_audio(audio_samples)
    # # audio_samples = percussive_audio
    #
    #
    # '''
    #  Calculate the envelope of the audio with a 'sample and hold' style function.
    #  This ensures that the minimum attack time is not limited by low-pass filtering,
    #  a common method of obtaining the envelope.
    # '''
    # decay_time = 0.2
    # envelope = sample_and_hold_envelope_calculation(audio_samples, fs)
    #
    # # prepend 512 zeros to allow for initial onset to be detected
    # envelope = np.lib.pad(envelope, (512, 0), 'constant', constant_values=(0, 0))
    # audio_samples = np.lib.pad(audio_samples, (512, 0), 'constant', constant_values=(0, 0))
    #
    # # get bandwidth
    # nperseg = 512
    # overlap_step = 32
    # bandwidth, t, f = get_bandwidth_array(audio_samples, fs, nperseg=512, overlap_step=32)
    #
    # thd_corrected_onsets = calculate_onsets(audio_samples, envelope, fs, look_back_time=10,
    #                                         hysteresis_time=200, hysteresis_percent=5, onset_in_noise_threshold=10)
    #
    # onset_strength = librosa.onset.onset_strength(audio_samples, fs)
    # # onset_strength = 1
    # # get the onsets in the format for the bandwidth (since bandwidth is from a spectrogram, doesn't have the same time scale)
    # bandwidth_onsets = np.array(thd_corrected_onsets) - nperseg
    # bandwidth_onsets[bandwidth_onsets<0] = 0
    # bandwidth_onsets = np.array(bandwidth_onsets) / overlap_step
    #
    #
    # '''
    #  Preallocate arrays for storing variables
    # '''
    # # feature storage
    # all_attack_time = []
    # all_attack_gradient = []
    # all_temporal_centroid = []
    # all_attack_centroid = []
    # all_attack_spread = []
    # # all_spread = []
    # # all_spread_difference = []
    # all_unitless_attack_centroid = []
    # all_highpass_attack_time = []
    # all_highpass_attack_gradient = []
    # all_highpass_temporal_centroid = []
    # all_perceptual_centroid = []
    # all_perceptual_spread = []
    # all_perceptual_unitless_centroid = []
    # all_spread_gradient = []
    # all_octaveband_attack_time = []
    # all_octaveband_attack_gradient = []
    # all_octaveband_temporal_centroid = []
    # all_bandwidth_gradient = []
    #
    # centre_frequencies = [125, 250, 500, 1000, 2000, 4000, 8000]
    # for band_no in range(len(centre_frequencies)):
    #     all_octaveband_attack_time.append([])
    #     all_octaveband_attack_gradient.append([])
    #     all_octaveband_temporal_centroid.append([])
    #
    #
    # '''
    #   Filter the audio into octave bands
    # '''
    # octave_band_envelope = []
    # # plotcounter = 3
    # for f0 in centre_frequencies:
    #     hold_audio = filter_audio_bandpass(audio_samples, f0, noct=1, fs=fs, order=2)
    #     hold_envelope = sample_and_hold_envelope_calculation(hold_audio, fs)
    #     octave_band_envelope.append(hold_envelope)
    #
    # '''
    #   Highpass filter the audio and calculate the attack time
    # '''
    # crossover_frequency = 125
    # hold_audio = filter_audio_highpass(audio_samples, crossover_frequency, fs, order=2)
    # highpass_envelope = sample_and_hold_envelope_calculation(hold_audio, fs)
    #
    #
    # '''
    #  Calculate the metrics for each onset.
    # '''
    # if thd_corrected_onsets[0] >= 0:
    #     thresh_no = 8  # number of thresholds to calculate
    #     centroid_int_time = 0.125  # 125ms maximum integration time for spectral centroid
    #     centroid_int_samples = int(centroid_int_time * fs)  # number of samples for attack time integration
    #
    #     for i in range(len(thd_corrected_onsets)):
    #         # get the segment containing the attack, both the envelope and audio
    #         start_idx = (thd_corrected_onsets[i])
    #         # if it's the last element
    #         if thd_corrected_onsets[i] == thd_corrected_onsets[-1]:
    #             segment = envelope[start_idx:]
    #             highpass_segment = highpass_envelope[start_idx:]
    #             audio_segment = audio_samples[start_idx:]
    #             bandwidth_segment = bandwidth[bandwidth_onsets[i]:]
    #         else:
    #             end_idx = thd_corrected_onsets[i + 1]
    #             segment = envelope[start_idx:end_idx]
    #             highpass_segment = highpass_envelope[start_idx:end_idx]
    #             audio_segment = audio_samples[start_idx:end_idx]
    #             bandwidth_segment = bandwidth[bandwidth_onsets[i]:bandwidth_onsets[i+1]]
    #
    #         bandwidth_gradient = calculate_bandwidth_gradient(bandwidth_segment, t)
    #         if bandwidth_gradient != False:
    #             all_bandwidth_gradient.append(bandwidth_gradient)
    #
    #         attack_metrics = calculate_attack_time(segment, fs)
    #         if attack_metrics:
    #             all_attack_time.append(attack_metrics[0])
    #             all_attack_gradient.append(attack_metrics[1])
    #             th_start_idx = attack_metrics[2]
    #             all_temporal_centroid.append(attack_metrics[3])
    #
    #             '''
    #               Calculate the spectral centroid of the attack
    #             '''
    #             if (centroid_int_samples + th_start_idx) <= len(segment):
    #                 post_attack_segment = audio_segment[int(th_start_idx):int(th_start_idx + centroid_int_samples)]
    #             else:
    #                 post_attack_segment = audio_segment[int(th_start_idx):len(segment)]
    #
    #             spectral_features_hold = get_spectral_features(post_attack_segment, fs)
    #             if spectral_features_hold != 0:
    #                 all_attack_centroid.append(np.log10(spectral_features_hold[0]))
    #                 all_attack_spread.append(spectral_features_hold[1])
    #                 all_unitless_attack_centroid.append(spectral_features_hold[0])
    #
    #         '''
    #           Calculate highpass attack features
    #         '''
    #         highpass_attack_metrics = calculate_attack_time(highpass_segment, fs)
    #         if highpass_attack_metrics:
    #             all_highpass_attack_time.append(highpass_attack_metrics[0])
    #             all_highpass_attack_gradient.append(highpass_attack_metrics[1])
    #             all_highpass_temporal_centroid.append(highpass_attack_metrics[3])
    #
    #         '''
    #           From here until the end of the loop, new metrics will be calculated.
    #           These include:
    #             perceptual spectral metrics
    #
    #             attack spectral slope - yet to be implemented
    #
    #             bandwidth gradient
    #         '''
    #
    #         '''
    #           Perceptual centroid, spread, and unitless centroid
    #         '''
    #         perceptual_spectral_features_hold = get_spectral_features(post_attack_segment, fs, scale='mel')
    #         if perceptual_spectral_features_hold != 0:
    #             all_perceptual_centroid.append(perceptual_spectral_features_hold[0])
    #             all_perceptual_spread.append(perceptual_spectral_features_hold[0])
    #             all_perceptual_unitless_centroid.append(perceptual_spectral_features_hold[0])
    #
    #         '''
    #           bandwidth gradient
    #         '''
    #         # get a spectrogram of the signal
    #         f, t, spec = spectrogram(audio_segment,fs)
    #         spread_over_time = []
    #
    #         for t_seg in range(len(t)):
    #             spec_hold = spec[:,t_seg]
    #
    #             if np.sum(spec_hold) == 0:
    #                 spread_hold = 0
    #             else:
    #                 hold_centroid = np.sum(spec_hold * f) / np.sum(spec_hold)
    #                 deviation = np.abs(f - hold_centroid)
    #                 spread_hold = np.sqrt(np.sum(spec_hold * (deviation ** 2)))
    #
    #             spread_over_time.append(spread_hold)
    #
    #
    #         max_idx = np.argmax(spread_over_time)
    #         if max_idx > 0:
    #             min_idx = np.argmin(spread_over_time[:max_idx+1])
    #         else:
    #             min_idx = 0
    #
    #         time_change = t[max_idx] - t[min_idx]
    #         spread_change = spread_over_time[max_idx] - spread_over_time[min_idx]
    #
    #         all_spread_gradient.append(spread_change/time_change)
    #
    #
    #         '''
    #           Octave band attack time
    #         '''
    #         for band_no in range(len(octave_band_envelope)):
    #             octave_env = octave_band_envelope[band_no]
    #
    #             # if it's the last element
    #             if thd_corrected_onsets[i] == thd_corrected_onsets[-1]:
    #                 octave_env = octave_env[start_idx:]
    #             else:
    #                 octave_env = octave_env[start_idx:end_idx]
    #
    #             attack_metrics = calculate_attack_time(octave_env, fs)
    #             if attack_metrics:
    #                 all_octaveband_attack_time[band_no].append(attack_metrics[0])
    #                 all_octaveband_attack_gradient[band_no].append(attack_metrics[1])
    #                 all_octaveband_temporal_centroid[band_no].append(attack_metrics[3])
    #
    # '''
    #  Get mean values of metrics for all onsets
    # '''
    # # all_spread_gradient = []
    # # all_octaveband_attack_time = []
    # # all_octaveband_attack_gradient = []
    # # all_octaveband_temporal_centroid = []
    # # for band_no in range(len(centre_frequencies)):
    # #     all_octaveband_attack_time.append([])
    # #     all_octaveband_attack_gradient.append([])
    # #     all_octaveband_temporal_centroid.append([])
    #
    #
    #
    #
    #
    #
    #
    # mean_attack_time = np.nanmean(all_attack_time)
    # mean_attack_gradient = np.nanmean(all_attack_gradient)
    # mean_temporal_centroid = np.nanmean(all_temporal_centroid)
    #
    # mean_attack_centroid = np.nanmean(all_attack_centroid)
    # mean_attack_spread = np.mean(all_attack_spread)
    # mean_unitless_attack_centroid = np.mean(all_unitless_attack_centroid)
    #
    # mean_highpass_attack_time = np.mean(all_highpass_attack_time)
    # mean_highpass_attack_gradient = np.mean(all_highpass_attack_gradient)
    # mean_highpass_temporal_centroid = np.mean(all_highpass_temporal_centroid)
    #
    # mean_perceptual_centroid = np.mean(all_perceptual_centroid)
    # mean_perceptual_spread = np.mean(all_perceptual_spread)
    # mean_perceptual_unitless_centroid = np.nanmean(all_perceptual_unitless_centroid)
    #
    # mean_spread_gradient = np.nanmean(all_spread_gradient)
    #
    # mean_octaveband_attack_time = []
    # mean_octaveband_attack_gradient = []
    # mean_octaveband_temporal_centroid = []
    #
    # for band_no in range(len(all_octaveband_attack_time)):
    #     mean_octaveband_attack_time.append(np.nanmean(all_octaveband_attack_time[band_no]))
    #     mean_octaveband_attack_gradient.append(np.nanmean(all_octaveband_attack_gradient[band_no]))
    #     mean_octaveband_temporal_centroid.append(np.nanmean(all_octaveband_temporal_centroid[band_no]))
    #
    # mean_bandwidth_gradient = np.mean(all_bandwidth_gradient)
    # mean_onset_strength = np.mean(onset_strength)
    #
    # '''
    #  Apply regression model
    # '''
    #
    # all_metrics = np.ones(8)
    #
    # all_metrics[0] = mean_attack_time
    # all_metrics[1] = mean_attack_gradient
    # all_metrics[2] = mean_attack_centroid
    # all_metrics[3] = np.array(mean_attack_time) * np.array(mean_attack_gradient)
    # all_metrics[4] = np.array(mean_attack_time) * np.array(mean_attack_centroid)
    # all_metrics[5] = np.array(mean_attack_gradient) * np.array(mean_attack_centroid)
    # all_metrics[6] = np.array(mean_attack_time) * np.array(mean_attack_gradient) * np.array(mean_attack_centroid)
    #
    # coefficients = np.array([-365.232930949, -3276.19905108, 426.807386961, -624.176789287, 104.143569144,
    #                          988.622652051, 189.821027679, -1407.5701818])
    #
    # hardness = np.sum(all_metrics * coefficients)
    #
    # # return hardness
    # # this is for testing the function
    # return mean_attack_time, mean_attack_gradient, mean_temporal_centroid, mean_attack_centroid, mean_attack_spread, \
    #        mean_unitless_attack_centroid, mean_highpass_attack_time, mean_highpass_attack_gradient, \
    #        mean_highpass_temporal_centroid, mean_perceptual_centroid, mean_perceptual_spread, \
    #        mean_perceptual_unitless_centroid, mean_spread_gradient, mean_octaveband_attack_time, \
    #        mean_octaveband_attack_gradient, mean_octaveband_temporal_centroid, mean_bandwidth_gradient, mean_onset_strength
    #
    # old return version
    # return mean_attack_time, mean_attack_gradient, mean_attack_centroid, hardness, mean_attack_bandwidth, \
    #        mean_all_spread, mean_bandwidth_difference, mean_unitless_attack_centroid, mean_highpass_all_attack_time