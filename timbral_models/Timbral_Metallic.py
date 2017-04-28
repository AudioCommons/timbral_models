import numpy as np
from scipy.signal import hilbert, butter, lfilter
import Timbral_Roughness as rough
import soundfile as sf
import librosa


def butter_lowpass(cutoff, fs, order=2):
    """
     This function calculates the butterworth filter coefficients for a given cutoff frequency and order.

     :param cutoff: Cutoff frequency as a proportion of the Nyquist frequency (Freq (Hz) / (Sample rate / 2))
     :param fs:     Sample rate of signal to be filtered
     :param order:  Filter order, defaults to second order filter, as required by timbral_metallic
     :return: returns the coefficients of a Butterworth filter
    """
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=2):
    """
     This function designs a Butterworth lowpass filter and applies it to the data 
    
    :param data:    Audio data to be lowpass filtered
    :param cutoff:  Cutoff frequency in Hz
    :param fs:      Samplerate of audio
    :param order:   Order of the filter, defaults to second order filter as required by timbral_metallic
    :return:        Returns the filtered signal
    """
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y


def return_loop(onset_loc, envelope, function_time_thresh, hist_threshold, hist_time_samples):
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
        current_sample = envelope[onset_loc]
        # get the previous 10ms worth of samples
        if onset_loc - function_time_thresh >= 0:
            evaluation_array = envelope[onset_loc - function_time_thresh - 1:onset_loc - 1]
        else:
            evaluation_array = envelope[:onset_loc - 1]

        if min(evaluation_array) - current_sample < 0:
            '''
             If the minimum value within previous 10ms is less than current sample,
             move to the start position to the minimum value and look again.
            '''
            min_idx = np.argmin(evaluation_array)
            new_onset_loc = min_idx + onset_loc - function_time_thresh - 1

            if new_onset_loc > 512:
                onset_loc = new_onset_loc
            else:
                ''' Current index is close to start of the envelope, so exit with the idx as 0 '''
                return 0

        else:
            '''
             If the minimum value within previous 10ms is greater than current sample,
             introduce the time and level hysteresis to check again.
            '''
            # get the array of 200ms previous to the current onset idx
            if (onset_loc - hist_time_samples - 1) > 0:
                hyst_evaluation_array = envelope[onset_loc - hist_time_samples - 1:onset_loc - 1]
            else:
                hyst_evaluation_array = envelope[:onset_loc - 1]

            # values less than the current sample
            all_match = np.where(hyst_evaluation_array < envelope[onset_loc])

            # if no minimum was found within the extended time, exit with current onset idx
            if len(all_match[0]) == 0:
                return onset_loc

            # get the idx of the closest value which is lower than the current onset idx
            last_min = all_match[0][-1]
            last_idx = int(onset_loc - len(hyst_evaluation_array) + last_min - 1)

            # get the dynamic range
            segment_dynamic_range = max(hyst_evaluation_array[last_min:]) - current_sample

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
                onset_loc = last_idx


def get_attack_envelope(time_samples, fs, decay_time=0.2, hold_time=0.01):
    """
     This is the causal function implemented from Timbral_Hardness, returning the envelope of the signal.
     
    :param time_samples:    Audio array 
    :param fs:              Sample rate of audio file
    :param decay_time:      Decay time for the causal function, set by default to 200ms
    :param hold_time:       Hold time of the causal function, set by default to 100ms. 
    :return:                Approximated envelope of the audio array
    """
    abs_samples = abs(time_samples)
    attack_envelope = []

    # set parameters for envelope function
    decay = max(abs_samples) / (decay_time * fs)  # decay rate relative to peak level of audio signal
    hold_samples = hold_time * fs  # number of samples to hold before decay
    hold_counter = 0
    previous_sample = 0.0

    # perform the sample, hold, and decay function to obtain envelope
    for sample in abs_samples:
        if sample >= previous_sample:
            attack_envelope.append(sample)
            previous_sample = sample
            hold_counter = 0
        else:
            # check hold length
            if hold_counter < hold_samples:
                hold_counter += 1
                attack_envelope.append(previous_sample)
            else:
                out = previous_sample - decay
                if out > sample:
                    attack_envelope.append(out)
                    previous_sample = out
                else:
                    attack_envelope.append(sample)
                    previous_sample = sample

    # convert to numpy array
    attack_envelope = np.array(attack_envelope)
    return attack_envelope


def get_spectral_features(audio, fs, lf_limit=20):
    """
     This function calculates the spectral centroid and spectral spread of an audio array.
     
     :param audio:      Audio array 
     :param fs:         Sample rate of audio file
     :param lf_limit:   Low frequency limit, in Hz, to be analysed.  Defaults to 20Hz.
     :return:           Returns the spectral centroid and spectral spread
    """
    # use a hanning window
    window = np.hanning(len(audio))
    next_pow_2 = int(pow(2, np.ceil(np.log2(len(window)))))
    # get frequency domain representation
    spectrum = np.fft.fft((window * audio), next_pow_2)
    spectrum = np.absolute(spectrum[0:int(len(spectrum) / 2) + 1])
    freq = np.arange(0, len(spectrum), 1) * (fs / (2.0 * (len(spectrum) - 1)))

    # find lowest frequency index, zeros used to unpack result
    lf_limit_idx = np.where(freq >= lf_limit)[0][0]
    spectrum = spectrum[lf_limit_idx:]
    freq = freq[lf_limit_idx:]

    # calculate centroid and spread
    centroid = sum(spectrum * freq) / float(sum(spectrum))
    spread = np.sqrt(sum(((freq - centroid) ** 2) * spectrum) / sum(spectrum))

    return centroid, spread


def timbral_metallic(fname):
    """
     This function calculates the probability of an audio file sounding metallic.
     Currently the function only works on the left channel of multi-channel audio files.

     :param fname: Audio filename to be analysed, including full file path and extension.
     :return: returns the a probability, from 0.0 to 1.0, representing the probability of a signal sounding metallic.
    """
    # use pysoundfile to read audio
    audio_samples, fs = sf.read(fname, always_2d=False)

    # get left channel only
    num_channels = np.shape(audio_samples)
    if len(num_channels) > 1:
        audio_samples = audio_samples[:, 0]

    # normalise the audio
    audio_samples *= 1.0 / max(abs(audio_samples))

    # calculate the roughness of the audio file
    roughness = rough.timbral_roughness(fname)

    '''
     Get the envelope of the audio file using the Hilbert transform and filtering method
    '''
    # get analytical signal
    h = np.absolute(hilbert(audio_samples))

    # Filter requirements.
    order = 2
    cutoff = 50  # desired cutoff frequency of the filter, Hz

    # low-pass filter the signal
    envelope = butter_lowpass_filter(h, cutoff, fs, order)

    # normalise the envelope
    envelope *= 1.0 / max(envelope)

    '''
     Get the envelope using the Timbral_Hardness method
    '''
    decay_time = 0.2  # in seconds, 200 ms
    hold_time = 0.01  # time to hold samples before decay, in seconds
    attack_env = get_attack_envelope(audio_samples, fs, decay_time, hold_time)

    # normalise the attack envelope
    attack_env = attack_env / max(attack_env)

    '''
     Zero-pad audio and envelopes with 512 samples 
    '''
    # prepend 512 zeros to allow for initial onset to be detected
    audio_samples = np.lib.pad(audio_samples, (512, 0), 'constant', constant_values=(0, 0))
    attack_env = np.lib.pad(attack_env, (512, 0), 'constant', constant_values=(0, 0))
    envelope = np.lib.pad(envelope, (512, 0), 'constant', constant_values=(0, 0))
    logenvelope = np.log10(envelope+1.0)  # +1 to prevent log(0)
    logenvelope = logenvelope / max(logenvelope)  # re-normalise the envelope

    # calculate the onsets using librosa
    onsets = librosa.onset.onset_detect(audio_samples, fs)

    '''
     Get correct onset locations and assess against the dynamic range
    '''
    # set defaults for assessing the onsets
    corrected_onsets = []
    time_thresh = int(0.01 * fs)  # 10 ms time threshold for backwards analysis of rising
    hist_percent = 0.05  # percentage of dynamic range to allow variation
    hist_thresh = (max(attack_env) - min(attack_env)) * hist_percent
    hysteresis_samples = int(0.2 * fs)  # number of samples for 200ms

    # get more precise onset locations
    for onset_idx in onsets:
        if onset_idx <= 1:
            # onset is too close to the start on the audio file, so set as the start
            corrected_onsets.append(0)
        else:
            onset_loc = onset_idx * 512  # actual onset location in samples (since librosa uses 512 samples by default)
            onset_hold = return_loop(onset_loc, attack_env, time_thresh, hist_thresh,
                                     hysteresis_samples)  # do the step back in time stuff

            # append onsets
            corrected_onsets.append(onset_hold)

    # remove any zeros from corrected_onsets not at the start
    zero_loc = np.where(np.array(corrected_onsets) == 0)[0]
    if list(zero_loc):
        if zero_loc[0] == 0:
            zero_loc = zero_loc[1:]
            corrected_onsets = np.delete(corrected_onsets, zero_loc)

    # remove duplicates
    corrected_onsets = (list(set(corrected_onsets)))
    corrected_onsets.sort()

    # set defauls for analysis and initialise store lists
    time_decay = int(fs * 0.1)  # number of samples for 100 ms
    time_count = np.arange(0, time_decay, 1)
    decays = []
    all_spread = []
    all_decay = []
    all_attack = []
    all_centroid = []

    for count in range((len(corrected_onsets))):
        current_onset = corrected_onsets[count]
        if current_onset == corrected_onsets[-1]:
            attack_segment = attack_env[current_onset:]
            decay_segment = logenvelope[current_onset:]
            audio_segment = audio_samples[current_onset:]
        else:
            attack_segment = attack_env[current_onset:corrected_onsets[count + 1]]
            decay_segment = logenvelope[current_onset:corrected_onsets[count + 1]]
            audio_segment = audio_samples[current_onset:corrected_onsets[count + 1]]

        # only if the segment is long enough
        if len(attack_segment) > 1:
            ''' 
             Calculate attack time
            '''
            peak_idx = np.argmax(attack_segment)
            if peak_idx > 0:
                if not peak_idx == len(attack_segment)-1:
                    attack_subsegment = attack_segment[:peak_idx+1]
                else:
                    attack_subsegment = attack_segment

                attack_dyn_range = max(attack_subsegment) - min(attack_subsegment)

                # only perform the attack time analysis if dynamic range is greater than 10 %
                if attack_dyn_range >= 0.1:
                    min_idx = np.argmin(attack_subsegment)
                    attack_subsegment = attack_subsegment[min_idx:]
                    attack_subsegment = attack_subsegment - min(attack_subsegment)
                    attack_subsegment *= 1.0 / max(attack_subsegment)

                    start_idx = np.where(attack_subsegment > 0.1)[0][0]
                    end_idx = np.where(attack_subsegment > 0.9)[0][0]

                    # calculate attack time
                    if not end_idx == start_idx:
                        attack_time = (end_idx - start_idx) / float(fs)
                        all_attack.append(np.log10(attack_time))
                    else:
                        all_attack.append(np.log10(1.0/fs))

            '''
             Get spectral features
            '''
            centroid, spread = get_spectral_features(audio_segment, fs)
            all_spread.append(spread)
            all_centroid.append(centroid)

            '''
             Calculate normalised decay time
            '''
            peak_idx = np.argmax(decay_segment)
            decay_subsegment = decay_segment[peak_idx:]
            decay_dyn_range = max(decay_subsegment) - min(decay_subsegment)
            log_dyn_range = max(logenvelope) - min(logenvelope)

            # only perform decay analysis if decay is greater than 10% total dynamic range
            if decay_dyn_range >= (log_dyn_range * 0.1):
                # get the decay
                count2 = 0
                MSE = []
                r = []
                if count2 + time_decay <= len(decay_subsegment):
                    while count2 + time_decay < len(decay_subsegment):
                        # get the decay segment to be evaluated
                        eval_segment = decay_subsegment[count2:count2 + time_decay]

                        # get Pearson's r, test test for negative slope
                        r_hold = np.corrcoef(time_count, eval_segment)[1][0]
                        r.append(r_hold)

                        if r_hold < 0:
                            hold_poly_vals = np.polyfit(time_count, eval_segment, 1)
                            hold_best_fit = hold_poly_vals[1] + hold_poly_vals[0] * time_count
                            this_MSE = np.mean((eval_segment - hold_best_fit) ** 2)
                            MSE.append(this_MSE)
                        else:
                            MSE.append(float("inf"))

                        # step back 100 samples, code takes too long to run stepping back a single sample
                        count2 += 100

                    # identify most linear portion of decay
                    if np.min(MSE) != float("inf"):
                        best_idx = np.argmin(MSE)
                        eval_segment = decay_subsegment[best_idx * 100:(best_idx * 100) + time_decay]

                    # get linear predcition of decay
                    hold_poly_vals = np.polyfit(time_count, eval_segment, 1)
                    hold_best_fit = hold_poly_vals[1] + hold_poly_vals[0] * time_count

                    # calculate the decay time
                    decay = (hold_best_fit[-1] - hold_best_fit[0]) * 10
                    decays.append(decay)

                    # get normalised decay
                    normalised_dacay = decay / centroid
                    all_decay.append(normalised_dacay)

    '''
     Get mean values for parameters
    '''
    if all_centroid:
        mean_centroid = np.mean(all_centroid)
    else:
        mean_centroid = 0.0

    if all_decay:
        mean_decay = np.mean(all_decay)
    else:
        mean_decay = 0.0
    if all_attack:
        mean_attack = np.mean(all_attack)
        if mean_attack > 0:
            mean_attack = np.log10(mean_attack)
    else:
        mean_attack = 0.0
    if all_spread:
        mean_spread = np.mean(all_spread)
    else:
        mean_spread = 0.0

    ''' Implementation of a logistic regression model to calculate metallic probability '''
    # regression coefficients
    coefficients = [528.406157139, -0.237151712738, 0.000118653519107, 0.00170027628632, -1.08383115685]

    # apply linear coefficients
    attributes = [mean_decay, mean_attack, mean_spread, roughness, 1.0]
    logit_model = np.sum(np.array(coefficients) * np.array(attributes))

    # apply inverse of Logit function to obtain probability
    probability = np.exp(logit_model) / (1.0 + np.exp(logit_model))

    return probability
    # maintained for testing and development of model
    # return mean_decay, mean_attack, mean_spread, roughness, mean_centroid, probability
