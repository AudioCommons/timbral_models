import numpy as np
from scipy.signal import hilbert, butter, lfilter
import Timbral_Roughness as rough
import soundfile as sf
import librosa


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def return_loop(onset_loc, envelope, function_time_thresh, wiggle_threshold, hist_time_samples):
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
            evaluation_array = envelope[onset_loc-function_time_thresh-1:onset_loc-1]
        else:
            evaluation_array = envelope[:onset_loc-1]

        if min(evaluation_array) - current_sample < 0:
            min_idx = np.argmin(evaluation_array)

            new_onset_loc = min_idx + onset_loc - function_time_thresh - 1

            if new_onset_loc > 512:
                onset_loc = new_onset_loc
            else:
                return 0

        else:
            ''' Introduce the time and level hysteresis'''
            # check we can look back enough
            if (onset_loc - hist_time_samples - 1) > 0:
                hyst_evaluation_array = envelope[onset_loc-hist_time_samples-1:onset_loc - 1]
            else:
                hyst_evaluation_array = envelope[:onset_loc-1]

            # find the first value that's less than current sample
            all_match = np.where(hyst_evaluation_array < envelope[onset_loc])

            # if no minimum was found within the extended time
            if len(all_match[0]) == 0:
                return onset_loc

            last_min = all_match[0][-1]
            # last_idx = int(onset_loc-hist_time_samples-1 + last_min)
            last_idx = int(onset_loc-len(hyst_evaluation_array) + last_min)


            # get the dynamic range
            segment_dynamic_range = max(hyst_evaluation_array[last_min:]) - min(hyst_evaluation_array[last_min:])

            if segment_dynamic_range >= wiggle_threshold:
                # this is a seperate audio event
                return onset_loc
            else:
                onset_loc = last_idx
            # return onset_loc


def get_attack_envelope(time_samples, fs, decay_time=0.2, hold_time=0.1):
    """
     Get the attack envelope in the same way as the timbral hardness model
    
    :param hold_time: 
    :param decay_time: 
    :param fs: 
    :param time_samples: 
    :return: 
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


def timbral_metallic(fname):
    # use pysoundfile instead
    audio_samples, fs = sf.read(fname, always_2d=False)

    num_channels = np.shape(audio_samples)
    if len(num_channels) > 1:
        # take just the left channel
        audio_samples = audio_samples[:, 0]
    audio_samples *= 1.0 / max(abs(audio_samples))

    roughness = rough.timbral_roughness(fname)

    h = np.absolute(hilbert(audio_samples))

    # Filter requirements.
    order = 2
    cutoff = 50  # desired cutoff frequency of the filter, Hz

    envelope = butter_lowpass_filter(h, cutoff, fs, order)
    envelope *= 1.0 / max(envelope)

    # calculate the onsets using librosa
    onsets = librosa.onset.onset_detect(audio_samples, fs)
    onset_loc = onsets * 512

    # set defaults for backwards stepping
    corrected_onsets = []
    previous_onset_idx = -1
    time_thresh = int(0.01 * fs)  # 10 ms time threshold for backwards analysis of rising
    wiggle_percent = 5  # percentage of dynamic range to allow variation
    wiggle_thresh = (max(envelope) - min(envelope)) * wiggle_percent
    hysteresis_samples = int(0.2 * fs)

    for onset_idx in onsets:
        if onset_idx <= 1:
            # can't correct for this
            corrected_onsets.append(0)
        else:
            onset_loc = onset_idx * 512  # actual onset location in samples
            onset_hold = return_loop(onset_loc, envelope, time_thresh, wiggle_thresh,
                                     hysteresis_samples)  # do the step back in time stuff

            corrected_onsets.append(onset_hold)
            if onset_hold > 0:
                previous_onset_idx = onset_hold - 1  # -1 is to allow the code to return to previous values in the loop

    # remove any zeros from corrected_onsets
    zero_loc = np.where(np.array(corrected_onsets) == 0)[0]
    if list(zero_loc):
        if zero_loc[0] == 0:
            zero_loc = zero_loc[1:]
    corrected_onsets = np.delete(corrected_onsets, zero_loc)

    # remove duplicates
    prev_onset = -1
    onset_hold = []

    for onset in corrected_onsets:
        if onset > prev_onset:
            onset_hold.append(onset)
            prev_onset = onset
    corrected_onsets = onset_hold

    logenvelope = np.log10(envelope)
    time_decay = int(fs * 0.1)  # number of samples for 100 ms
    time_count = np.arange(0, time_decay, 1)
    decays = []
    all_spread = []
    all_alpha = []
    all_attack = []
    all_centroid = []

    for count in range(len(corrected_onsets)):
        current_onset = corrected_onsets[count]
        if current_onset == corrected_onsets[-1]:
            eval_segment = logenvelope[current_onset:]
            time_segment = audio_samples[current_onset:]
        else:
            eval_segment = logenvelope[current_onset:corrected_onsets[count+1]]
            time_segment = audio_samples[current_onset:corrected_onsets[count+1]]

        # get the centroid and spread
        window = np.hanning(len(time_segment))
        spectrum = np.fft.fft(window * time_segment)
        spectrum = np.absolute(spectrum[0:int(len(spectrum)/2)+1])
        freq = np.arange(0, len(spectrum), 1) * fs / (2.0 * (len(spectrum) - 1))
        centroid = sum(spectrum*freq) / float(sum(spectrum))
        all_centroid.append(centroid)
        spectral_spread = np.sqrt(sum(((freq - centroid) ** 2) * spectrum) / sum(spectrum))
        all_spread.append(spectral_spread)

        # get the decay
        count2 = 0
        MSE = []
        r = []
        if count2+time_decay <= len(eval_segment):
            while count2+time_decay < len(eval_segment):
                eval_subsegment = eval_segment[count2:count2+time_decay]

                r_hold = np.corrcoef(time_count, eval_subsegment)[1][0]
                r.append(r_hold)

                if r_hold < 0:
                    hold_poly_vals = np.polyfit(time_count, eval_subsegment, 1)
                    hold_best_fit = hold_poly_vals[1] + hold_poly_vals[0] * time_count
                    this_MSE = np.mean((eval_subsegment - hold_best_fit) ** 2)
                    MSE.append(this_MSE)
                else:
                    MSE.append(float("inf"))

                count2 += 100

            if np.min(MSE) != float("inf"):
                best_idx = np.argmin(MSE)
                eval_subsegment = eval_segment[best_idx*100:(best_idx*100)+time_decay]

                hold_poly_vals = np.polyfit(time_count, eval_subsegment, 1)
                hold_best_fit = hold_poly_vals[1] + hold_poly_vals[0] * time_count

                decay = (hold_best_fit[-1] - hold_best_fit[0]) * 10
                decays.append(decay)

                # calculate normalised decay
                alpha = decay / float(centroid)
                all_alpha.append(alpha)






















    ''' Get attack time'''
    # get the attack envelope
    decay_time = 0.2  # in seconds, 200 ms
    hold_time = 0.01  # time to hold samples before decay, in seconds
    attack_env = get_attack_envelope(audio_samples, fs, decay_time, hold_time)

    # prepend 512 zeros to allow for initial onset to be detected
    attack_env = np.lib.pad(attack_env, (512, 0), 'constant', constant_values=(0, 0))

    # Calculate onset idx using librosa
    onsets = librosa.onset.onset_detect(audio_samples, fs)

    if np.size(onsets):

        # identify actual attack start times
        time_thresh = int(0.01 * fs)  # 10 ms time look-back time, in samples
        hysteresis_samples = int(decay_time * fs)  # hysteresis time, in samples
        hist_percent = 5  # percentage of dynamic range to allow variation
        hist_thresh = (max(attack_env) - min(attack_env)) * hist_percent

        # empty array for storing exact onset idxs
        corrected_onsets = []

        for onset_idx in onsets:
            # if the onset is 1 or 0, it's too close to the start to be corrected (1 is here due to zero padding)
            if onset_idx <= 1:
                corrected_onsets.append(0)
            else:
                # actual onset location in samples (librosa uses 512 window size by default)
                onset_loc = onset_idx * 512
                # use return loop to identify actual onset idx
                onset_hold = return_loop(onset_loc, attack_env, time_thresh, hist_thresh, hysteresis_samples)

                # store new onset in list
                corrected_onsets.append(onset_hold)

        # remove duplicates
        corrected_onsets = list(set(corrected_onsets))

        # zero is returned from return_loop if no valid onset identified
        # remove zeros (except the first)
        zero_loc = np.where(np.array(corrected_onsets) == 0)[0]
        # ignore if the first value is zero
        if list(zero_loc):
            if zero_loc[0] == 0:
                zero_loc = zero_loc[1:]
        corrected_onsets = np.delete(corrected_onsets, zero_loc)

        '''
         Compare onset segments against the dynamic range to remove erroneous onsets in noise.  
         If the onset segment (samples between adjacent onsets) has a dynamic range less than 
         10% of total dynamic range, remove this onset.
        '''
        last_value = corrected_onsets[-1]
        envelope_dyn_range = max(attack_env) - min(attack_env)
        threshold = 0.1 * envelope_dyn_range

        for i in range(len(corrected_onsets)):
            if corrected_onsets[i] == last_value:
                segment = attack_env[corrected_onsets[i]:]
            else:
                segment = attack_env[corrected_onsets[i]:corrected_onsets[i+1]]

            # only conduct if the segment if greater than 1 sample long
            if len(segment) > 1:
                # find attack portion SNR
                peak_idx = np.argmax(segment)
                if peak_idx > 0:
                    # get the dynamic range of the attack portion
                    attack_segment = segment[:peak_idx]
                    seg_dyn_range = max(attack_segment) - min(attack_segment)
                    if seg_dyn_range >= threshold:
                        # thd_corrected_onsets.append(corrected_onsets[i])

                        min_idx = np.argmin(attack_segment)
                        attack_segment = attack_segment[min_idx:]
                        attack_segment = attack_segment - min(attack_segment)
                        attack_segment *= 1.0 / max(attack_segment)

                        start_idx = np.where(attack_segment > 0.1)[0][0]
                        end_idx = np.where(attack_segment > 0.9)[0][0]

                        if not end_idx == start_idx:
                            attack_time = (end_idx - start_idx) / float(fs)
                            all_attack.append(attack_time)
                        else:
                            all_attack.append(0)



    mean_centroid = np.mean(centroid)

    if all_alpha:
        mean_alpha = np.mean(all_alpha)
    else:
        mean_alpha = 0
    if all_attack:
        mean_attack = np.mean(all_attack)
        if mean_attack > 0:
            mean_attack = np.log10(mean_attack)
    else:
        mean_attack = 0
    if all_spread:
        mean_spread = np.mean(all_spread)
    else:
        mean_spread = 0


    ''' get the metallic probability from logistic regression model '''

    coefficients = [130.633579243, -0.409398377632, 0.000211216295667, 0.001432156288, -1.70566517363]
    # coefficients if using spectral centroid
    # coefficients = [129.539857931, -0.405972586334, 0.000197263145611, 0.00118344866937, 1.6049278962e-05,-1.70277588538]

    attributes = [mean_alpha, mean_attack, mean_spread, roughness, 1.0]
    # attributes = [mean_alpha, mean_attack, mean_spread, roughness, mean_centroid, 1.0]
    logit_model = np.sum(np.array(coefficients) * np.array(attributes))
    probability = np.exp(logit_model) / (1.0 + np.exp(logit_model))

    return mean_alpha, mean_attack, mean_spread, roughness, mean_centroid, probability
    # return probability
