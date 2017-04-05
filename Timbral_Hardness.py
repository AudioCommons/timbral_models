import numpy as np
import librosa
import soundfile as sf


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
            evaluation_array = envelope[onset_loc-function_time_thresh-1:onset_loc-1]
        else:
            evaluation_array = envelope[:onset_loc-1]

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
                hyst_evaluation_array = envelope[onset_loc-hist_time_samples-1:onset_loc - 1]
            else:
                hyst_evaluation_array = envelope[:onset_loc-1]

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
                onset_loc = last_idx


def timbral_hardness(fname):
    """
     This function calculates the apparent hardness of an audio file.
     Currently the function only works on the left channel of multi-channel audio files.

      :param fname: Audio filename to be analysed, including full file path and extension.
      :return: returns the hardness as a float representing the relative hardness.
    """
    # use pysoundfile to read audio
    audio_samples, fs = sf.read(fname, always_2d=False)

    # get left channel only
    num_channels = np.shape(audio_samples)
    if len(num_channels) > 1:
        audio_samples = audio_samples[:, 0]

    '''
     Calculate the envelope of the audio with a 'sample and hold' style function.
     This ensures that the minimum attack time is not limited by low-pass filtering,
     a common method of obtaining the envelope.
    '''
    # rectify the audio signal
    abs_samples = abs(audio_samples)
    envelope = []

    # set parameters for envelope function
    decay_time = 0.2  # in seconds, 200 ms
    decay = max(abs_samples) / (decay_time * fs)  # decay rate relative to peak level of audio signal
    hold_time = 0.01  # time to hold samples before decay, in seconds
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
    envelope = np.array(envelope)

    # prepend 512 zeros to allow for initial onset to be detected
    envelope = np.lib.pad(envelope, (512, 0), 'constant', constant_values=(0, 0))
    audio_samples = np.lib.pad(audio_samples, (512, 0), 'constant', constant_values=(0, 0))

    '''
     Calculate onset idx using librosa
    '''
    onsets = librosa.onset.onset_detect(audio_samples, fs)

    '''
     Preallocate arrays for storing variables
    '''
    all_attack_time = []
    all_attack_gradient = []
    all_attack_centroid = []

    '''
     Identify actual onset times using return_loop
    '''
    # set values for return_loop method
    time_thresh = int(0.01 * fs)  # 10 ms time look-back time, in samples
    hysteresis_samples = int(decay_time * fs)  # hysteresis time, in samples
    hist_percent = 5  # percentage of dynamic range to allow variation
    hist_thresh = (max(envelope) - min(envelope)) * hist_percent

    # only conduct analysis if there are onsets detected
    if np.size(onsets):
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
                onset_hold = return_loop(onset_loc, envelope, time_thresh, hist_thresh, hysteresis_samples)

                # store new onset in list
                corrected_onsets.append(onset_hold)

        # zero is returned from return_loop if no valid onset identified
        # remove zeros (except the first)
        zero_loc = np.where(np.array(corrected_onsets) == 0)[0]
        # ignore if the first value is zero
        if list(zero_loc):
            if zero_loc[0] == 0:
                zero_loc = zero_loc[1:]
        corrected_onsets = np.delete(corrected_onsets, zero_loc)

        '''
         Remove repeated onsets and compare onset segments against the dynamic range
         to remove erroneous onsets in noise.  If the onset segment (samples between
         adjacent onsets) has a dynamic range less than 10% of total dynamic range,
         remove this onset.
        '''
        thd_corrected_onsets = []
        last_value = corrected_onsets[-1]
        envelope_dyn_range = max(envelope) - min(envelope)
        threshold = 0.1 * envelope_dyn_range

        for i in range(len(corrected_onsets)):
            if corrected_onsets[i] == last_value:
                next_onset_idx = len(envelope) - 1  # index of the last envelope value
            else:
                next_onset_idx = corrected_onsets[i + 1]

            # check for duplicates looking forwards
            if not corrected_onsets[i] == next_onset_idx:
                # get the onset segment
                segment = envelope[corrected_onsets[i]:next_onset_idx]

                # only conduct if the segment if greater than 1 sample long
                if len(segment) > 1:
                    # find attack portion SNR
                    peak_idx = np.argmax(segment)
                    if peak_idx > 0:
                        # get the dynamic range of the attack portion
                        seg_dyn_range = max(segment) - min(segment[:peak_idx])
                        if seg_dyn_range >= threshold:
                            thd_corrected_onsets.append(corrected_onsets[i])

        '''
         Calculate the attack time using the minimum effort method.
        '''
        thresh_no = 8  # number of thresholds to calculate
        centroid_int_time = 0.125  # 125ms maximum integration time for spectral centroid
        centroid_int_samples = int(centroid_int_time * fs)  # number of samples for attack time integration

        for i in range(len(thd_corrected_onsets)):
            # get the segment containing the attack, both the envelope and audio
            start_idx = (thd_corrected_onsets[i])
            # if it's the last element
            if thd_corrected_onsets[i] == thd_corrected_onsets[-1]:
                segment = envelope[start_idx:]
                audio_segment = audio_samples[start_idx:]
            else:
                end_idx = thd_corrected_onsets[i + 1]
                segment = envelope[start_idx:end_idx]
                audio_segment = audio_samples[start_idx:end_idx]

            # normalise the segments
            segment /= float(max(segment))
            audio_segment /= float(max(abs(audio_segment)))

            # identify pre-attack segment
            peak_idx = np.argmax(segment)
            min_pre_peak_idx = np.argmin(segment[:peak_idx])

            attack_segment = segment[min_pre_peak_idx:peak_idx]

            # get threshold time array (we ignore 0% and 100%)
            threshold_step = 1 / float(thresh_no + 2)
            seg_dyn_range = max(attack_segment) - min(attack_segment)
            thresh_level = np.linspace(threshold_step, (1 - threshold_step), thresh_no+1)
            thresh_level = (thresh_level * seg_dyn_range) + min(attack_segment)
            threshold_idxs = np.zeros(thresh_no+1)

            # get indexes for when threshold is crossed
            for j in range(len(thresh_level)):
                threshold_hold = np.argmax(attack_segment >= thresh_level[j])
                threshold_idxs[j] = threshold_hold + min_pre_peak_idx

            # calculate effort values (distances between thresholds)
            effort = np.zeros(len(threshold_idxs) - 1)
            for j in range(len(effort)):
                effort[j] = threshold_idxs[j + 1] - threshold_idxs[j]

            # get the mean effort value
            effort_mean = np.mean(effort)
            m = 3  # value used for computation, set as 3 in CUIDADO project
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

            '''
             Will only output when the step is greater than 1 sample
            '''
            if th_start_idx == th_end_idx:
                attack_time = 1.0 / float(fs)
            else:
                attack_time = (th_end_idx - th_start_idx + 1) / float(fs)

            attack_time = np.log10(attack_time)
            all_attack_time.append(attack_time)

            '''
             Calculate the attack gradient
            '''
            # find the level of the first sample used
            start_level = segment[int(th_start_idx)]
            # find the level of the last sample used
            end_level = segment[int(th_end_idx)]
            if start_level == end_level:
                if th_start_idx > 0:
                    start_level = segment[int(th_start_idx)-1]
                else:
                    # set start level to zero if onset is at the first sample (indicating a step function at time zero)
                    start_level = 0

            level_change = end_level - start_level

            gradient = level_change / attack_time
            all_attack_gradient.append(gradient)

            '''
             Calculate the spectral centroid of the attack
            '''
            if (centroid_int_samples + th_start_idx) <= len(segment):
                post_attack_segment = audio_segment[int(th_start_idx):int(th_start_idx + centroid_int_samples)]
            else:
                post_attack_segment = audio_samples[int(th_start_idx):len(segment)]

            attack_spectrum = np.fft.fft(post_attack_segment)
            attack_spectrum = abs(attack_spectrum[:int(len(attack_spectrum) / 2)])  # take the real half od the spectrum
            fls = np.linspace(0, fs/2, len(attack_spectrum))

            # calculate the centroid
            centroid = np.log10(np.sum(attack_spectrum * fls) / np.sum(attack_spectrum))
            all_attack_centroid.append(centroid)

    '''
     Get mean values for all onsets
    '''
    mean_attack_time = np.mean(all_attack_time)
    mean_attack_gradient = np.mean(all_attack_gradient)
    mean_attack_centroid = np.mean(all_attack_centroid)

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

    coefficients = np.array([-367.825120662, -3337.86806786, 432.934655169, -628.173834828, 105.07915434,
                             1008.30167456, 191.819736608, -1427.29046921])

    hardness = np.sum(all_metrics * coefficients)

    return hardness
    # this is for testing the function
    # return mean_attack_time, mean_attack_gradient, mean_attack_centroid, hardness
