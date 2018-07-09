import numpy as np
import librosa
import soundfile as sf
from scipy.signal import spectrogram


def db2mag(dB):
    mag = 10 ** (dB / 20.0)
    return mag


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


def calculate_attack_time(envelope_samples, fs, calculate_attack_segment=True, thresh_no=8, normalise=True, m=3,
                          calculation_type='min_effort', gradient_calulation_type='all', return_descriptive_data=False,
                          max_attack_time = -1):
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
                print('unable to calculate attack gradient with the \'mean\' method, reverting to \'all\' method.')

        if gradient_calulation_type == 'mean':
            # calculate the gradient based on the weighted mean of each attack
            threshold_step = dyn_range / (thresh_no + 2)

            gradient_thresh_array = np.arange(start_level, end_level+(threshold_step*dyn_range),
                                              (threshold_step*dyn_range))
            cross_threshold_times = np.zeros(len(gradient_thresh_array))
            cross_threshold_values = np.zeros(len(gradient_thresh_array))
            gradient_envelope_segment = envelope_samples[th_start_idx:th_end_idx+1]

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
        thresholds_to_return = [calculation_type, th_start_idx+min_pre_peak_idx, th_end_idx+min_pre_peak_idx,
                                threshold_idxs+min_pre_peak_idx]

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
        thresholds_to_return = [calculation_type, th_start_idx+min_pre_peak_idx, th_end_idx+min_pre_peak_idx]

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
    hold_env = envelope_samples[int(th_start_idx):int(th_end_idx)+1]
    t = np.arange(0, len(hold_env)) / float(fs)
    temp_centroid = np.sum(t * hold_env) / np.sum(hold_env)
    temp_centroid /= float(len(hold_env))

    if return_descriptive_data:
        return attack_time, attack_gradient, int(th_start_idx+min_pre_peak_idx), temp_centroid, thresholds_to_return
    else:
        return attack_time, attack_gradient, int(th_start_idx+min_pre_peak_idx), temp_centroid


def calculate_onsets(audio_samples, envelope_samples, fs, look_back_time=20, hysteresis_time=300, hysteresis_percent=10,
                     onset_in_noise_threshold=10, threshold_correction='onset_strength',
                     minimum_onset_time_separation=100, nperseg=512):
    '''
     Calculate onset idx using librosa

     I'm still adding new features to this so remember to update the definitions of all the input features
    '''
    onsets = librosa.onset.onset_detect(audio_samples, fs, backtrack=True, units='samples')

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
                onset_loc = np.array(onset_idx).astype('int')

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
                for i in range(len(corrected_onsets)-1):
                    # are the last two onsets too close?
                    if abs(corrected_onsets[i+1] - corrected_onsets[i]) < minimum_onset_time_separation_samples:
                        onsets_to_remove.append(i+1)

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

        strength_onset_times = np.array(corrected_onsets) // 512
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
                strength_onset_times = np.delete(strength_onset_times,onset_idx)
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
        seg = spec[:,time_count]
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
                rollon.append(f[rollon_counter-1])
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
                    bandwidth.append(f[rolloff_idx] - f[rollon_counter-1])
            else:
                bandwidth.append(0)
        else:
            bandwidth.append(0)

        if tpower > 0.05:
            centroid.append(np.sum(seg*f) / np.sum(seg))
            centroid_power.append(tpower)
    if return_centroid:
        return bandwidth, t, f, np.average(centroid, weights=centroid_power)
    else:
        return bandwidth, t, f


def timbral_hardness(fname, dev_output=False, max_attack_time=0.3, bandwidth_thresh_db=-50, phase_correction=False):
    """
     This function calculates the apparent hardness of an audio file.
     Currently the function only works on the left channel of multi-channel audio files.

      :param fname: Audio filename to be analysed, including full file path and extension.
      :return: returns the hardness as a float representing the relative hardness.
    """
    import scipy.stats
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
                audio_samples = audio_samples[:,0]  #[:,1] *= -1.0
            else:
                audio_samples = np.sum(audio_samples, axis=1)
        else:
            audio_samples = np.sum(audio_samples, axis=1)

    # normalise audio
    audio_samples /= max(abs(audio_samples))

    # zero pad the signal
    audio_samples = np.lib.pad(audio_samples, (nperseg+1, 0), 'constant', constant_values=(0.0, 0.0))

    # calculate the envelope of the signal
    envelope = sample_and_hold_envelope_calculation(audio_samples, fs, decay_time=0.1)
    envelope_time = np.arange(len(envelope)) / float(fs)

    # calculate the onsets
    original_onsets = calculate_onsets(audio_samples, envelope, fs, nperseg=nperseg)
    onset_strength = librosa.onset.onset_strength(audio_samples, fs)

    '''
      Calculate the spectrogram so that the bandwidth can be created
    '''
    mag = db2mag(bandwidth_thresh_db)
    bandwidth, t, f = get_bandwidth_array(audio_samples, fs, nperseg=nperseg, overlap_step=32, rolloff_thresh=mag)
    all_bandwidth_max = []
    all_attack_time = []
    all_max_strength = []
    all_max_strength_bandwidth = []

    '''
      Get bandwidth onset times and max bandwidth
    '''
    # If onsets don't exist, set it to time zero
    if not original_onsets:
        original_onsets = [0]

    if original_onsets:
        onsets = np.array(original_onsets) - nperseg
        onsets[onsets < 0] = 0
        bandwidth_onset = np.array(onsets) // 32 # overlap_step=32

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

            bandwidth_fs = fs / float(32.0) #overlap_step=32
            if max(bandwidth_seg) > 0:
                bandwidth_attack = calculate_attack_time(bandwidth_seg, bandwidth_fs, calculation_type='fixed_threshold',
                                                         max_attack_time=max_attack_time)
            else:
                bandwidth_attack = []

            # attack_time, attack_gradient, th_start_idx, temp_centroid
            if bandwidth_attack:
                start_idx = bandwidth_attack[2]
                if max_attack_time > 0:
                    max_attack_time_samples = int(max_attack_time * bandwidth_fs)
                    if len(hold_bandwidth_seg[start_idx:]) > start_idx+max_attack_time_samples:
                        all_bandwidth_max.append(max(hold_bandwidth_seg[start_idx:start_idx+max_attack_time_samples]))
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
                strength_seg = np.array(onset_strength[(onset // 512):])
            else:
                attack_seg = np.array(envelope[onset:original_onsets[onset_count + 1]])
                strength_seg = np.array(onset_strength[(onset // 512):(original_onsets[onset_count+1] // 512)])

            attack_time = calculate_attack_time(attack_seg, fs, max_attack_time=max_attack_time)
            all_attack_time.append(attack_time[0])

            all_max_strength.append(max(strength_seg))
            if bandwidth_attack:
                all_max_strength_bandwidth.append(max(strength_seg))



        mean_bandwidth_max = np.mean(all_bandwidth_max)
        mean_attack_time = np.mean(all_attack_time)

        # weight all metrics prior to taking the mean
        mean_weighted_bandwidth_max = np.average(all_bandwidth_max, weights=all_max_strength_bandwidth)
        mean_weighted_attack_time = np.average(all_attack_time, weights=all_max_strength)

    else:

        mean_bandwidth_max = 0
        mean_attack_time = 0

        # weight all metrics prior to taking the mean
        mean_weighted_bandwidth_max = 0
        mean_weighted_attack_time = 0


    # return the values
    if dev_output:
        return mean_bandwidth_max, mean_attack_time, mean_weighted_bandwidth_max, mean_weighted_attack_time
    else:
        '''
         Apply regression model
        '''

        all_metrics = np.ones(3)
        all_metrics[0] = mean_weighted_bandwidth_max
        all_metrics[1] = mean_attack_time

        coefficients = np.array([0.00200588674443, -4.98150425805, 27.5863066923])

        hardness = np.sum(all_metrics * coefficients)
        return hardness


