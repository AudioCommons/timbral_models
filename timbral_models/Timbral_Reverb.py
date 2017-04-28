import numpy as np
from scipy.signal import spectrogram
import soundfile as sf


def strictly_decreasing(eval_list):
    """
     Function which assess if an array is constantly decreasing
    """
    return all(x > y for x, y in zip(eval_list, eval_list[1:]))


def find_most_linear(full_eval_array, find_dyn_range, time_step):
    """
    This function identifies the most linear section which has the dynamic range of at least find_dyn_range
    """
    # test array is long enough
    if len(full_eval_array) < 3:
        return 0

    st = 0
    en = st + 1
    find_st_idx = []
    find_en_idx = []
    MSE = []
    while st < (len(full_eval_array) - 1):
        hold_array = full_eval_array[st:en]
        # check array dyn range
        if (max(hold_array) - min(hold_array)) >= find_dyn_range:
            find_st_idx.append(st)
            find_en_idx.append(en)
            hold_time_array = np.arange(0, len(hold_array), dtype=np.float)
            hold_time_array *= time_step
            hold_poly_vals = np.polyfit(hold_time_array, hold_array, 1)
            hold_best_fit = hold_poly_vals[1] + hold_poly_vals[0] * hold_time_array
            MSE.append(np.mean((hold_array - hold_best_fit) ** 2))
        if en >= (len(full_eval_array) - 1):
            st += 1
            en = st + 1
        else:
            en += 1
    if MSE:
        choose_idx = np.argmin(MSE)
        choose_array = full_eval_array[find_st_idx[choose_idx]:find_en_idx[choose_idx]]
        return choose_array
    else:
        return 0


def timbral_reverb(fname):
    """
     This function estimates the reverberation time (RT60) of an audio sample.
     A value of 0 is returned if the RT60 could not be estimated by the algorithm.

     The algorithm used is an adaptation of Prego et al. [2015].
     This algorithm identifies free decay regions (FDR), areas where the energy decays in adjacent time frames, for
     each frequency subband within the spectrogram of the signal.  A decay region is then identified by finding the
     best linear fit within this decay, and extrapolating the gradient to estimate the RT60.

     Prego, T., Lima, A., Zambrano-Lopez, R., and Netto, S. 'Blind estimators for reverberation time and
     direct-to-reverberant energy ratio using subband speech decomposition.' IEEE Workshops on Applications of signal
     processing to audio and acoustics, October, 2015.

    :param fname:   File name of the audio
    :return:        Estimated RT60 of the signal, returns 0 if this cannot be calculated.
    """
    # use pysoundfile to read in audio
    audio_samples, fs = sf.read(fname, always_2d=False)

    num_channels = np.shape(audio_samples)
    if len(num_channels) > 1:
        # take just the left channel
        audio_samples = audio_samples[:, 0]

    '''
      Get a spectrogram of the signal
    '''
    # normalise audio
    audio_samples *= (1.0 / max(audio_samples))
    nfft = 2048
    hop_size = nfft / 4
    freq, time, spec = spectrogram(audio_samples, fs, 'hamming', nfft,
                                   hop_size, nfft, 'constant', True, 'spectrum')

    '''
      Set limits of spectrogram
      The original algorithm suggests only looking at up to 4kHz
    '''
    min_freq = 200
    max_freq = 4000
    min_idx = np.argmax(freq > min_freq)
    max_idx = np.argmax(freq > max_freq)

    '''
      Do the subband detection
    '''
    Llim = int((0.5 * fs) / (nfft - hop_size))
    all_RT = []

    # for each subband in the spectrogram
    for freq_bin in range(min_idx, max_idx):
        this_freq_RT = []
        true_values = []
        L_values = []
        Llim_var = Llim

        # identify the free decay ranges (FDR)
        while Llim_var > 3 and not true_values:
            for i in range(0, (len(spec[0, :])-Llim_var)):
                evaluation_segment = spec[freq_bin, i:i+Llim_var]

                if strictly_decreasing(evaluation_segment):
                    true_values.append(i)
                    L_values.append(Llim_var)

            # if we found it on the first pass
            if Llim_var == Llim and true_values:
                # for each appropriate decrease
                for i in range(len(true_values)):
                    start_idx = true_values[i]
                    # can we increase the decay length?
                    decrease = True
                    new_Llim = Llim_var + 1
                    while decrease and (start_idx+Llim_var) <= len(spec[0, :]):
                        new_evaluation_segment = spec[freq_bin, start_idx:start_idx+new_Llim]
                        if strictly_decreasing(new_evaluation_segment):
                            L_values[i] += 1
                            new_Llim += 1
                        else:
                            decrease = False

            Llim_var -= 1

        # check for overlapping free decay regions
        overlap = []
        for i in range(1, len(true_values)):
            if true_values[i] < (L_values[i-1] + true_values[i-1]):
                overlap.append(i)

        # remove the overlaps
        if overlap:
            overlap.reverse()
            for i in overlap:
                del true_values[i]
                del L_values[i]

        # for each appropriate free decay range (FDR)
        time_step = (time[1] - time[0])
        for i in range(len(true_values)):
            evaluation_segment = spec[freq_bin, true_values[i]:true_values[i]+L_values[i]]
            evaluation_segment += 1.0 / 2**16  # add one bit depth (at 16 bit) to avoid divide by zero

            # get schroeder integral
            tpower = sum(evaluation_segment)
            sum_power = np.zeros(len(evaluation_segment))

            for j in range(len(evaluation_segment)):
                sum_power[j] = sum(evaluation_segment[j:])

            # schroeder inegration of the SFDR
            c = 10 * np.log10(sum_power / float(tpower))

            # check sufficient dynamic range
            if (c[0] - c[-1]) > 15:

                # find -5 dB point
                start_idx = np.argmax(c <= -5.0)

                # is the first decay greater than 10 dB
                if (c[start_idx] - c[start_idx+1] > 10) and (start_idx <= len(c)-2):
                    # estimate the RT60
                    decay_range = c[start_idx] - c[start_idx+1]
                    decay_time = time_step

                    RT60 = (60.0 / decay_range) * decay_time * 1000
                    this_freq_RT.append(RT60)

                # do not do this if the array is too short
                elif start_idx <= len(c)-3:
                    # calculate the MSE for each pair
                    MSE = []
                    test_arrays = []

                    for count in range(3, len(c[start_idx:])+1):
                        test_array = c[start_idx:start_idx+count]
                        test_arrays.append(test_array)
                        time_array = np.arange(0, count, dtype=np.float)
                        time_array *= time_step
                        poly_vals = np.polyfit(time_array, test_array, 1)
                        best_fit = poly_vals[1] + poly_vals[0] * time_array

                        MSE.append(np.mean((test_array - best_fit) ** 2))

                    # identify best fitting model
                    min_rms_idx = np.argmin(MSE)
                    test_array = test_arrays[min_rms_idx]
                    dyn_range = test_array[0] - test_array[-1]

                    if dyn_range <= -10:
                        # print 'best fit is greater than 10 dB'
                        test_array = c[start_idx:start_idx + min_rms_idx + 3]
                        time_array = np.arange(0, len(test_array), dtype=np.float)
                        time_array *= time_step

                        poly_vals = np.polyfit(time_array, test_array, 1)
                        best_fit = poly_vals[1] + poly_vals[0] * time_array

                        decay_range = best_fit[0] - best_fit[-1]
                        decay_time = len(best_fit) * time_step

                        RT60 = (60.0 / decay_range) * decay_time * 1000
                        this_freq_RT.append(RT60)

                    else:
                        # get the next best dynamic range
                        # what dynamic range is available?
                        this_dyn_range = c[start_idx] - c[-1]

                        if this_dyn_range >= 60:
                            # find the most linear 60dB
                            selected_array = find_most_linear(c[start_idx:], 60, time_step)
                        elif this_dyn_range >= 40:
                            # find the most linear 40 dB
                            selected_array = find_most_linear(c[start_idx:], 40, time_step)
                        elif this_dyn_range >= 20:
                            # find the most linear 20 dB
                            selected_array = find_most_linear(c[start_idx:], 20, time_step)
                        elif this_dyn_range >= 10:
                            # find the most linear 10 dB
                            selected_array = find_most_linear(c[start_idx:], 10, time_step)
                        else:
                            # there is no suitable output
                            selected_array = 0

                        if isinstance(selected_array, int):
                            RT60 = 0
                        else:
                            time_array = np.arange(0, len(selected_array), dtype=np.float)
                            time_array *= time_step

                            poly_vals = np.polyfit(time_array, selected_array, 1)
                            best_fit = poly_vals[1] + poly_vals[0] * time_array

                            decay_range = best_fit[0] - best_fit[-1]
                            decay_time = len(best_fit) * time_step

                            RT60 = (60.0 / decay_range) * decay_time * 1000
                            this_freq_RT.append(RT60)

        all_RT.append(this_freq_RT)

    '''
      Convert between subband RT and actual RT values
    '''
    values = []
    for RTs in all_RT:
        # is it a value?
        if RTs:
            values.append(np.median(RTs))

    if values:
        medianRT60 = np.median(values)

        # divide by 3 seems to give a reasonable estimate from the training dataset
        outputRt60 = medianRT60 / 3.0

        # the Prego et al. [2015] algorithm, does not produce reasonable results for the training dataset:
        # outputRt60 = 1000.0 * (3.4 * (medianRT60 / 1000.0) - 1.170)

    else:
        # unable to predict RT60
        medianRT60 = 0
        outputRt60 = 0

    return outputRt60




