import numpy as np
import soundfile as sf


def calcFrequencyScale(spectrum_length, fs, blockSize):
    """ calculate frequency scale """
    freq = np.array(range(0, spectrum_length)) * fs / float(blockSize)
    return freq


def calcLowFrequencyLimit(fls, noct, max_idx):
    """ Calculates an array containing the indexes of lower frequency
    bounds for the spectral smoothing. """
    # floats required due to integer division in Python 2.7
    f_lower = fls[0:max_idx] / (2.0 ** (1 / (2.0 * noct)))
    step_size = fls[1] - fls[0]
    approx_idx = f_lower / (1.0 * step_size)
    f_lower = np.round(approx_idx).astype(int)
    return f_lower


def calcUpperFrequencyLimit(fls, noct, max_idx):
    """ calculate an array containing the indexes of upper frequency
    bounds for the spectral smoothing. """
    # floats required due to integer division in Python 2.7
    f_upper = fls[0:max_idx] * (2.0 ** (1.0 / (2.0 * noct)))
    step_size = fls[1] - fls[0]
    approx_idx = f_upper / float(step_size)
    f_upper = np.round(approx_idx).astype(int)
    return f_upper


def calcMaxIDX(fls, noct):
    """ calculates the max index of the frequency scale with respect
    to the level of smoothing applied.  This ensures that only
    frequencies below the Nyquist limit are used for smoothing """
    freq_l = fls[-1] / (2.0 ** (1 / (2.0 * noct)))
    max_idx = np.array(abs(fls - freq_l)).argmin()
    return max_idx


def calcMinIDX(fls, minFreq):
    """ calculates the index of the minimum frequency used for
    smoothing.  This ensures that frequencies below the lower
    limit are not considered.  """
    min_idx = np.argmax(fls >= minFreq)
    return min_idx


def spectralSmoothing(spectrum, f_upper, f_lower):
    """ Performs the spectral smoothing on a step-by-step basis.
    For each valid frequency bin, the magnitude is averaged over
    1/Noct below and above the centre frequency.
    """
    smoothed_array = np.zeros(len(f_upper))
    for i in range(len(smoothed_array)):
        # if statement is required since clever indexing isn't that clever.
        if f_upper[i] == f_lower[i]:
            smoothed_array[i] = spectrum[f_lower[i]]
        else:
            smooth_values = spectrum[f_lower[i]:(f_upper[i] + 1)]
            smoothed_array[i] = np.mean(smooth_values)
    return smoothed_array


def calcCrossoverIDX(fls, crossover_freq):
    """ Calculates the index of the crossover frequency """
    crossover_idx = np.argmax(fls >= crossover_freq)
    return crossover_idx


def timbral_brightness(fname):
    """ The main process block where computations are completed with VAMP.
     Data is accepted in the time domain, where a Hamming window is applied,
     and the FFT of this data taken.  The magnitude spectrum is then smoothed
     according to Noct (set to half-octave band smoothing by default),
     before calculating the spectral centroid above the crossover frequency,
     and the ratio of energy above the crossover frequency to all energy.
    """
    # use pysoundfile instead
    audio_samples, fs = sf.read(fname, always_2d=False)

    num_channels = np.shape(audio_samples)
    if len(num_channels) > 1:
        # take just the left channel
        audio_samples = audio_samples[:, 0]

    # initialise default settings
    stepSize = 1024
    blockSize = 2048
    threshold = 0
    n_oct = 2
    centroid_list = []
    crossover = 3000
    crossover_idx = 0
    minFreq = 20
    mag_hi_list = []
    mag_all_list = []
    updated = False
    window = np.hamming(blockSize)

    i = 0
    # split the audio into blocks of audio (ignore last block like matlab
    while (i + blockSize) < len(audio_samples):
        eval_audio = audio_samples[i:i + blockSize]
        complex_spectrum = np.fft.fft(eval_audio * window)
        magnitude_spectrum = np.absolute(complex_spectrum[0:1 + len(complex_spectrum) / 2])

        if sum(magnitude_spectrum) > 0:

            if not updated:
                fls = calcFrequencyScale(len(magnitude_spectrum), fs, blockSize)
                crossover_idx = calcCrossoverIDX(fls, crossover)
                minIDX = calcMinIDX(fls, minFreq)
                maxIDX = calcMaxIDX(fls, n_oct)
                if n_oct > 0:
                    f_upper = calcUpperFrequencyLimit(fls, n_oct, maxIDX)
                    f_lower = calcLowFrequencyLimit(fls, n_oct, maxIDX)

                updated = True

            if n_oct > 0:
                smoothed_spectrum = spectralSmoothing(magnitude_spectrum, f_upper, f_lower)
            else:
                smoothed_spectrum = magnitude_spectrum

            tpower = sum(smoothed_spectrum[minIDX:])

            # calculate the spectral centroid
            if tpower > threshold:
                upper_spectrum = smoothed_spectrum[crossover_idx:]
                upper_fls = fls[crossover_idx:(crossover_idx + len(upper_spectrum))]
                upper_power = sum(upper_spectrum)
                centroid = sum(upper_spectrum * upper_fls) / upper_power

                centroid_list.append(centroid)
                mag_all_list.append(tpower)
                mag_hi_list.append(upper_power)

        else:
            centroid_list.append(0)
            mag_all_list.append(0)
            mag_hi_list.append(0)

        i += stepSize

    if sum(mag_all_list) == 0:
        return 0

    mean_centroid = np.mean(np.array(centroid_list))
    mean_mag_hi = np.mean(np.array(mag_hi_list))
    mean_mag_all = np.mean(np.array(mag_all_list))

    # float required for float division in Python 2.7
    ratio = mean_mag_hi / float(mean_mag_all)

    # equation taken directly from Pearce [2016]
    bright = -25.8699 + (64.0127 * (np.log10(ratio) + (0.44 * np.log10(mean_centroid))))

    return bright
