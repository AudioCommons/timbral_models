from __future__ import division
import numpy as np
import soundfile as sf
from . import timbral_util


def boominess_calculate(loudspec):
    """
      Calculates the Booming Index as described by Hatano, S., and Hashimoto, T. "Booming index as a measure for
      evaluating booming sensation", The 29th International congress and Exhibition on Noise Control Engineering, 2000.
    """

    # loudspec from the loudness_1991 code results in values from 0.1 to 24 Bark in 0.1 steps
    z = np.arange(0.1, 24.05, 0.1)  #0.1 to 24 bark in 0.1 steps
    f = 600 * np.sinh(z / 6.0)  # convert these bark values to frequency
    FR = [25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500,
          3150, 4000, 5000, 6300, 8000, 10000, 12500] # get the centre frequencies of 3rd octave bands

    # now I need to convert f onto the FR scale
    logFR = np.log10(FR)
    FR_step = logFR[1] - logFR[0]  # get the step size on the log scale
    FR_min = logFR[0]  # get the minimum value of the logFR

    logf = np.log10(f)  # get the log version of estimated frequencies
    # estimate the indexes of the bark scale on the 3rd octave scale
    estimated_index = ((logf - FR_min) / float(FR_step)) + 1

    # weighting function based from the estimated indexes
    Weighting_function = 2.13 * np.exp(-0.151 * estimated_index)

    # change the LF indexes to roll off
    Weighting_function[0] = 0.8  # this value is estimated
    Weighting_function[1] = 1.05
    Weighting_function[2] = 1.10
    Weighting_function[3] = 1.18

    # identify index where frequency is less than 280Hz
    below_280_idx = np.where(f >= 280)[0][0]

    I = loudspec * Weighting_function
    loudness = np.sum(loudspec)
    Ll = np.sum(loudspec[:below_280_idx])

    Bandsum = timbral_util.log_sum(I)
    BoomingIndex = Bandsum * (Ll / loudness)

    return BoomingIndex


def timbral_booming(fname, dev_output=False, phase_correction=False):
    """
     This is an implementation of the hasimoto booming index feature.
     There are a few fudge factors with the code to convert between the internal representation of the sound using the
     same loudness calculation as the sharpness code.  The equation for calculating the booming index is not
     specifically quoted anywhere so I've done the best i can with the code that was presented.

     Shin, SH, Ih, JG, Hashimoto, T., and Hatano, S.: "Sound quality evaluation of the booming sensation for passenger
      cars", Applied Acoustics, Vol. 70, 2009.

     Hatano, S., and Hashimoto, T. "Booming index as a measure for
      evaluating booming sensation", The 29th International congress and Exhibition on Noise Control Engineering, 2000.

     This function calculates the apparent Boominess of an audio file.
     This version of timbral_booming relates to D5.7.

     Version 0.3

     Required parameter
      :param fname:                   string, audio filename to be analysed, including full file path and extension.

     Optional parameters
      :param dev_output:              bool, when False return the warmth, when True return all extracted features.
                                      Defaults to False.
      :param phase_correction:        bool, if the inter-channel phase should be estimated when performing a mono sum.
                                      Defaults to False.

      :return                         float, apparent boominess of the audio file.

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
    audio_samples = timbral_util.channel_reduction(audio_samples, phase_correction=phase_correction)

    # window the audio file into 4096 sample sections
    windowed_audio = timbral_util.window_audio(audio_samples, window_length=4096)

    windowed_booming = []
    windowed_rms = []
    for i in range(windowed_audio.shape[0]):
        samples = windowed_audio[i, :]  # the current time window
        # get the rms value and append to list
        windowed_rms.append(np.sqrt(np.mean(samples * samples)))

        # calculate the specific loudness
        N_entire, N_single = timbral_util.specific_loudness(samples, Pref=100.0, fs=fs, Mod=0)

        # calculate the booming index is contains a level
        if N_entire > 0:
            # boom = boominess_calculate(N_single)
            BoomingIndex = boominess_calculate(N_single)
        else:
            BoomingIndex = 0

        windowed_booming.append(BoomingIndex)

    # get the rms-weighted average
    rms_boom = np.average(windowed_booming, weights=windowed_rms)

    if dev_output:
        return [rms_boom]
    else:
        return rms_boom
