from __future__ import division
import numpy as np
import librosa
import soundfile as sf
from scipy.signal import spectrogram
from . import timbral_util


def timbral_hardness(fname, dev_output=False, phase_correction=False, clip_output=False, max_attack_time=0.1,
                     bandwidth_thresh_db=-50):
    """
     This function calculates the apparent hardness of an audio file.

     Version 0.2

      Required parameters
      :param fname:                 Audio filename to be analysed, including full file path and extension.

      Optional parameters
      :param dev_output:            Bool, when False return the hardness, when True return all extracted features.
      :param phase_correction:      If the inter-channel phase should be estimated when performing a mono sum.
                                    Defaults to False.
      :param max_attack_time:       Set the maximum attack time, in seconds.  Defaults to 0.1.
      :param bandwidth_thresh_db:   Set the threshold for calculating the bandwidth, Defaults to -50dB.


      :return:                      Apparent hardness of audio file, float (dev_output = False/default).
                                    With dev_output set to True returns the weighted mean bandwidth,
                                    mean attack time, harmonic-percussive ratio, and unitless attack centroid.

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

    # reduce to mono
    audio_samples = timbral_util.channel_reduction(audio_samples, phase_correction)

    # normalise audio
    audio_samples /= max(abs(audio_samples))

    '''
      Calculate the harmonic-percussive ratio pre zero-padding the signal
    '''
    HP_ratio = timbral_util.get_percussive_audio(audio_samples, return_ratio=True)

    # zero pad the signal
    nperseg = 4096  # default value for spectrogram analysis
    audio_samples = np.lib.pad(audio_samples, (nperseg+1, 0), 'constant', constant_values=(0.0, 0.0))

    # calculate the envelope of the signal
    envelope = timbral_util.sample_and_hold_envelope_calculation(audio_samples, fs, decay_time=0.1)
    envelope_time = np.arange(len(envelope)) / fs

    # calculate the onsets
    original_onsets = timbral_util.calculate_onsets(audio_samples, envelope, fs, nperseg=nperseg)
    onset_strength = librosa.onset.onset_strength(audio_samples, fs)

    '''
      Calculate the spectrogram so that the bandwidth can be created
    '''
    mag = timbral_util.db2mag(bandwidth_thresh_db)  # calculate threshold in linear from dB
    bandwidth, t, f = timbral_util.get_bandwidth_array(audio_samples, fs, nperseg=nperseg, overlap_step=32,
                                                       rolloff_thresh=mag)

    '''
      Set all parameters for holding data per onset
    '''
    all_bandwidth_max = []
    all_attack_time = []
    all_max_strength = []
    all_max_strength_bandwidth = []
    all_attack_unitless_centroid = []

    '''
      Get bandwidth onset times and max bandwidth
    '''
    # If onsets don't exist, set it to time zero
    if not original_onsets:
        original_onsets = [0]

    onsets = np.array(original_onsets) - nperseg
    onsets[onsets < 0] = 0
    bandwidth_onset = np.array(onsets / 32.0).astype('int')  # overlap_step=32

    # itterate through each onset
    for onset_count in range(len(bandwidth_onset)):
        '''
          Calculate the bandwidth max for the attack portion of the onset
        '''
        # get the section of the bandwidth array between onsets
        onset = bandwidth_onset[onset_count]
        if onset == bandwidth_onset[-1]:
            bandwidth_seg = np.array(bandwidth[onset:])
        else:
            next_onset = bandwidth_onset[onset_count + 1]
            bandwidth_seg = np.array(bandwidth[onset:next_onset])

        # making a copy of the bandqwidth segment to avoid array changes
        hold_bandwidth_seg = list(bandwidth_seg)

        # bandwidth sample rate
        bandwidth_fs = fs / 32.0  # fs due to spectrogram step size

        # calculate onset of the attack in the bandwidth array
        if max(bandwidth_seg) > 0:
            bandwidth_attack = timbral_util.calculate_attack_time(bandwidth_seg, bandwidth_fs,
                                                                  calculation_type='fixed_threshold',
                                                                  max_attack_time=max_attack_time)
        else:
            bandwidth_attack = []

        # calculate the badiwdth max for the attack portion
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
          Calculate the attack time
        '''
        onset = original_onsets[onset_count]
        if onset == original_onsets[-1]:
            attack_seg = np.array(envelope[onset:])
            strength_seg = np.array(onset_strength[int(onset/512):])
            audio_seg = np.array(audio_samples[onset:])
        else:
            attack_seg = np.array(envelope[onset:original_onsets[onset_count + 1]])
            strength_seg = np.array(onset_strength[int(onset / 512):int(original_onsets[onset_count+1] / 512)])
            audio_seg = np.array(audio_samples[onset:original_onsets[onset_count + 1]])

        attack_time = timbral_util.calculate_attack_time(attack_seg, fs, max_attack_time=max_attack_time)
        all_attack_time.append(attack_time[0])

        '''
          Get the attack strength for weighting the bandwidth max
        '''
        all_max_strength.append(max(strength_seg))
        if bandwidth_attack:
            all_max_strength_bandwidth.append(max(strength_seg))

        '''
          Get the spectral centroid of the attack
        '''
        # identify the start of the attack
        th_start_idx = attack_time[2]
        # define how long the attack time can be
        centroid_int_samples = int(max_attack_time * fs)  # number of samples for attack time integration

        # start of attack section from attack time calculation
        if th_start_idx + centroid_int_samples >= len(audio_seg):
            audio_seg = audio_seg[th_start_idx:]
        else:
            audio_seg = audio_seg[th_start_idx:th_start_idx + centroid_int_samples]

        # get all spectral features for this attack section
        spectral_features_hold = timbral_util.get_spectral_features(audio_seg, fs)

        # store unitless attack centroid if exists
        if spectral_features_hold:
            all_attack_unitless_centroid.append(spectral_features_hold[2])


    # calculate mean values of parameters
    mean_attack_time = np.mean(all_attack_time)
    mean_attack_unitless_centroid = np.mean(all_attack_unitless_centroid)

    # get the weighted mean
    mean_weighted_bandwidth_max = np.average(all_bandwidth_max, weights=all_max_strength_bandwidth)

    if dev_output:
        # return the raw features
        return mean_weighted_bandwidth_max, mean_attack_time, HP_ratio, mean_attack_unitless_centroid
    else:
        '''
         Apply regression model
        '''
        all_metrics = np.ones(5)
        all_metrics[0] = mean_weighted_bandwidth_max
        all_metrics[1] = mean_attack_time
        all_metrics[2] = HP_ratio
        all_metrics[3] = mean_attack_unitless_centroid

        coefficients = np.array([0.0018878390275408967, -3.909514490829272, 8.633397251921414, 0.0422155929894954,
                                 27.045935123823075])

        hardness = np.sum(all_metrics * coefficients)

        if clip_output:
            hardness = timbral_util.output_clip(hardness)

        return hardness
