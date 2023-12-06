from __future__ import division
import soundfile as sf
import numpy as np
import six
from . import timbral_util, timbral_hardness, timbral_depth, timbral_brightness, timbral_roughness, timbral_warmth, \
    timbral_sharpness, timbral_booming, timbral_reverb

def timbral_extractor(fname, fs=0, dev_output=False, phase_correction=False, clip_output=False, output_type='dictionary', verbose=True):
    """
      The Timbral Extractor will extract all timbral attribute sin one function call, returning the results as either
      a list or dictionary, depending on input definitions.

      Version 0.4

      Simply calls each function with

      Required parameter
      :param fname:             string or numpy array
                                string, audio filename to be analysed, including full file path and extension.
                                numpy array, array of audio samples, requires fs to be set to the sample rate.

     Optional parameters
      :param fs:                int/float, when fname is a numpy array, this is a required to be the sample rate.
                                Defaults to 0.
      :param phase_correction:  bool, perform phase checking before summing to mono.  Defaults to False.
      :param dev_output:        bool, when False return the depth, when True return all extracted
                                features.  Default to False.
      :param clip_output:             bool, force the output to be between 0 and 100.
      :param output_type:       string, defines the type the output should be formatted in.  Accepts either
                                'dictionary' or 'list' as parameters.  Default to 'dictionary'.

      :return: timbre           the results from all timbral attributes as either a dictionary or list, depending
                                on output_type.

      Copyright 2019 Andy Pearce, Institute of Sound Recording, University of Surrey, UK.
    """
    '''
      Check output_type before calculating anything
    '''
    if output_type != 'dictionary' and output_type != 'list':
        raise ValueError('output_type must be \'dictionary\' or \'list\'.')

    '''
      Basic audio reading
    '''
    if isinstance(fname, six.string_types):
        # read audio file only once and pass arrays to algorithms
        try:
            audio_samples, fs = sf.read(fname)
            # making an array again for copying purposes
            multi_channel_audio = np.array(audio_samples)
        except:
            print('Soundfile failed to load: ' + str(fname))
            raise TypeError('Unable to read audio file.')
    elif hasattr(fname, 'shape'):
        if fs==0:
            raise ValueError('If giving function an array, \'fs\' must be specified')
        audio_samples = fname
        multi_channel_audio = np.array(fname)
    else:
        raise ValueError('Input must be either a string or a numpy array.')

    # channel reduction
    audio_samples = timbral_util.channel_reduction(audio_samples)

    # resample audio file if sample rate is less than 44100
    audio_samples, fs = timbral_util.check_upsampling(audio_samples, fs)

    # functions can be given audio samples as well
    if verbose:
        print('Calculating hardness...')
    hardness = timbral_hardness(audio_samples, fs=fs,
                                dev_output=dev_output,
                                phase_correction=phase_correction,
                                clip_output=clip_output)
    if verbose:
        print('Calculating depth...')
    depth = timbral_depth(audio_samples, fs=fs,
                          dev_output=dev_output,
                          phase_correction=phase_correction,
                          clip_output=clip_output)
    if verbose:
        print('Calculating brightness...')
    brightness = timbral_brightness(audio_samples, fs=fs,
                                    dev_output=dev_output,
                                    phase_correction=phase_correction,
                                    clip_output=clip_output)
    if verbose:
        print('Calculating roughness...')
    roughness = timbral_roughness(audio_samples, fs=fs,
                                  dev_output=dev_output,
                                  phase_correction=phase_correction,
                                  clip_output=clip_output)
    if verbose:
        print('Calculating warmth...')
    warmth = timbral_warmth(audio_samples, fs=fs,
                            dev_output=dev_output,
                            phase_correction=phase_correction,
                            clip_output=clip_output)
    if verbose:
        print('Calculating sharpness...')
    sharpness = timbral_sharpness(audio_samples, fs=fs,
                                  dev_output=dev_output,
                                  phase_correction=phase_correction,
                                  clip_output=clip_output)
    if verbose:
        print('Calculating boominess...')
    boominess = timbral_booming(audio_samples, fs=fs,
                                dev_output=dev_output,
                                phase_correction=phase_correction,
                                clip_output=clip_output)
    if verbose:
        print('Calculating reverb...')
    # reverb calculated on all channels
    reverb = timbral_reverb(multi_channel_audio, fs=fs)

    '''
      Format output
    '''
    if output_type=='dictionary':
        timbre = {
            'hardness': hardness,
            'depth': depth,
            'brightness': brightness,
            'roughness': roughness,
            'warmth': warmth,
            'sharpness': sharpness,
            'boominess': boominess,
            'reverb': reverb
        }
    elif output_type == 'list':
        timbre = [hardness, depth, brightness, roughness, warmth, sharpness, boominess, reverb]
    else:
        raise ValueError('output_type must be \'dictionary\' or \'list\'.')


    return timbre



