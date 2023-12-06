from __future__ import division
import numpy as np
import soundfile as sf
import six
from scipy.signal import spectrogram
from . import timbral_util

def timbral_reverb(fname, fs=0, dev_output=False, phase_correction=False, clip_output=False):
    """
     This function classifies the audio file as either not sounding reverberant.

     This is based on the RT60 estimation algirhtm documented in:
     Jan, T., and Wang, W., 2012: "Blind reverberation time estimation based on Laplace distribution",
     EUSIPCO. pp. 2050-2054, Bucharest, Romania.

     Version 0.4

     Required parameter
       :param fname:                 string or numpy array
                                     string, audio filename to be analysed, including full file path and extension.
                                     numpy array, array of audio samples, requires fs to be set to the sample rate.

     Optional parameters
      :param fs:                     int/float, when fname is a numpy array, this is a required to be the sample rate.
                                     Defaults to 0.
      :param phase_correction:       Has no effect on the code.  Implemented for consistency with other timbral
                                     functions.
      :param dev_output:             Has no effect on the code.  Implemented for consistency with other timbral
                                     functions.
      :param clip_output:            Has no effect on the code.  Implemented for consistency with other timbral
                                     functions.

      :return:                       predicted reverb of audio file.  1 represents the files osunds reverberant, 0
                                     represents the files does not sound reverberant.

      Copyright 2019 Andy Pearce, Institute of Sound Recording, University of Surrey, UK.
    """
    # needs to accept the input as audio file
    raw_audio_samples, fs = timbral_util.file_read(fname, fs=fs, phase_correction=False, mono_sum=False, loudnorm=False)

    # check for mono file
    if len(raw_audio_samples.shape) < 2:
        # it's a mono file
        mean_RT60 = estimate_RT60(raw_audio_samples, fs)
    else:
        # the file has channels, estimate RT for the first two and take the mean
        l_RT60 = estimate_RT60(raw_audio_samples[:, 0], fs)
        r_RT60 = estimate_RT60(raw_audio_samples[:, 1], fs)

        mean_RT60 = np.mean([l_RT60, r_RT60])

    '''
      need to develop a logistic regression model to test this.
    '''
    probability = reverb_logistic_regression(mean_RT60)

    if dev_output:
        return mean_RT60, probability
    else:
        if probability < 0.5:
            return 0
        else:
            return 1


def estimate_RT60(audio_samples, fs):

    ''' No chanel rediuction, perform on each channel '''

    # function[rt_est, par] = RT_estimation_my(y, fs)

    # performs blind RT estimation
    # INPUT
    # y: reverberant speech
    # fs: sampling frequency
    #
    # OUTPUT
    # rt_est: estimated RT
    # par: struct with parameters used to execute the function
    # rt_estimate_frame_my.m
    #
    # Codes were adapted from the original codes by Heinrich Loellmann, IND, RWTH Aachen
    #
    # Authors: Tariqullah Jan, moderated by Wenwu Wang, University of Surrey(2012)

    ''' 
     Initialisation 
    '''
    # ---------------------------------------------

    par = init_rt_estimate_e(fs)  # struct with all parameters and buffers for frame-wise processing
    BL = par['N'] * par['down']  # to simplify notation

    Laudio = len(audio_samples)

    # check audio file is long enough for analysis
    if BL < Laudio:
        rt_est = [] #np.zeros(int(round(Laudio / par['N_shift'])))
        RT_final = []

        '''
         frame-wise processing in the time-domain
        '''
        # ---------------------------------------------

        k = 0
        n_array = np.arange(0, Laudio - BL + 1, par['N_shift'])
        for n in n_array:
            k += 1  # frame counter
            ind = np.arange(n, n + BL) # indices of current frame

            # Actual RT estimation
            RT, par, finalrt = rt_estimate_frame_my(audio_samples[ind[np.arange(0, len(ind), par['down'])]], par)

            rt_est.append(RT) # store estimated value
            RT_final.append(finalrt)
    else:
        # audio too short for analysis, for returning smallest Rt value
        return par['Tquant'][0]

    RT_final = np.clip(RT_final, 0, max(RT_final))
    aaa = RT_final[np.where(RT_final>0)]


    RT_temp_new = []
    for i in range(1, len(aaa)):
        RT_temp_new.append(0.49 * aaa[i - 1] + (1 - 0.49) * np.max(aaa))


    if aaa.size:
        RTfinal_value = np.min(aaa)

        RT_temp_new = []
        for i in range(1, len(aaa)):
            RT_temp_new.append(0.49 * aaa[i - 1] + (1 - 0.49) * np.max(aaa))

    else:
        RTfinal_value = par['Tquant'][0]

    rt_est = np.array(rt_est)
    rt_est = rt_est[np.where(RT_final>0)]
    if rt_est.size:
        return np.mean(rt_est)
    else:
        return par['Tquant'][0]


def init_rt_estimate_e(fs=24000):
    '''
     par = init_rt_estimate_e(fs)
     executes initialization for the function
     rt_estimate_frame.m to perform a blind estimation of the reverberation time
     (RT) by frame-wise processing in the time-domain.

     INPUT
     fs: sampling frequency(default=24 kHz)

     OUTPUT
     par: struct containing all parameters and buffer for executing the
     function rt_estimate_frame.m

     author: Heiner Loellmann, IND, RWTH Aachen University

     created: August 2011

     general paraemters
    '''
    par = {"fs":fs}
    no = par['fs'] / 24000.0 # correction factor to account for different sampling frequency

    # pararmeters for pre - selection of suitable segments
    if par['fs'] > 8e3:
        par['down'] = 2  # rate for downsampling applied before RT estimation to reduce computational complexity
    else:
        par['down'] = 1

    par['N_sub'] = int(round(no * 700 / par['down']))  # sub-frame length(after downsampling)
    par['N_shift'] = int(round(no * 200 / par['down']))  # frame shift(before downsampling)
    par['nos_min'] = 3  # minimal number of subframes to detect a sound decay
    par['nos_max'] = 7  # maximal number of subframes to detect a sound decay
    par['N'] = int(par['nos_max'] * par['N_sub'])  # maximal frame length(after downsampling)

    # parameters for ML - estimation
    Tmax = 1.1  # max RT being considered
    Tmin = 0.2  #min RT being considered
    par['bin'] = 0.1  # step-size for RT estimation
    par['Tquant'] = np.arange(Tmin, Tmax+par['bin']/2, par['bin']) # set of qunatized RTs considered for maximum search
    par['a'] = np.exp(-3.0 * np.log(10) / ( par['Tquant'] * (par['fs'] / par['down'])))  # corresponding decay rate factors
    par['La'] = len(par['a'])  # num of considered decay rate factors( = no of.RTs)

    # paramters for histogram - based approach to reduce outliers (order statistics)
    par['buffer_size'] = int(round(no * 800 / par['down']))  # buffer size
    par['buffer'] = np.zeros(par['buffer_size']) # buffer with previous indices to update histogram
    par['no_bins'] = int(par['La'])  # no. of histogram bins
    par['hist_limits'] = np.arange(Tmin - par['bin'] / 2.0,  Tmax + par['bin'], par['bin'])  # limits of histogram bins
    par['hist_rt'] = np.zeros(par['no_bins'])  # histogram with ML estimates
    par['hist_counter'] = 0  # counter increased if histogram is updated

    # paramters for recursive smoothing of final RT estimate
    par['alpha'] = 0.995  # smoothing factor
    par['RT_initial'] = 0.3  # initial RT estimate
    par['RT_last'] = par['RT_initial']  # last RT estimate
    par['RT_raw'] = par['RT_initial']  # raw RT estimate obtained by histogram - approach

    return par


def rt_estimate_frame_my(frame, par):
    '''
     performs an efficient blind estimation of the reverberation time(RT) for frame-wise
     processing based on Laplacian distribution.

     INPUT
     frame: (time-domain) segment with reverberant speech
     par: struct with all parameters and buffers created by the function
     init_binaural_speech_enhancement_e.m

     OUTPUT
     RT: estimated RT
     par: struct with updated buffers to enable a frame-wise processing
     RT_pre: raw RT estimate(for debugging and analysis of the algorithm)

     Reference:  LAllmann, H.W., Jeub, M., Yilmaz, E., and Vary, P.:
     An Improved Algorithm for Blind Reverberation Time Estimation, a
     International Workshop on Acoustic Echo and Noise Control(IWAENC), Tel Aviv, Israel, Aug. 2010.

     Tariqullah Jan and Wenwu Wang:
     Blind reverberation time estimation based on Laplacian distribution
     European Signal Processing Conference(EUSIPCO), 2012.

     The codes were adapted based on the original codes by Heinrich Loellmann, IND, RWTH Aachen

     Authors: Tariqullah Jan, moderated by Wenwu Wang, University of Surrey(2012)
    '''
    if len(np.shape(np.squeeze(frame))) > 1:
        raise ValueError('Something went wrong...')

    cnt = 0  # sub-frame counter for pre - selection of possible sound decay
    RTml = -1  # default RT estimate (-1 indicates no new RT estimate)

    # calculate variance, minimum and maximum of first sub-frame
    seg = frame[:par['N_sub']]

    var_pre = np.var(seg)
    min_pre = np.min(seg)
    max_pre = np.max(seg)

    for k in range(2, par['nos_max']):
        # calculate variance, minimum and maximum of succeding sub-frame
        seg = frame[(k-1) * par['N_sub'] : k * par['N_sub']+1]
        var_cur = np.var(seg)
        max_cur = max(seg)
        min_cur = min(seg)

        #-- Pre-Selection of suitable speech decays --------------------
        if (var_pre > var_cur) and (max_pre > max_cur) and (min_pre < min_cur):
            # if variance, maximum decraease, and minimum increase
            # = > possible sound decay detected

            cnt += 1

            # current values becomes previous values
            var_pre = var_cur
            max_pre = max_cur
            min_pre = min_cur

        else:
            if cnt >= par['nos_min']:
                # minimum length for assumed sound decay achieved?
                # -- Maximum Likelihood(ML) Estimation of the RT
                RTml, _ = max_loglf(frame[:cnt*par['N_sub']], par['a'], par['Tquant'])

            break


        if k == par['nos_max']:
            # maximum frame length achieved?
            RTml, _ = max_loglf(frame[0:cnt * par['N_sub']], par['a'], par['Tquant'])

    # end of sub-frame loop

    if RTml >= 0:  # new ML estimate calculated

        # apply order statistics to reduce outliers
        par['hist_counter'] += 1

        for i in range(par['no_bins']):

            # find index corresponding to the ML estimate
            # find index corresponding to the ML estimate
            if (RTml >= par['hist_limits'][i]) and (RTml <= par['hist_limits'][i+1]):

                index = i
                break

        # update histogram with ML estimates for the RT
        par['hist_rt'][index] += 1

        if par['hist_counter'] > par['buffer_size'] + 1:
            # remove old values from histogram
            par['hist_rt'][int(par['buffer'][0])] = par['hist_rt'][int(par['buffer'][0])] - 1

        par['buffer'] = np.append(par['buffer'][1:], index)  # % update buffer with indices
        idx = np.argmax(par['hist_rt'])  # find index for maximum of the histogram

        par['RT_raw'] = par['Tquant'][idx]  # map index to RT value


    # final RT estimate obtained by recursive smoothing
    RT = par['alpha'] * par['RT_last'] + (1 - par['alpha']) * par['RT_raw']
    par['RT_last'] = RT

    RT_pre = RTml  # intermediate ML estimate for later analysis

    return RT, par, RT_pre


def max_loglf(h, a, Tquant):
    '''
     [ML, ll] = max_loglf(h, a, Tquant)

     returns the maximum of the log-likelihood(LL) function and the LL
     function itself for a finite set of decay rates

     INPUT
     h: input frame
     a: finite set of values for which the max.should be found
     T: corresponding RT values for vector a

     OUTPUT
     ML: ML estimate for the RT
     ll: underlying LL - function
    '''

    N = len(h)
    n = np.arange(0, N)  # indices for input vector
    ll = np.zeros(len(a))

    # transpose?
    h_square = h.transpose()

    for i in range(len(a)):
        sum1 = np.dot((a[i] ** (-1.0 * n)), np.abs(h_square))
        sum2 = np.sum(np.abs(h_square))
        sigma = (1 / N) * sum1
        ll[i] = -N * np.log(2) - N * np.log(sigma) - np.sum(np.log(a[i] ** n)) - (1 / sigma) * sum1


    idx = np.argmax(ll)  # maximum of the log-likelihood function
    ML = Tquant[idx]  # corresponding ML estimate for the RT

    return ML, ll


def reverb_logistic_regression(mean_RT60):
    """
      Logistic regression function to determine if the file sound reverberant or not.
    :param mean_RT60:
    :return:
    """
    # apply linear coefficients
    coefficients = [2.97126461]
    intercept = -1.45082989
    attributes = [mean_RT60]
    logit_model = np.sum(np.array(coefficients) * np.array(attributes)) + intercept

    # apply inverse of Logit function to obtain probability
    probability = np.exp(logit_model) / (1.0 + np.exp(logit_model))

    return probability