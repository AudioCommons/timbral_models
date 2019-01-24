function thresholddB = ath_terhardt(f)
%ATH_TERHARDT Absolute Threshold of Hearing (ATH) by Terhardt (1979)
%   Frequency-dependent Absolute Threshold of Hearing
%   Equation given in 1979 Terhardt paper, Hearing Research (Eq. 1)
%   f: frequency (Hz)
%   thresholddB: ATH (dB SPL)

thresholddB = 3.64*(f/1000).^(-0.8) - 6.5*exp(-0.6*((f/1000-3.3).^2)) + ...
    1e-3*((f/1000).^4);


end

