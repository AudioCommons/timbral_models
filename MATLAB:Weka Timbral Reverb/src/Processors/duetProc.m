classdef duetProc < Processor
%duetProc Processor providing weighted histogram of level/phase differences.
%
% DUET is a binaural blind source seperation algorithm based on a histogram
% of differences in level and phase of the time-frequency representation of
% the mixture signals. We use an intermediate result of the algorithm, the
% weighted histogram of R/L level and time differences, which shows clear
% peaks for individual sources. This method relies on the W-disjunct orthogonality
% of the sources (e.g. only one source is active per time-frequency bin). This
% is assumption is (mildly) violated for echoic mixtures.
%
% Actually it seems like phase information is severely distorted in the noisy reverberant
% case.
%
% [1] Rickard, S. (2007). The DUET Blind Source Separation Algorithm.
%     In S. Makino, H. Sawada, & T.-W. Lee (Eds.),
%     Blind Speech Separation (pp. 217–241). inbook,
%     Dordrecht: Springer Netherlands.
%     http://doi.org/10.1007/978-1-4020-6479-1_8
%
    
    %% Properties
    properties (Dependent = true)
        wSTFTSec;   % window size in seconds for STFT
        wDUETSec;   % window size in seconds for DUET histogram
        hDUETSec;   % window shift for duet historgram
        binsAlpha;  % number of histogram bins for alpha dimension
        binsDelta;  % number of histogram bins for delta dimension
        maxAlpha;   % masking threshold for alpha dimension
        maxDelta;   % masking threshold for delta dimension
    end
    
    properties (GetAccess = protected)
        % STFT properties
        wSFTT;          % STFT window size in samples
        hSTFT;          % STFT window shift in samples
        nSTFT;          % STFT FFT frequency bins
        windowSTFT;     % STFT window vector
        FsOutSTFT;      % sampling frequency STFT output
        % DUET properties
        wDUET;          % number of STFT frames per DUET histogram.
        hDUET;          % number of STFT frames per DUET histogram.
        % buffers
        bufferT_l;      % unframed remainder of time domain input, left
        bufferT_r;      % unframed remainder of time domain input, right
        bufferTF_l;     % unframed remainder of time-frequency domain input, left
        bufferTF_r;     % unframed remainder of time-frequency domain input, right
    end
    
    %% Methods
    methods
        
        function pObj = duetProc(fs, parObj)            
            if nargin<2 || isempty(parObj); parObj = Parameters; end
            if nargin<1; fs = []; end
            pObj = pObj@Processor(fs, fs, 'duetProc', parObj);            
            
            if nargin>0
                pObj.reset();
            end            
        end
        
        function out = processChunk(pObj, in_l, in_r)            
            % prepend time buffer contents
            if ~isempty(pObj.bufferT_l)
                in_l = [pObj.bufferT_l; in_l];
            end
            if ~isempty(pObj.bufferT_r)
                in_r = [pObj.bufferT_r; in_r];
            end
            if max(size(in_l)~=size(in_r))
                error('Buffered inputs should be of same dimension for both ears')
            end            
            
            % prepare time domain frames
            nSamplesT = size(in_l, 1);
            nFramesT = floor((nSamplesT-(pObj.wSFTT-pObj.hSTFT))/pObj.hSTFT);
            stft_l = zeros(nFramesT, pObj.nSTFT);
            stft_r = zeros(nFramesT, pObj.nSTFT);            
            % STFT over time domain frames
            for ii = 1:nFramesT
                % frame start and end indices
                fr_start = (ii-1)*pObj.hSTFT+1;
                fr_end = (ii-1)*pObj.hSTFT+pObj.wSFTT;
                % frame data
                frame_l = pObj.windowSTFT.*in_l(fr_start:fr_end);
                frame_r = pObj.windowSTFT.*in_r(fr_start:fr_end);
                % compute frame FFT
                stft_l(ii,:) = fft(frame_l, pObj.nSTFT);
                stft_r(ii,:) = fft(frame_r, pObj.nSTFT);
            end
            % time domain buffer
            pObj.bufferT_l = in_l(nFramesT*pObj.hSTFT+1:end,:);
            pObj.bufferT_r = in_l(nFramesT*pObj.hSTFT+1:end,:);
            % time-frequency domain buffer
            pObj.bufferTF_l = [pObj.bufferTF_l; stft_l];
            pObj.bufferTF_r = [pObj.bufferTF_r; stft_r];
            
            % prepare time-frequency domain frames
            nSamplesTF = size(pObj.bufferTF_l, 1);
            nFramesTF = floor((nSamplesTF-(pObj.wDUET-pObj.hDUET))/pObj.hDUET);
            out = zeros(nFramesTF, pObj.binsAlpha, pObj.binsDelta);            
            % DUET over time-frequency domain frames
            for ii = 1:nFramesTF
                % frame start and end indices
                fr_start = (ii-1)*pObj.hDUET+1;
                fr_end = (ii-1)*pObj.hDUET+pObj.wDUET;
                % frame data
                frame_l = pObj.bufferTF_l(fr_start:fr_end,:);
                frame_r = pObj.bufferTF_r(fr_start:fr_end,:);
                % compute frame FFT
                out(ii,:,:) = pObj.computeDUEThistogram(...
                    frame_l, frame_r,...
                    pObj.binsAlpha, pObj.binsDelta,...
                    pObj.maxAlpha, pObj.maxDelta);
            end
            % prune non-causal time-frequency input from buffer
            pObj.bufferTF_l(1:nFramesTF*pObj.hDUET,:) = [];
            pObj.bufferTF_r(1:nFramesTF*pObj.hDUET,:) = [];
        end
        
        function reset(pObj)
            pObj.bufferT_l = [];
            pObj.bufferT_r = [];
            pObj.bufferTF_l = [];
            pObj.bufferTF_r = [];
        end
        
    end
    
    %% "Overridden" methods
    methods (Hidden = true)
        
        function prepareForProcessing(pObj)            
            % STFT params: even window size with periodic hamming window q=2
            pObj.wSFTT = 2*round(pObj.wSTFTSec*pObj.FsHzIn/2);
            pObj.hSTFT = pObj.wSFTT/2;
            pObj.windowSTFT = hamming(pObj.wSFTT,'periodic');
            pObj.nSTFT = 2^nextpow2(2*pObj.wSFTT-1);
            pObj.FsOutSTFT = 2/pObj.wSTFTSec;
            % DUET internal parameters
            pObj.wDUET = floor(pObj.wDUETSec * pObj.FsOutSTFT);
            pObj.hDUET = floor(pObj.hDUETSec * pObj.FsOutSTFT);
            % output sampling rate
            pObj.FsHzOut = pObj.FsOutSTFT/pObj.hDUET;
        end
        
        function addInput(pObj,dependency)
            pObj.Input{1,1} = dependency{1}.Output{1};
            pObj.Input{1,2} = dependency{1}.Output{2};
        end
        
    end
    
    %% "Getter" methods
    methods
        
        function rval = get.wSTFTSec(pObj)
            rval = pObj.parameters.map('duet_wSTFTSec');
        end
        
        function rval = get.wDUETSec(pObj)
            rval = pObj.parameters.map('duet_wDUETSec');
        end
        
        function rval = get.hDUETSec(pObj)
            rval = pObj.parameters.map('duet_hDUETSec');
        end
        
        function rval = get.binsAlpha(pObj)
            rval = pObj.parameters.map('duet_binsAlpha');
        end
        
        function rval = get.binsDelta(pObj)
            rval = pObj.parameters.map('duet_binsDelta');
        end
        
        function rval = get.maxAlpha(pObj)
            rval = pObj.parameters.map('duet_maxAlpha');
        end
        
        function rval = get.maxDelta(pObj)
            rval = pObj.parameters.map('duet_maxDelta');
        end
        
    end
    
    %% Static methods
    methods (Static)
        
        function dep = getDependency()            
            dep = 'time';     
        end
        
        function [names, defaultValues, descriptions] = getParameterInfo()            
            names = {...
                'duet_wSTFTSec',...
                'duet_wDUETSec',...
                'duet_hDUETSec',...
                'duet_binsAlpha',...
                'duet_binsDelta',...
                'duet_maxAlpha',...
                'duet_maxDelta'};
                        
            descriptions = {...
                'STFT window duration (s)',...
                'DUET window duration (s)',...
                'DUET window shift (s)',...
                'number of histogram bins for alpha dimension',...
                'number of histogram bins for delta dimension',...
                'max alpha threshold',...
                'max delta threshold'};
            
            defaultValues = {...
                20E-3,...
                1/2,...
                1/6.,...
                35,...
                51,...
                0.7,...
                3.6};
        end
        
        function pInfo = getProcessorInfo()            
            pInfo = struct;            
            pInfo.name = 'DUET';
            pInfo.label = 'DUET histogram';
            pInfo.requestName = 'duet';
            pInfo.requestLabel = 'DUET histogram based on STFT';
            pInfo.outputType = 'DuetHistogramSignal';
            pInfo.isBinaural = true;    
        end
        
        function rval = computeDUEThistogram(...
                tf_l, tf_r,...
                bins_alpha, bins_delta, max_alpha, max_delta)
            %computeDUEThistogram   builds a duet histogram over the input
            %   returns the weighted, smoothed histogram
            %
            % INPUTS:
            %   tf_l,tf_r           : left and right ear time-frequency representation
            %   resampleFactor      : if ~= 1, down sample the output by that factor
            %   bins and thresholds : bins per axis and respective masking thresholds
            %
            
            % init
            if nargin < 2
                error('method needs at least the two time-frequency mixtures as input!');
            else
                tf_l = double(tf_l);
                tf_r = double(tf_r);
            end
            % default parameters from [1]
            if nargin < 3
                bins_alpha = 35;
            end
            if nargin < 4
                bins_delta = 51;
            end
            if nargin < 5
                max_alpha = 0.7;
            end
            if nargin < 6
                max_delta = 3.6;
            end
            
            
            % preprocess tf data
            [nFrames, nFFT] = size(tf_l);
            tf_l(:,1) = [];
            tf_r(:,1) = [];
            freqInfo = [ (1:nFFT/2) -((nFFT/2-1):-1:1) ] * (2*pi/nFFT);
            fmat = freqInfo(ones(nFrames,1),:);
            
            % DUET: estimation of alpha (power) and delta (delay) mixture parameters
            tf_rl = (tf_r + eps)./(tf_l + eps);
            % following projection is called symmetric attenuation in [1]
            alpha = abs(tf_rl) - 1./abs(tf_rl);
            % positive alpha means more power on the R channel
            % negative alpha means more power on the L channel
            delta = -imag(log(tf_rl))./fmat;
            % delta is the phase difference in radiants, means R phase earlier than L
            % negative delta means L phase earlier than R
            
            % DUET: calculate weighted histogram
            % weighting powers according to [1]
            % p=0; q=0; % simple counting
            % p=1; q=0; % only symetric attenuation
            % p=1; q=2; % emphasis on delays
            % p=2; q=0; % reducing bias on the delay estimator
            % p=2; q=2; % low SRN and speech mixtures
            % we settle for p=2 and q=0
            p=2; q=2;
            tf_weight = (abs(tf_l).*abs(tf_r)).^p.*abs(fmat).^q; % weights vector
            % mask tf-points yielding estimates in bounds
            mask = (abs(alpha)<max_alpha) & (abs(delta)<max_delta);
            v_alpha = alpha(mask);
            v_delta = delta(mask);
            tf_weight = tf_weight(mask);
            % determine histogram indices
            idx_alpha = double(round(1+(bins_alpha-1)*(v_alpha+max_alpha)/(2*max_alpha)));
            idx_delta = double(round(1+(bins_delta-1)*(v_delta+max_delta)/(2*max_delta)));
            % full sparse trick to create 2d weighted histogram
            duet_hist_raw = full(sparse(idx_alpha, idx_delta, tf_weight, bins_alpha, bins_delta));
            % DEBUG:
            %mesh(linspace(-max_delta, max_delta, bins_delta),...
            %     linspace(-max_alpha, max_alpha, bins_alpha),...
            %     duet_hist_raw);
            
            % salvitzky-golay smoothing filter
            duet_hist = sgolayfilt(duet_hist_raw, 3, 5);
            duet_hist = sgolayfilt(duet_hist', 3, 5)';
            % cut filter artifacts and use max scaling
            duet_hist(duet_hist < 0) = 0;
            duet_hist = duet_hist ./ max(max(duet_hist));
            rval = duet_hist;
        end
        
    end
    
end
