classdef stftProc < Processor
    %stftProc Processor providing STFT of a time-domain ear signal.
    %
    % Returns the complex STFT time frequency spectrum per frame. To be more cache
    % friendly the output will have single precision. There are options to remove
    % the negative frequency bins (symmetry for real inputs).
    %
    
    %% Properties
    properties (Dependent = true)
        wSizeSec        % window size in seconds
        isPruned        % if true, remove DC and negative frequencies from output
    end
    
    properties (GetAccess = protected)
        wSize       % window size in samples
        hSize       % window shift in samples
        nFFT        % # of FFT frequency bins
        win         % window vector
        buffer      % unframed remainder of input
    end
    
    properties (SetAccess = protected)
        cfHz        % fft frequencies for output annotation
    end
    
    %% Methods
    methods
        
        function pObj = stftProc(fs, parObj)            
            if nargin<2 || isempty(parObj); parObj = Parameters; end
            if nargin<1; fs = []; end
            pObj = pObj@Processor(fs, fs, 'stftProc', parObj);            
            
            if nargin>0
                pObj.buffer = [];
            end            
        end
        
        function out = processChunk(pObj, in)            
            % prepend buffer content to data
            if ~isempty(pObj.buffer)
                in = [pObj.buffer; in];
            end
            
            % prepare frames and storage
            nSamples = size(in, 1);
            nFrames = floor((nSamples-(pObj.wSize-pObj.hSize))/pObj.hSize);
            out = zeros(nFrames, pObj.nFFT);
            
            % STFT over frames
            for ii = 1:nFrames
                % frame start and end indices
                fr_start = (ii-1)*pObj.hSize+1;
                fr_end = (ii-1)*pObj.hSize+pObj.wSize;
                % frame data with windowing function
                frame = pObj.win.*in(fr_start:fr_end);
                % compute frame FFT
                out(ii,:) = fft(frame, pObj.nFFT);
            end
            if pObj.isPruned
                out = out(:,1:pObj.nFFT/2);
            end
            % make sure output is single precision
            out = single(out);
            
            % save remaining time samples to buffer
            pObj.buffer = in(nFrames*pObj.hSize+1:end,:);            
        end
        
        function reset(pObj)
            pObj.buffer = [];
        end
        
    end
    
    %% "Overridden" methods
    methods (Hidden = true)
        
        function prepareForProcessing(pObj)            
            pObj.wSize = 2*round(pObj.parameters.map('stft_wSizeSec')*pObj.FsHzIn/2);
            pObj.hSize = pObj.wSize/2;
            pObj.win = hamming(pObj.wSize,'periodic');
            pObj.nFFT = 2^nextpow2(2*pObj.wSize-1);
            pObj.FsHzOut = 2/(pObj.wSizeSec);
            
            % freqz
            if pObj.isPruned
                % returns nFFT/2-1 bins, pruning symmetry and DC
                pObj.cfHz = (1:pObj.nFFT/2) * pObj.FsHzIn/pObj.nFFT;
            else
                % returns all frequencies
            	pObj.cfHz = [ (0:pObj.nFFT/2) ((pObj.nFFT+1):-1:1) ] * pObj.FsHzIn/pObj.nFFT;
            end            
        end
        
    end
    
    %% "Getter" methods
    methods
        
        function rval = get.wSizeSec(pObj)
            rval = pObj.parameters.map('stft_wSizeSec');
        end
        
        function rval = get.isPruned(pObj)
            rval = pObj.parameters.map('stft_isPruned');
        end
        
    end
    
    %% Static methods
    methods (Static)
        
        function dep = getDependency()            
            dep = 'time';            
        end
        
        function [names, defaultValues, descriptions] = getParameterInfo()            
            names = {...
                'stft_wSizeSec',...
                'stft_isPruned'};
                        
            descriptions = {...
                'Window duration (s)',...
                'if true, remove DC and negative frequencies from output'};
            
            defaultValues = {...
                20E-3,...
                true};
        end
        
        function pInfo = getProcessorInfo()            
            pInfo = struct;            
            pInfo.name = 'STFT';
            pInfo.label = 'STFT';
            pInfo.requestName = 'stft';
            pInfo.requestLabel = 'Short-Time-Fourier-Transform';
            pInfo.outputType = 'STFTSignal';
            pInfo.isBinaural = false;            
        end
        
    end
    
end
