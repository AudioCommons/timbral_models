classdef modulationProc < Processor
%MODULATIONPROC Amplitude modulation spectogram processor.
%   The Amplitude Modulation Spectrogram is derived by analizying the Inner
%   Hair-Cell representation at each frequency channel with a bank of
%   modulation filters, which mimics the envelope fluctuation detection of
%   the human auditory system.
%
%   MODULATIONPROC properties (note that the input parameters have different names):
%        modCfHz         - Modulation filter center frequencies (Hz)
%        filterType      - Filterbank type ('lin' [1,2] vs. 'log' [3,4])
%        lowFreqHz       - Lowest modulation center frequency 
%        highFreqHz      - Highest modulation center frequency 
%        winName         - Window name
%        stepSec         - Window step size in seconds
%        blockSec        - Window duration in seconds
%        dsRatio         - Down-sampling ratio of the IHC representation
%        nAudioChan      - Number of (IHC rep) audio frequency channels
%        nModChan        - Number of modulation filters
%
%   See also: Processor, ihcProc
%
%REFERENCES:
% 
%   [1] Kim, G., Lu, Y, Hu, Y. and Loizou, P. C. (2009), "An algorithm that
%       improves speech intelligibility in noise for normal-hearing
%       listeners," Journal of the Acoustical Society of America 126(3),
%       pp. 1486-1494. 
% 
%   [2] May, T. and Dau, T. (2014), "Requirements for the evaluation of 
%       computational speech segregation systems," Journal of the 
%       Acoustical Society of America 136(6), pp. EL398-EL404.
% 
%   [3] May, T. and Dau, T. (2014), "Computational speech segregation based
%      on an auditory-inspired modulation analysis," Journal of the
%      Acoustical Society of America 136(6), pp. 3350-3359. 
% 
%   [4] Ewert, S. D. and Dau, T. (2000), "Characterizing frequency 
%       selectivity for envelope fluctuations," Journal of the Acoustical 
%       Society of America 108(3), pp. 1181-1196.

    properties (SetAccess = private)
        modCfHz         % Modulation filters center frequencies
    end

    properties (Dependent = true)
        filterType      % 'lin' vs. 'log'
        
        winName         % Window
        stepSec         % Step size in seconds
        blockSec        % Block size in seconds
        
        dsRatio         % Down-sampling ratio
        nAudioChan      % Number of audio frequency channels
    end
    
    properties %(GetAccess = private)
        buffer          % Buffered input (for fft-based)
        
        lowFreqHz       % Lowest modulation center frequency 
        highFreqHz      % Highest modulation center frequency 
        
        overlap         % Overlap in samples
        blockSize       % Block size in samples

        % Downsampling
        dsProc          % Downsampling processor
        fs_ds           % Input sampling frequency after downsampling
        
        % For fft-based processing and framing in filter-based processing
        nfft            % FFT size
        win             % Window vector
        
        % For fft-based processing
        wts             % Sparse matrix for spectrogram mixing
        mn              % Index of lowest bin
        mx              % Index of higher bin
        
        % For filter-based processing
        b               % Cell array of filter coefficients
        a               % Cell array of filter coefficients
        bw              % Filters bandwidths (necessary?)
        Filters         % Cell array of filter objects
    end
    
    methods
        
        function pObj = modulationProc(fs,parObj)
		%modualationProc   Construct an amplitude modulation extraction processor
        %
        % USAGE:
        %   pObj = modulationProc(fs, parObj)
        %
        % INPUT ARGUMENTS:
        %     fs : Input sampling frequency (Hz)
        % parObj : Parameter object instance
        %
        % OUTPUT ARGUMENTS:
        %   pObj : Processor instance
        %
        % NOTE: Parameter object instance, parObj, can be generated using genParStruct.m
        % User-controllable parameters for this processor and their default values can be
        % found by browsing the script parameterHelper.m
        %
        % See also: genParStruct, parameterHelper, Processor
            
            
            % TODO:
            % - make standalone? (currently uses melbankm.m and createFB_Mod.m)
            % - signal normalization with long-term rms (how to integrate
            % it in online-processing?)
            % - envelope normalization?
            
            
            % Checking input parameter
            if nargin<2||isempty(parObj); parObj = Parameters; end
            if nargin<1; fs = []; end
            
            % Call superconstructor
            pObj = pObj@Processor(fs,[],'modulationProc',parObj);
            
            if nargin>0
                pObj.buffer = [];
            end
            
        end
        
        function out = processChunk(pObj,in)
            %processChunk       Requests the processing for a new chunk of
            %                   signal
            %
            %USAGE:
            %    out = processChunk(in)
            %
            %INPUT ARGUMENTS:
            %   pObj : Processor instance
            %     in : Input chunk
            %
            %OUTPUT ARGUMENT:
            %    out : Processor output for that chunk
            
            % Down-sample the input if needed
            if pObj.dsRatio > 1
                in = pObj.dsProc.processChunk(in);
            end
            
            if strcmp(pObj.filterType,'lin')
            
                % 1- Append the buffer to the input
                in = [pObj.buffer;in];  % Time spans the first dimension

                % 2- Initialize the output
                nbins = max(floor((size(in,1)-pObj.overlap)/(pObj.blockSize-pObj.overlap)),0);
                out = zeros(nbins,size(in,2),size(pObj.modCfHz,2));

                % 3- Process each frequency channel and store remaining buffer
                
                % Process if the input is long enough (spectrogram returns
                % an error if the input is shorter than one window)
                if nbins > 0

                    for ii = 1:size(in,2)

                        % Calculate the modulation pattern for this filter
                        ams = spectrogram(in(:,ii),pObj.win,pObj.overlap,pObj.nfft);

                        % Normalize spectrogram
                        ams = ams / pObj.nfft;
                        
                        % Restrain the pattern to the required mod. frequency bins
                        output = pObj.wts*abs(ams(pObj.mn:pObj.mx,:));

                        % Store it appropriately in the output
                        out(:,ii,:) = permute(output,[2 3 1]);

                        % Initialize a buffer for the first frequency channel
                        if ii == 1
                            % Initialize a temporary buffer for that chunk
                            % Buffer size might change between chunks, hence need
                            % to re-initialize it

                            % Indexes in the input of buffer start and end
                            bstart = size(output,2)*(length(pObj.win)-pObj.overlap)+1;
                            bend = size(in,1);

                            % Initialize a temporary buffer
                            temp_buffer = zeros(bend-bstart+1,size(in,2));
                        end

                        % Store the buffered input for that channel
                        temp_buffer(:,ii) = in(bstart:bend,ii);
                    end

                % If not, then buffer the all input signal    
                else
                    temp_buffer = in;
                end
                    
                % 4- Update the buffer from buffers collected in step 3
                pObj.buffer = temp_buffer;
                
            
            elseif strcmp(pObj.filterType,'log')
                
                % Initialize the output
                nbins = floor(((size(in,1)+size(pObj.buffer,1))-pObj.overlap)/(pObj.blockSize-pObj.overlap));
                out = zeros(nbins,size(in,2),size(pObj.modCfHz,2));

                % Process each frequency channel
                for ii = 1:size(in,2)
                    
                    % Calculate the modulation pattern for this audio filter
                    % Loop over number of modulation filter
                    for jj = 1:numel(pObj.modCfHz)

                        % Calculate AMS pattern of jj-th filter
                        currAMS = pObj.Filters((ii-1)*numel(pObj.modCfHz)+jj).filter(in(:,ii));
                            
                        % Append the buffer to the ams (TODO: might want to
                        % move the isempty test out of the loops)
                        if ~isempty(pObj.buffer)
                            currAMS = [pObj.buffer(:,(ii-1)*numel(pObj.modCfHz)+jj);currAMS];
                        end
                        
                        % Frame-based analysis...
                        out(:,ii,jj) = mean(abs(frameData(currAMS,pObj.blockSize,pObj.blockSize-pObj.overlap,pObj.win,false)));

                        % Initialize a temporary buffer
                        if (ii==1) && (jj==1)
                            bstart = size(out,1)*(length(pObj.win)-pObj.overlap)+1;
                            bend = size(currAMS,1);
                            temp_buffer = zeros(bend-bstart+1,size(in,2)*numel(pObj.modCfHz));
                        end
                        
                        % Update the buffer for this audio and modulation
                        % frequencies
                        temp_buffer(:,(ii-1)*numel(pObj.modCfHz)+jj)=currAMS(bstart:bend);
                    end
                    
                end
                
                % Store the buffer
                pObj.buffer = temp_buffer;
                
            end
                
        end
        
        function reset(pObj)
            %reset      Resets the internal states of the processor
            %
            %USAGE:
            %    pObj.reset
            % 
            %INPUT ARGUMENT:
            %    pObj : Processor instance
           
            % Reset the buffer
            pObj.buffer = [];
            
            % Reset the filters states if needed
            if strcmp(pObj.filterType,'filter')
                for ii = 1:size(pObj.Filters)
                    pObj.Filters(ii).reset;
                end
            end
            
        end
        
    end
    
    methods (Access=protected)
        
        function verifyParameters(pObj)
           
            % Check inputs
            if mod(pObj.dsRatio,1)~=0 || pObj.dsRatio < 1
                error('The down-sampling ratio should be a positive integer')
            end
            
            if ~strcmp(pObj.filterType,'lin') && ~strcmp(pObj.filterType,'log')
                warning('%s is an invalid argument for modulation filterbank instantiation, switching to ''log''.')
                pObj.filterType = 'log';
            end
            
            if pObj.highFreqHz > pObj.fs_ds/2
                error('Highest modulation center frequency is above Nyquist frequency. Either reduce the downsampling ratio or decrease the upper frequency limit. ')
            end
        end
        
    end
    
    methods (Hidden = true)
        
        function prepareForProcessing(pObj)
            
            fs = pObj.FsHzIn;
            
            % Instantiate a down-sampler if needed
            if pObj.dsRatio > 1
                pObj.dsProc = downSamplerProc(fs,pObj.dsRatio,1);
            end

            % Compute internal parameters
            pObj.fs_ds  = fs / pObj.dsRatio;
            pObj.blockSize = 2 * round(pObj.blockSec * pObj.fs_ds / 2);

            pObj.lowFreqHz = pObj.parameters.map('ams_lowFreqHz');
            pObj.highFreqHz = pObj.parameters.map('ams_highFreqHz');

            % Get filterbank properties
            if strcmp(pObj.filterType,'lin')

                if isempty(pObj.parameters.map('ams_cfHz'))
                    % FFT-size
                    fftFactor = 2;  
                    pObj.nfft = pow2(nextpow2(fftFactor*pObj.blockSize));

%                         if isempty(lowFreqHz);
%                             lowFreqHz = 0;
%                         end
%                         if isempty(highFreqHz);
%                             highFreqHz = 400;
%                         end
%                         if isempty(nFilters)
%                             nFilters = 15;
%                         end

                    % Normalized lower and upper frequencies of the mod. filterbank
                    fLow  = pObj.lowFreqHz  / pObj.fs_ds;
                    fHigh = pObj.highFreqHz / pObj.fs_ds;
                    [pObj.wts,pObj.modCfHz,pObj.mn,pObj.mx] = ...
                        melbankm(pObj.parameters.map('ams_nFilters'), pObj.nfft, ...
                                 pObj.fs_ds,fLow,fHigh,'fs');
                else
                    error('The specification of center frequencies is not supported by the FFT-based method')
                end

            elseif strcmp(pObj.filterType,'log')

                % Get center frequencies
                if isempty(pObj.parameters.map('ams_cfHz'))

                    pObj.modCfHz = createFreqAxisLog(pObj.lowFreqHz, ...
                                                     pObj.highFreqHz, ...
                                                 pObj.parameters.map('ams_nFilters'));
                else
                    % Overwrite frequency range
                    pObj.modCfHz = pObj.parameters.map('ams_cfHz');

                    pObj.lowFreqHz  = min(cfHz);
                    pObj.highFreqHz = max(cfHz);
                end

                % Hard-coded filterbank properties
                Q = 1;              % Q-factor
                use_lp = true;      % Use low-pass filter as lowest filter
                use_hp = false;     % Use high-pass for highest filter   

                % Implement second-order butterworth modulation filterbank
                try
                    [pObj.b,pObj.a,pObj.bw] = createFB_Mod(pObj.fs_ds, pObj.modCfHz, ...
                                                           Q, use_lp, use_hp);
                catch
                    % Try/catch is used here to avoid an error when generating a dummy
                    % processor for the hasParameters method.
                end

                % Get bandwidths in hertz
                pObj.bw = pObj.bw*(pObj.fs_ds/2);
                
                % Instantiate the filters
                pObj.Filters = pObj.populateFilters;
                
            end

            % Compute framing parameters
            stepSizeSamples = round(pObj.blockSize / (pObj.blockSec/pObj.stepSec));
            pObj.overlap    = pObj.blockSize - stepSizeSamples;
            pObj.win = window(pObj.winName,pObj.blockSize);

            % Output sampling frequency (input was downsampled, and framed)
            pObj.FsHzOut = pObj.fs_ds/stepSizeSamples;
            
        end
        
    end
    
    % "Getter" methods
    methods
        
        function filterType = get.filterType(pObj)
            filterType = pObj.parameters.map('ams_fbType');
        end
        
        function winName = get.winName(pObj)
            winName = pObj.parameters.map('ams_wname');
        end
        
        function stepSec = get.stepSec(pObj)
            stepSec = pObj.parameters.map('ams_hSizeSec');
        end
        
        function blockSec = get.blockSec(pObj)
            blockSec = pObj.parameters.map('ams_wSizeSec');
        end
        
        function dsRatio = get.dsRatio(pObj)
            dsRatio = pObj.parameters.map('ams_dsRatio');
        end
        
        function nAudioChan = get.nAudioChan(pObj)
            nAudioChan = size(pObj.getDependentParameter('fb_cfHz'),2);
        end
    end
    
    methods (Static)
        
        function dep = getDependency()
            dep = 'innerhaircell';
        end
        
        function [names, defaultValues, descriptions] = getParameterInfo()
            %getParameterInfo   Returns the parameter names, default values
            %                   and descriptions for that processor
            %
            %USAGE:
            %  [names, defaultValues, description] =  ...
            %                           gammatoneProc.getParameterInfo;
            %
            %OUTPUT ARGUMENTS:
            %         names : Parameter names
            % defaultValues : Parameter default values
            %  descriptions : Parameter descriptions
            
            
            names = {'ams_fbType',...
                    'ams_nFilters',...
                    'ams_lowFreqHz',...
                    'ams_highFreqHz',...
                    'ams_cfHz',...
                    'ams_dsRatio',...
                    'ams_wSizeSec',...
                    'ams_hSizeSec',...
                    'ams_wname'};
            
            descriptions = {'Filterbank type (''lin'' or ''log'')',...
                    'Requested number of modulation filters (integer)',...
                    'Lowest modulation center frequency (Hz)',...
                    'Highest modulation center frequency (Hz)',...
                    'Vector of channels'' center frequencies in Hz',...
                    'Downsampling ratio of the envelope',...
                    'Window duration (s)',...
                    'Window step size (s)',...
                    'Window name'};
            
            defaultValues = {'log',...
                            15,...
                            4,...
                            1024,...
                            [],...
                            4,...
                            32E-3,...
                            16E-3,...
                            'rectwin'};
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'Amplitude modulation spectrogram extraction';
            pInfo.label = 'Amplitude modulation';
            pInfo.requestName = 'amsFeatures';
            pInfo.requestLabel = 'Amplitude modulation spectrogram';
            pInfo.outputType = 'ModulationSignal';
            pInfo.isBinaural = false;
            
        end
        
    end
    
    methods (Access = private)
        function obj = populateFilters(pObj)
            % This function is a workaround to assign an array of objects
            % as one of the processor's property, should remain private

            % Total number of filters
            nFilter = numel(pObj.modCfHz)*pObj.nAudioChan;
            
            % Preallocate memory by instantiating last filter
            obj(1,nFilter) = genericFilter(pObj.b{end},pObj.a{end},pObj.fs_ds);
            
            % Instantiating remaining filters
            for ii = 0:pObj.nAudioChan-1
                for jj = 1:numel(pObj.modCfHz)
                    obj(1,ii*numel(pObj.modCfHz)+jj) = genericFilter(pObj.b{jj},pObj.a{jj},pObj.fs_ds);
                end
            end                        
            
        end
    end
    
    
end