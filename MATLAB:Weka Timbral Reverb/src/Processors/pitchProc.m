classdef pitchProc < Processor
%PITCHPROC Pitch processor.
%   Based on the Summary Auto-Correlation Function representation, this processor
%   detects the most salient peak within the given plausible pitch frequency range 
%   for each time frame in order to obtain an estimation of the fundamental 
%   frequency.
%
%   PITCHPROC properties:
%        pitchRangeHz    - Range in Hz for valid pitch estimation
%        confThresPerc   - Threshold for pitch condidence measure (re. 1)
%        orderMedFilt    - Median order filter for pitch smoothing
%        lags            - Vector of auto-correlation lags
%
%   See also: Processor, autocorrelationProc

    properties (Dependent = true)
        pitchRangeHz    % Range in Hz for valid pitch estimation
        confThresPerc   % Threshold for pitch condidence measure (re. 1)
        orderMedFilt    % Median order filter for pitch smoothing
    end
    
    properties (SetAccess = protected)
        lags            % Vector of auto-correlation lags
    end
    
    properties (Access = protected)
        bValidLags      % Which lags are in the pitch range
%         pitchBuffer     % Buffer for online Median filtering
        maxConf         % Buffer interface for maximum confidence
        maxConfBuf      % Circular buffer for maximum confidence
        
    end
    
    
    methods
        function pObj = pitchProc(fs,parObj)
		%pitchProc   Construct a pitch extraction processor
        %
        % USAGE:
        %   pObj = pitchProc(fs, parObj)
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
            
            % Checking input parameter
            if nargin<2||isempty(parObj); parObj = Parameters; end
            if nargin<1; fs = []; end
            
            % Call superconstructor
            pObj = pObj@Processor(fs,fs,'pitchProc',parObj);
            
            
        end
        
        function out = processChunk(pObj,in)
            %processChunk       Apply the processor to a new chunk of input
            %                   signal
            %
            %USAGE
            %   out = pObj.processChunk(in)
            %
            %INPUT ARGUMENT
            %    in : New chunk of input data
            %
            %OUTPUT ARGUMENT
            %   out : Corresponding output
            
            
            % Compute summary ACF
            sacf = squeeze(mean(in,2));
            
            % Input size
            [nFrames,nLags] = size(sacf);
            
            % Restrict lags to plausible pitch range (only for first call)
            if isempty(pObj.bValidLags)
                rangeLagSec = 1./pObj.pitchRangeHz;
                pObj.bValidLags = (pObj.lags >= min(rangeLagSec)) & ...
                    (pObj.lags <= min(max(rangeLagSec),nLags));
            end
            
            % Restrict lags to predefined pitch range
            sacf = sacf(:,pObj.bValidLags);
            lagsSec = pObj.lags(pObj.bValidLags);
            
            %% DETECT PITCH CANDIDATES
            % 
            % 
            % Allocate memory
            pitchHzRaw = zeros(nFrames,1);
            confidence = zeros(nFrames,1);

            % Loop over number of frames
            for ii = 1 : nFrames

                % Detect local peaks
                [peakIdx,peakVal] = findpeaks_VB(sacf(ii,:));

                % Find maximum peak position and confidence value
                [maxVal,maxIdx] = max(peakVal);

                % Confidence value
                if isempty(maxVal)
                    maxVal = 0;
                else
                    confidence(ii) = maxVal;
                    pObj.maxConfBuf.append(max(maxVal,0));
                end

                % Only accept pitch estimate if confidence value is above 0
                if maxVal > 0
                    % Pitch estimate in Hertz
                    pitchHzRaw(ii) = 1/lagsSec(peakIdx(maxIdx));
                end
            end
            
            %% POST-PROCESSING
            % 
            % 
            % Floor confidence value
            confidence = max(confidence,0);

            % Compute threshold
%             pObj.confThresPerc = max(confidence,[],1) * pObj.confThresPerc;
            confThresPerc = max(pObj.maxConf(:)) * pObj.confThresPerc;

            % Apply confidence threshold
            bSetToZero = confidence < repmat(confThresPerc,[nFrames 1]);

            % Set unreliable pitch estimates to zero
            pitchHz = pitchHzRaw; pitchHz(bSetToZero) = 0;
            
            % Apply median filtering to reduce octave errors
            Npre = floor(pObj.orderMedFilt/2);      % Past samples considered in median
            Npost = ceil(pObj.orderMedFilt/2)-1;    % Future samples considered in median
            
            % We have to fiddle a bit with online compatibility...
            % For a filter of order N,
            %   - Last (N-1)/2 (N odd) or N/2-1 (N even) sample cannot be obtained
            %   - Corresponding ones from the previous chunk can now be obtained
            
            % Append buffer
%             pitchHzBuf = [pObj.pitchBuffer ; pitchHz];
            
%             pitchHzFilt = medfilt1(pitchHzBuf,pObj.orderMedFilt);
            filtPitch = medfilt1(pitchHz,pObj.orderMedFilt);
            
            % Discard the last Npost samples (accurate median estimation would need next 
            % input chunk) and the first Npre samples (used only to compute
            % Npre+1:Npre+Npost samples)
%             filtPitch = pitchHzFilt(Npre+1:end-Npost);

            % NB: The above statement leads to a pitch estimate which is shifted in time
            % by Npost samples. Could we correct for that?
            
            % Update buffer: the last Npost samples that were discarded need to be
            % computed at next chunk, so we need Npost+Npre extra samples in the buffer
%             pObj.pitchBuffer = pitchHzBuf(end-Npost-Npre:end);
            
            % Replace all zeros with NANs
            filtPitch(filtPitch==0) = NaN;
            
            % TODO: Solve the following
%             filtPitch = [NaN;filtPitch;NaN];
            
            % Generate the output
            out = [filtPitch pitchHzRaw confidence];
%             out = filtPitch;
%             out = pitchHzFilt;
%             out = bSetToZero;
%             out = pitchHzRaw;
%             out = confidence;
            
        end
        
        function reset(~)
%             pObj.pitchBuffer = [];
        end
        
        function output = instantiateOutput(pObj,dObj)
            %INSTANTIATEOUTPUT  Instantiate the output signal for this processor
            %
            %NB: This method is overloaded here from the master Processor class, as
            %feature signals need additional input argument to construct
            
            featureNames = {'pitch','rawPitch','confidence'};
            
            sig = feval(pObj.getProcessorInfo.outputType, ...
                        pObj, ...
                        dObj.bufferSize_s, ...
                        pObj.Channel,...
                        [],...
                        featureNames);
            
            dObj.addSignal(sig);
            
            output = {sig};
            
        end
        
    end
    
    methods (Access=protected)
        
        function verifyParameters(pObj)
            
            % The median filter order should be an integer
            if mod(pObj.orderMedFilt,1)~=0
                pObj.parameters.map('pi_medianOrder') = ...
                        round(pObj.parameters.map('pi_medianOrder'));
                warning('Median filter order should be an integer, using %i instead',...
                    pObj.parameters.map('pi_medianOrder'))
            end
            
            
        end
        
    end
    
    methods (Hidden = true)
        
        function prepareForProcessing(pObj)
            
            % Compute internal parameters
            bufferDurSec = 5;   % Maximum confidence is taken in the past 5 seconds
            pObj.maxConfBuf = circVBuf(bufferDurSec*pObj.FsHzIn,1);
            pObj.maxConf = circVBufArrayInterface(pObj.maxConfBuf);
            pObj.lags = pObj.getDependentProperty('lags');
            
            % Reset the valid lags vector
            pObj.bValidLags = [];
            
        end
        
    end
    
    % "Getter" methods
    methods
        function pitchRangeHz = get.pitchRangeHz(pObj)
            pitchRangeHz = pObj.parameters.map('pi_rangeHz');
        end
        
        function confThresPerc = get.confThresPerc(pObj)
            confThresPerc = pObj.parameters.map('pi_confThres');
        end
        
        function orderMedFilt = get.orderMedFilt(pObj)
            orderMedFilt = pObj.parameters.map('pi_medianOrder');
        end
    end
    
    methods (Static)
        
        function dep = getDependency()
            dep = 'autocorrelation';
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
            
            
            names = {'pi_rangeHz',...
                    'pi_confThres',...
                    'pi_medianOrder'};
            
            descriptions = {'Range in Hz for valid pitch estimation',...
                    'Threshold for pitch condidence measure (re. 1)',...
                    'Median order filter for pitch smoothing (integer)'};
            
            defaultValues = {[80 400],...
                            0.7,...
                            3};
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'Pitch';
            pInfo.label = 'Pitch';
            pInfo.requestName = 'pitch';
            pInfo.requestLabel = 'Pitch estimation';
            pInfo.outputType = 'FeatureSignal';
            pInfo.isBinaural = false;
            
        end
        
    end
    
    
end