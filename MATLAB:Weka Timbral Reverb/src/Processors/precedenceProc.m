classdef precedenceProc < Processor
%PRECEDENCEPROC Precedence effect processor.
%   This processor estimates the interaural time and level differences (ITD
%   and ILD) of a binaural input signal in the presence of a reflection, 
%   by mimicking the precedence effect. The implementation is based on the
%   model developed by Braasch [1]. The input signal is expressed in 
%   time-frequency domain, corresponding to the basilar membrane output,
%   and the output ITD and ILD are calculated in time domain through 
%   weighted summation across frequency channels. 
%
%   PRECEDENCEPROC properties:
%       wSizeSec        - Window duration in seconds
%       hSizeSec        - Window step size
%       maxDelaySec     - Maximum delay in correlation computation in seconds
%
%   Note that this processor needs the following functions stored in the 
%   /src/Tools/ folder:
%       prec_acmod.m
%       prec_anaone.m
%       prec_de_conv.m
%       prec_peakratio.m
%       prec_reconHW.m
%       prec_rev_conv.m
%       calcXCorr.m 
%
%   See also: Processor, gammatoneProc
%
%   Reference:
%   [1] Jonas Braasch (2013), A precedence effect model to simulate 
%       localization dominance using an adaptive, stimulus parameter-based 
%       inhibition process, J. Acoust. Soc. Am., Vol. 134, No. 1

%% Properties

    properties (Dependent = true)
        % external properties of processor, i.e., those associated with
        % a user-controllable parameter
        % The values of these properties are not stored here, but instead read via a
        % "getter" method implemented below
        wSizeSec            % windowlength for temporal steps
        hSizeSec            % Window step size 
        lags                % Vector of lags at which cross-correlation is computed
        maxDelaySec         % Maximum delay in correlation computation (s)      
    end
    
    properties (GetAccess = private)
        % internal properties of processor that do not need public access
        wSize
        hSize
        win         % Window vector
        buffer_l    % Buffered input signals (left ear)
        buffer_r    % Buffered input signals (right ear)
        
        stateStore     % data stored from previous calculation
                       % includes cumulative ac, cc, ITD, ILD, and frame
                       % counter

    end

    
%% Methods
    methods
        function pObj = precedenceProc(fs,parObj)
        %precedenceProc   Construct a precedence effect processor
        %
        % USAGE:
        %   pObj = precedenceProc(fs, parObj)
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
            
            if nargin<2||isempty(parObj); parObj = Parameters; end
            if nargin<1; fs = []; end
            
            % Call super-constructor
            pObj = pObj@Processor(fs,fs,'precedenceProc',parObj);
            
            if nargin>0 && ~isempty(fs)
                % additional initialization steps which need only be called
                % once at instantiation. Additional initialization steps that should be
                % performed again after receiving feedback are implemented below
                                
                % initialise stateStore 
                pObj.stateStore = struct('ac1', [], 'ac2', [], ...
                    'cc1', [], 'cc2', [], 'cc3', [], 'cc4', [], ...
                    'ICCintT', [], 'ILDint', [], 'Eint', [], 'lastFrame', []);

            end
        end
        
        function [itd, ild, cc] = processChunk(pObj,in_l,in_r)
            %processChunk       Processing method of your processor
            %
            %USAGE
            %       out = processChunk(pObj,in)
            %       out = pObj.processChunk(in)
            %
            %INPUT ARGUMENTS
            %      pObj : precedenceProc processor object
            %      in_l,
            %      in_r : time-frequency signal (time (row) x frequency (column))
            %
            %OUTPUT ARGUMENTS
            %       itd : ITD output in given time frame
            %             (timeframe) x (frequency)
            %       ild : ILD output in given time frame
            %             (timeframe) x (frequency)
            %       cc  : CC output in given time frame
            %             (timeframe) x (frequency) x (lags)
            %
            %SEE ALSO:
            %       precedenceProc.m
                      
            % Append provided input to the buffer
            if ~isempty(pObj.buffer_l)
                in_l = [pObj.buffer_l;in_l];
                in_r = [pObj.buffer_r;in_r];
            end
            
            % Quick control of dimensionality
            if max(size(in_l)~=size(in_r))
                error('Buffered inputs should be of same dimension for both ears')
            end
            
            [nSamples,nChannels] = size(in_l);
            
            % How many frames are in the buffered input?
            nFrames = floor((nSamples-(pObj.wSize-pObj.hSize))/pObj.hSize);
            
            % Determine maximum lag
            maxLag = ceil(pObj.maxDelaySec*pObj.FsHzIn);
            
            ITD = zeros(nChannels, nFrames);    % (freq) x (time frame)
            ILD = zeros(nChannels, nFrames);    % because of acmod output format
            % Allocate memory also for cc output
            CC = zeros(max(0,nFrames),nChannels,maxLag*2+1);
            
            for ii = 1:nFrames
                % Get start and end indices for the current frame
                n_start = (ii-1)*pObj.hSize+1;
                n_end = (ii-1)*pObj.hSize+pObj.wSize;
               
                % Extract current frame (this will be in time-frequency domain)
                frame_l = repmat(hanning(pObj.wSize), 1, nChannels) .* in_l(n_start:n_end, :);
                frame_r = repmat(hanning(pObj.wSize), 1, nChannels) .* in_r(n_start:n_end, :);

                [CC(ii, :, :),ITD(:, ii),ILD(:, ii),lagL,lagR,BL,BR,LLARmode,pObj.stateStore] = ...
                    prec_acmod(frame_l, frame_r, pObj.FsHzIn, ...
                    pObj.LowerDependencies{1}.cfHz, maxLag, pObj.stateStore);   
                
            end
                       
            itd = ITD.';          % prec_acmod function used to return ITD in ms!!
            ild = ILD.';          % output itd, ild are (time)x(freq) signals
            cc = CC;
            
            % Update the buffer: the input that was not extracted as a
            % frame should be stored            
            pObj.buffer_l = in_l(nFrames*pObj.hSize+1:end,:);
            pObj.buffer_r = in_r(nFrames*pObj.hSize+1:end,:);
        end
        
        function reset(pObj)
            %reset          Order the processor to reset its internal
            %               states, e.g., when some critical parameters in
            %               the processing have been changed
            %USAGE
            %       pObj.reset()
            %       reset(pObj)
            %
            %INPUT ARGUMENT
            %       pObj : Processor object
                        
            % If any of the stateStore structure is not empty clean the
            % whole structure
            if any(~structfun(@isempty, pObj.stateStore))
                pObj.stateStore = struct('ac1', [], 'ac2', [], ...
                    'cc1', [], 'cc2', [], 'cc3', [], 'cc4', [], ...
                    'ICCintT', [], 'ILDint', [], 'Eint', [], 'lastFrame', []);
            end
            
        end
        
        % Uncomment and implement the following method only if your processor is one of
        % many alternative for a given processing step
        % See, e.g., gammatoneProc.m and drnlProc.m
        
%         function bInBranch = isSuitableForRequest(pObj)
%             % Look into pObj.parameters to determine if this processor is the right one
%             % considering the parameters entered by the user, return 1 if it is, or 0 if
%             % it is not.
%         end
        
    end
    
    methods (Access=protected)
        
        % Uncomment and implement the following method only if you should perform a check
        % on user-provided parameters or solve potential conflicts between parameters.
        % See, e.g., gammatoneProc.m or ihcProc.m
        
%         function verifyParameters(pObj)
%                         
%             % Solve any potential conflicts and/or control the validity of parameters
%             % here.
%             % Invalid parameter values are usually replaced with their corresponding
%             % default value.
%             
%         end
        
    end
    
    % "Overridden" methods
    methods (Hidden = true)
        % NOTE: All the methods in this block are for processors with "special" behavior.
        % Examples were provided when available.
        % The methods below are approximately ordered from the most likely to the least 
        % likely to be overridden.
        
        
        %% Uncomment and implement only if some internal parameters or properties of your
        % processor are subject to change following user feedback
        % All initialization steps that are potentially part of a "re-initialization"
        % following a change of parameter should go in here. If an internal parameter as a
        % value dependent on other processors this one depends on, it should be
        % initialized here, using e.g., the getDependentProperty or getDependentParameter
        % method of your processor.
        % See, e.g., ratemapProc.m
        
        function prepareForProcessing(pObj)
            
            % Finalize initialization of processor 
            pObj.wSize = 2*round(pObj.parameters.map('prec_wSizeSec')*pObj.FsHzIn/2);
            pObj.hSize = round(pObj.parameters.map('prec_hSizeSec')*pObj.FsHzIn);
            % Output sampling frequency
            pObj.FsHzOut = 1/(pObj.hSizeSec);
        end
        

        %% Uncomment and modify only if your processor generates multiple outputs 
        % (including left and right channels of a same type, e.g., pre-processor) or if
        % the signal type it generates needs additional arguments in its constructor
        % (e.g., featureSignal)
        % See, e.g., preProc.m and spectralFeaturesProc.m

        function output = instantiateOutput(pObj,dObj)
            %INSTANTIATEOUTPUT  Instantiate the output signal for this processor
            %
            % USAGE:
            %  sObj = pObj.instantiateOutput(dObj)
            %
            % INPUT ARGUMENTS:
            %  pObj : Processor instance
            %  dObj : Data object instance
            %
            % OUTPUT ARGUMENTS:
            %  sObj : Handle to the instantiated output signal 
            %
            % NB: This method should be overridden in children processor where output 
            % differs from standard (e.g., multiple output or additional inputs to the 
            % signal constructor).
            
            % Specify two outputs (itd and ild) for the received dObj
            sig_itd = feval(pObj.getProcessorInfo.outputType{1}, ...
                        pObj, ...
                        dObj.bufferSize_s, ...
                        'mono', ...
                        []);
            sig_ild = feval(pObj.getProcessorInfo.outputType{1}, ...
                        pObj, ...
                        dObj.bufferSize_s, ...
                        'mono', ...
                        []);            
            % adding CorrelationSignal to have cc output as well
            sig_cc = feval(pObj.getProcessorInfo.outputType{2}, ...
                        pObj, ...
                        dObj.bufferSize_s, ...
                        'mono', ...
                        []);        
            
            dObj.addSignal(sig_itd);
            dObj.addSignal(sig_ild);
            dObj.addSignal(sig_cc);
            
            output = {sig_itd sig_ild sig_cc};
            
        end

        %% Uncomment and modify only if your processor is not a single or binaural input 
        % and single output processor.
        % See, e.g., preProc.m

        function initiateProcessing(pObj)
            %INITIATEPROCESSING    Wrapper calling the processChunk method and routing I/O
            % Main purpose is to allow overloading of I/O routing in processors with
            % "unusual" number of input/outputs.
            %
            % Two cases considered here, monaural and binaural processors producing single
            % outputs. In other cases, the method should be overloaded in the particular
            % processor.
            %
            % USAGE:
            %   pObj.initiateProcessing
            %
            % INPUT ARGUMENT:
            %   pObj : Processor instance
            
            % Call processChunk function using the binaural inputs
            [itd, ild, cc] = pObj.processChunk( pObj.Input{1,1}.Data('new'),...
                pObj.Input{1,2}.Data('new'));
            % Append the [three] outputs
            pObj.Output{1}.appendChunk(itd);
            pObj.Output{2}.appendChunk(ild);
            pObj.Output{3}.appendChunk(cc);

        end


        %% Uncomment and modify only if your processor depends on multiple type of other
        % processors
        % NOTE: There are no such processor in the AFE at the moment, so you're on your
        % own!
        
%         function addUpperDependencies(pObj,dependentProcs)
%             %ADDUPPERDEPENDENCIES   Populate link to higher processors relying on this one
%             %
%             % USAGE:
%             %  pObj.addUpperDependencies(dependentProcs)
%             %
%             % INPUT ARGUMENTS:
%             %           pObj : Processor instance
%             % dependentProcs : Cell array of (correctly ordered) processors dependent on
%             %                  this processor
%             
%             pObj.UpperDependencies = [pObj.UpperDependencies dependentProcs];
%             
%         end


        %% Uncomment and modify only if your processor depends on multiple type of other
        % processors
        % NOTE: There are no such processor in the AFE at the moment, so you're on your
        % own!

%         function removeUpperDependency(pObj,upperProc)
%             %REMOVEUPPERDEPENDENCY   Removes a processor from the upper dependency list
%             %
%             % USAGE:
%             %    pObj.removeUpperDependency(upperProc)
%             %
%             % INPUT ARGUMENTS:
%             %       pObj : Processor instance
%             %  upperProc : Processor depending on pObj to remove from the dependency list
%             
%             for ii = 1:size(pObj.UpperDependencies,2)
%                 if pObj.UpperDependencies{ii} == upperProc
%                     pObj.UpperDependencies{ii} = [];
%                 end
%             end
%             
%             % Clean up empty elements
%             pObj.UpperDependencies = ...
%                 pObj.UpperDependencies( ~cellfun( @isempty, pObj.UpperDependencies ));
%             
%         end
%         


        %% Uncomment and implement only if your processors returns multiple outputs of
        % different types (i.e., not just 'left' and 'right' channels of a same type)

        function addOutput(pObj,sObj)
            %ADDOUTPUT  Adds a signal to the list of output of the processor
            %
            % USAGE:
            %   pObj.addOutput(sObj)
            %
            % INPUT ARGUMENTS:
            %   pObj : Processor instance
            %   sObj : Single instance or cell array of left and right channel instances
            %          of signal objects
            
            % Attach the outputs to the pObj
            for ii = 1:numel(sObj)
                pObj.Output{1, ii} = sObj{ii};
            end           
        end

    end
    
    %% "Getter" methods
    methods
        % All properties listed as "Dependent" in the beginning of the class definition
        % should have a corresponding "getter" method here.
        
        function wSizeSec = get.wSizeSec(pObj)
            wSizeSec = pObj.parameters.map('prec_wSizeSec');
        end
        
        function hSizeSec = get.hSizeSec(pObj)
            hSizeSec = pObj.parameters.map('prec_hSizeSec');
        end
                
        function lags = get.lags(pObj)
            maxLag = ceil(pObj.maxDelaySec*pObj.FsHzIn);
            lags = (-maxLag:maxLag).'./pObj.FsHzIn;
        end
        
        function maxDelaySec = get.maxDelaySec(pObj)
            maxDelaySec = pObj.parameters.map('prec_maxDelaySec');
        end
               
    end
    
    %% Static methods
    % These methods store all hard-coded information regarding the processor
    
    methods (Static)
        
        function dep = getDependency()
            dep = 'filterbank';
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
            
            
            names = {'prec_wSizeSec',...
                     'prec_hSizeSec',...
                     'prec_maxDelaySec'
                    };
            
            descriptions = {'window length for temporal steps (seconds)',...
                    'Window step size (s)',...
                    'Maximum delay in correlation computation (s)'
                    };
            
            defaultValues = {20e-3,...
                             10e-3,...
                             1e-3
                            };
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'Precedence effect';
            pInfo.label = 'Precedence effect';
            pInfo.requestName = 'precedence';
            pInfo.requestLabel = 'Precedence effect';
            pInfo.outputType = {'TimeFrequencySignal' 'CorrelationSignal'};
            pInfo.isBinaural = 1; % 0 or 1 (or 2 if it can do both, e.g., preProc.m)
            
        end
        
    end
        
end