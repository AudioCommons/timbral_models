classdef mocProc < Processor
%MOCPROC Medial Olivo-Cochlear (MOC) feedback processor.
%   The MOC feedback processor models the MOC efferent feedback onto the
%   basilar membrane, especially its nonlinear operation. Therefore it only
%   works in conjunction with the DRNL filterbank [1]. This processor takes
%   as the input a time frame-frequency representation from the ratemap
%   processor, and returns the output as the MOC attenuation in dB in the
%   same dimension. This is applied back as the nonlinear path gain of the
%   DRNL filterbank processor. The reflexive feedback is realised by means
%   of a fixed internal relationship between the ratemap output and the MOC
%   attenuation, such as that of Liberman [2], and the reflective feedback 
%   is realised by means of an additional set of parameters to externally 
%   adjust the MOC attenuation.
%   Note that this processor is currently "tuned" to reflect the Clark -
%   Liberman MOC-rate data [1] as closely as possible. The internal
%   parameters and the MOC attenuation derivation method can change if
%   other references are to be used.
%
%   MOCPROC properties:
%       extIpsi     - External Ipsilateral MOC feedback (as nonlinear gain factor)
%       extContra   - External Contralateral MOC feedback (as nonlinear gain factor)
%       thresholdRatedB      
%                   - threshold AN rate (in dB scale) over which the MOC 
%                       attenuation will be driven
%       maxAttenuationdB
%                   - Maximum MOC attenuation in dB
%
%   See also: Processor, drnlProc, ratemapProc
%
%   Reference:
%   [1] Clark, N. R., Brown, G. J., Jurgens, T., & Meddis, R. (2012). 
%    A frequency-selective feedback model of auditory efferent suppression 
%    and its implications for the recognition of speech in noise. 
%    The Journal of the Acoustical Society of America, 132(3), 1535-41. 
%   [2] Liberman, M. C. (1988). Response properties of cochlear efferent 
%    neurons: monaural vs. binaural stimulation and the effects of noise. 
%    Journal of Neurophysiology, 60(5), 1779?98. 

%% Properties
    properties (Dependent = true)

        % The values of these properties are not stored here, but instead read via a
        % "getter" method implemented below
        extIpsi                     % external ipsilateral moc 
        extContra                   % external contralateral moc
        thresholdRatedB          % threshold AN rate (in dB scale) over which the MOC attenuation will be driven
        maxAttenuationdB         % Maximum MOC attenuation in dB
        
    end
    
    properties (GetAccess = private)

        p1                  % used for curve fitting (cubic)
        p2
        p3
        p4
        
    end
            
%% Methods
    methods
        function pObj = mocProc(fs,parObj)
        %mocProc   Construct an MOC processor
        %
        % USAGE:
        %   pObj = mocProc(fs, parObj)
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
            pObj = pObj@Processor(fs,fs,'mocProc',parObj);
            
            if nargin>0 && ~isempty(fs)
                % Add here any additional initialization steps which need only be called
                % once at instantiation. Additional initialization steps that should be
                % performed again after receiving feedback are implemented below
                pObj.p1 = -0.00012411;   % coefficients for MOC-AN rate curve fitting
                pObj.p2 = -0.051322;
                pObj.p3 = -6.3528;
                pObj.p4 = -206.11;
                
            end
        end
        
        function out = processChunk(pObj,in)
            %processChunk       Processing method of MOC processor
            %
            %USAGE
            %       out = processChunk(pObj,in)
            %       out = pObj.processChunk(in)
            %
            %INPUT ARGUMENTS
            %      pObj : processor object
            %        in : rate-map representataion (time frame x frequency)
            %
            %OUTPUT ARGUMENTS
            %       out : MOC attenuation activity in dB
            %
            %SEE ALSO:
            %       templateProc.m
 
                % input level -> ratemap output -> MOC attenuation mapping,
                % based on Clark et al. 2012 paper (and Liberman 1988) data

                % Note that input can be zero - in that case log will give -inf
                % So firstly assume out will be 0 for zero input
                % set up the default zeros for out
                out = zeros(size(in));
            if ~isempty(in)

                % relationship found from curve fitting at 520Hz and 3980Hz
                % (Liberman data)
                % ratemap data should be in log scale first
                x = 20*log10(in);       % this relationship is monotonic
                % The [cubic] curve shape used for the fitting (below) makes it
                % irrelevant to use the "intermediate" x values lower than -180
                % so place a threshold for x here
                x(x<pObj.thresholdRatedB) = pObj.thresholdRatedB;

                out(in>0) = pObj.p1*x(in>0).^3 + pObj.p2*x(in>0).^2 + ...
                    pObj.p3*x(in>0) + pObj.p4;

                % set negative output to zero
                out(out<0) = 0;

                % set output saturation at 40 dB - this also corresponds to the
                % fitted curve shape
                out(out>pObj.maxAttenuationdB) = pObj.maxAttenuationdB;

                % apply the dB value as the "multiplication factor"
                % currently only use the very last (window) frame for the
                % feedback, 
                mocIpsiFeed = pObj.extIpsi .* (10.^(-out(1, :)/20));
                mocContraFeed = pObj.extContra;

                % Access the DRNL MOC parameters directly without
                % affecting anything else
                pObj.LowerDependencies{1}.LowerDependencies{1}.LowerDependencies{1}.parameters.map('fb_mocIpsi') = mocIpsiFeed;
                pObj.LowerDependencies{1}.LowerDependencies{1}.LowerDependencies{1}.parameters.map('fb_mocContra') = mocContraFeed;          
            
            end
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
            
            % Reset the internal states of your processor here, if any. 
            
        end
                
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
        
%         function prepareForProcessing(pObj)
%             
%             % Finalize the initialization of your processor. 
%             
%         end
        

        %% Uncomment and modify only if your processor generates multiple outputs 
        % (including left and right channels of a same type, e.g., pre-processor) or if
        % the signal type it generates needs additional arguments in its constructor
        % (e.g., featureSignal)
        % See, e.g., preProc.m and spectralFeaturesProc.m

%         function output = instantiateOutput(pObj,dObj)
%             %INSTANTIATEOUTPUT  Instantiate the output signal for this processor
%             %
%             % USAGE:
%             %  sObj = pObj.instantiateOutput(dObj)
%             %
%             % INPUT ARGUMENTS:
%             %  pObj : Processor instance
%             %  dObj : Data object instance
%             %
%             % OUTPUT ARGUMENTS:
%             %  sObj : Handle to the instantiated output signal 
%             %
%             % NB: This method should be overridden in children processor where output 
%             % differs from standard (e.g., multiple output or additional inputs to the 
%             % signal constructor).
%             
%             sig = feval(pObj.getProcessorInfo.outputType, ...
%                         pObj, ...
%                         dObj.bufferSize_s, ...
%                         pObj.Channel);
%             
%             dObj.addSignal(sig);
%             
%             output = {sig};
%             
%         end

        %% Uncomment and modify only if your processor is not a single or binaural input 
        % and single output processor.
        % See, e.g., preProc.m

%         function initiateProcessing(pObj)
%             %INITIATEPROCESSING    Wrapper calling the processChunk method and routing I/O
%             % Main purpose is to allow overloading of I/O routing in processors with
%             % "unusual" number of input/outputs.
%             %
%             % Two cases considered here, monaural and binaural processors producing single
%             % outputs. In other cases, the method should be overloaded in the particular
%             % processor.
%             %
%             % USAGE:
%             %   pObj.initiateProcessing
%             %
%             % INPUT ARGUMENT:
%             %   pObj : Processor instance
%             
%             if size(pObj.Input,1) > 1 || numel(pObj.Output) > 1
%                 % Then it is a multiple-input processor, return an error
%                 error(['Cannot initiate the processing for this ' ...
%                     'processor. Consider overloading this method in the children ' ... 
%                     'class definition.'])
%             
%             elseif size(pObj.Input,2) == 1
%                 % Then monaural processor
%                 pObj.Output{1}.appendChunk( ...
%                     pObj.processChunk( pObj.Input{1}.Data('new')));
%                 
%             elseif size(pObj.Input,2) == 2
%                 % Then binaural processor
%                 pObj.Output{1}.appendChunk( ...
%                     pObj.processChunk( ...
%                         pObj.Input{1}.Data('new'), pObj.Input{2}.Data('new')));
%                 
%                 
%             else
%                 % TODO: Remove after testing
%                 error('Something is wrong with the inputs of this processor, investigate.')
%             end
%             
%         end


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


        %% Uncomment and implement only if your processor has multiple inputs of different
        % types (i.e., not just 'left' and 'right' channel of a same type)
        % NOTE: There are no such processor in the AFE at the moment, so you're on your
        % own!

%         function addInput(pObj,dependency)
%             %ADDINPUT   Adds a signal to the list of input of the processor
%             % Will consider that the input signal is the output of the dependency. If
%             % the dependency has a left- and right-channel output, will pick the suitable
%             % one.
%             % Three cases are implemented here:
%             % 1- single dependency with single output
%             % 2- two dependencies (left and right channels) with each single output
%             % 3- single dependency with left and right outputs (e.g., pre-processor)
%             %
%             % Should the input attribution works in any other way for a given processor,
%             % this method should be overloaded for that specific children processor.
%             %
%             
%             
%         end


        %% Uncomment and implement only if your processors returns multiple outputs of
        % different types (i.e., not just 'left' and 'right' channels of a same type)

%         function addOutput(pObj,sObj)
%             %ADDOUTPUT  Adds a signal to the list of output of the processor
%             %
%             % USAGE:
%             %   pObj.addOutput(sObj)
%             %
%             % INPUT ARGUMENTS:
%             %   pObj : Processor instance
%             %   sObj : Single instance or cell array of left and right channel instances
%             %          of signal objects
%             
%             
%         end

    end
    
    %% "Getter" methods
    methods
        % All properties listed as "Dependent" in the beginning of the class definition
        % should have a corresponding "getter" method here.
        
        function extIpsi = get.extIpsi(pObj)
            extIpsi = pObj.parameters.map('moc_extIpsi');
        end
        
        function extContra = get.extContra(pObj)
            extContra = pObj.parameters.map('moc_extContra'); 
        end
        
        function thresholdRatedB = get.thresholdRatedB(pObj)
            thresholdRatedB = pObj.parameters.map('moc_thresholdRatedB'); 
        end

        function maxAttenuationdB = get.maxAttenuationdB(pObj)
            maxAttenuationdB = pObj.parameters.map('moc_maxAttenuationdB'); 
        end
        
        
    end
    
    %% Static methods
    % These methods store all hard-coded information regarding your processor, remember to
    % update them with actual information!
    methods (Static)
        
        function dep = getDependency()
            dep = 'ratemap';
        end
        
        function [names, defaultValues, descriptions] = getParameterInfo()
            %getParameterInfo   Returns the parameter names, default values
            %                   and descriptions for that processor
            %
            %USAGE:
            %  [names, defaultValues, description] =  ...
            %                           mocProc.getParameterInfo;
            %
            %OUTPUT ARGUMENTS:
            %         names : Parameter names
            % defaultValues : Parameter default values
            %  descriptions : Parameter descriptions
            
            
            names = {'moc_extIpsi', ...
                'moc_extContra', ...
                'moc_thresholdRatedB', ...
                'moc_maxAttenuationdB'};
            
            descriptions = {'External ipsilateral MOC feedback factor', ...
                            'External contralateral MOC feedback factor', ...
                            'Threshold ratemap value for MOC activation in dB', ...
                            'Maximum possible MOC attenuation in dB'};
            
            defaultValues = {1, 1, -180, 40};
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'MOC processor';
            pInfo.label = 'Medial Olivo-Cochlear feedback processor';
            pInfo.requestName = 'moc';
            pInfo.requestLabel = 'Medial Olivo-Cochlear feedback';
            pInfo.outputType = 'TimeFrequencySignal';
            pInfo.isBinaural = false; % 0 or 1 (or 2 if it can do both, e.g., preProc.m)
            
        end
        
    end
        
end