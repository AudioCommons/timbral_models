classdef preProc < Processor
%PREPROC Pre-processor.
%   Prior to computing any of the supported auditory representations, 
%   the input signal stored in the data object can be pre-processed with 
%   one of the following elements:
%       1. Direct current (DC) bias removal
%       2. Pre-emphasis
%       3. Root mean square (RMS) normalization [1] 
%       4. Level scaling to a pre-defiend sound pressure level (SPL) reference
%       5. Middle ear filtering [2] 
%
%   PREPROC properties:
%       bRemoveDC       - Flag to activate DC removal filter
%       cutoffHzDC      - Cutoff frequency in Hz of the high-pass filter
%         
%       bPreEmphasis    - Flag to activate pre-emphasis filter
%       coefPreEmphasis - Coefficient of first-order high-pass filter
%         
%       bNormalizeRMS   - Flag to activate RMS normalization
%       bBinauralRMS    - Flag to link RMS normalization across both ears
%       intTimeSecRMS   - Time constant used for RMS estimation
%         
%       bLevelScaling   - Flag to apply level sacling to given reference
%       refSPLdB        - Reference dB SPL to correspond to input RMS of 1
%         
%       bMiddleEarFiltering - Flag to apply middle ear filtering
%       middleEarModel      - Middle ear filter model
%
%   See also: Processor
%
%   Reference:
%   [1] Tchorz, J. and Kollmeier, B. (2003), "SNR estimation based on 
%       amplitude modulation analysis with applications to noise suppression," 
%       IEEE Transactions on Audio, Speech, and Language Processing 11(3),
%       pp. 184?192.
%   [2] Goode, R. L., Killion, M., Nakamura, K., and Nishihara, S. (1994), 
%       "New knowledge about the function of the human middle ear: 
%       development of an improved analog model." The American journal of 
%       otology 15(2), pp. 145?154.
    
    properties (SetAccess = protected, Dependent = true)
        bRemoveDC
        cutoffHzDC
        
        bPreEmphasis
        coefPreEmphasis
        
        bNormalizeRMS
        bBinauralRMS
        intTimeSecRMS
        
        bLevelScaling
        refSPLdB
        
        bMiddleEarFiltering
        middleEarModel
        
        bUnityComp
        
    end
    
    properties (Access = private)
        dcFilter_l
        dcFilter_r
        preEmphFilter_l
        preEmphFilter_r
        agcFilter_l
        agcFilter_r
        epsilon = 1E-8;
        midEarFilter_l
        midEarFilter_r
        meFilterPeakdB
    end
    
    
    methods
        function pObj = preProc(fs,parObj)
		%preProc   Construct a pre-processor
        %
        % USAGE:
        %   pObj = preProc(fs, parObj)
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
            pObj = pObj@Processor(fs,fs,'preProc',parObj);
            
            % This processor can take two inputs and two outputs
            pObj.isBinaural = true;
            pObj.hasTwoOutputs = true;
            
        end
        
        function [out_l, out_r] = processChunk(pObj,in_l,in_r)
            %processChunk       Apply the processor to a new chunk of input signal
            %
            %USAGE
            %   [out_l out_r] = pObj.processChunk(in_l,in_r)
            %
            %INPUT ARGUMENT
            %    in_l : New chunk of input data from left channel
            %    in_r : New chunk of input data from right channel
            %
            %OUTPUT ARGUMENT
            %   out_l : Corresponding output (left channel)
            %   out_r : Corresponding output (right channel)
            
            if nargin < 3 || isempty(in_r)
                in_r = [];
            end
            
            % 0- Initialization
            data_l = in_l;
            data_r = in_r;
            
            % 1- DC-removal filter
            if pObj.bRemoveDC
                data_l = pObj.dcFilter_l.filter(data_l);
                data_r = pObj.dcFilter_r.filter(data_r);
            end
            
            % 2- Pre-whitening
            if pObj.bPreEmphasis
                data_l = pObj.preEmphFilter_l.filter(data_l);
                data_r = pObj.preEmphFilter_r.filter(data_r);
            end
            
            % 3- Automatic gain control
            if pObj.bNormalizeRMS
                % Initialize the filter states if empty
                if ~pObj.agcFilter_l.isInitialized  %isempty(pObj.agcFilter_l.States)
                    % Mean square of input over the time constant
                    sm_l = mean(data_l(1:min(size(data_l,1),round(pObj.intTimeSecRMS*pObj.FsHzIn))).^2);
                    sm_r = mean(data_r(1:min(size(data_r,1),round(pObj.intTimeSecRMS*pObj.FsHzIn))).^2);
                    
                    % Initial filter states
                    s0_l = exp(-1/(pObj.intTimeSecRMS*pObj.FsHzIn))*sm_l;
                    s0_r = exp(-1/(pObj.intTimeSecRMS*pObj.FsHzIn))*sm_r;
                    
                    pObj.agcFilter_l.reset(s0_l)
                    pObj.agcFilter_r.reset(s0_r)
                end
                
                % Estimate normalization constants
                normFactor_l = sqrt(pObj.agcFilter_l.filter(data_l.^2))+pObj.epsilon;
                normFactor_r = sqrt(pObj.agcFilter_r.filter(data_r.^2))+pObj.epsilon;
                
                % Preserve multi-channel differences
                if ~isempty(normFactor_r) && pObj.bBinauralRMS
                    normFactor_l = max(normFactor_l,normFactor_r);
                    normFactor_r = normFactor_l;
                end
                
                % Apply normalization
                data_l = data_l./normFactor_l;
                if ~isempty(normFactor_r)
                    data_r = data_r./normFactor_r;
                else
                    data_r = [];
                end
                
            end
            
            if pObj.bLevelScaling
                current_dboffset = dbspl(1);
                if isscalar(pObj.refSPLdB)
                    data_l = gaindb(data_l, current_dboffset-pObj.refSPLdB);
                    data_r = gaindb(data_r, current_dboffset-pObj.refSPLdB);
                else
                    data_l = gaindb(data_l, current_dboffset-pObj.refSPLdB(1));
                    data_r = gaindb(data_r, current_dboffset-pObj.refSPLdB(2));
                end
            end
            
            if pObj.bUnityComp
                switch pObj.middleEarModel
                    case 'jepsen'
                        pObj.meFilterPeakdB = 55.9986;
                    case 'lopezpoveda'
                        pObj.meFilterPeakdB = 66.2888;
                end
            else
                pObj.meFilterPeakdB = 0;
            end

            if pObj.bMiddleEarFiltering
                data_l = pObj.midEarFilter_l.filter(data_l)* 10^(pObj.meFilterPeakdB/20);
                data_r = pObj.midEarFilter_r.filter(data_r)* 10^(pObj.meFilterPeakdB/20);
                
            end
            
            % Return the output
            out_l = data_l;
            out_r = data_r;
            
        end
          
        function reset(pObj)
            %reset     Resets the internal states of the pre-processor
            %
            %USAGE
            %      pObj.reset
            %
            %INPUT ARGUMENTS
            %  pObj : Pre-processor instance
            
            if pObj.bRemoveDC
                pObj.dcFilter_l.reset;
                pObj.dcFilter_r.reset;
            end
            if pObj.bPreEmphasis
                pObj.preEmphFilter_l.reset;
                pObj.preEmphFilter_r.reset;
            end
            if pObj.bNormalizeRMS
                pObj.agcFilter_l.reset;
                pObj.agcFilter_r.reset;
            end
            if pObj.bMiddleEarFiltering
                pObj.midEarFilter_l.reset;
                pObj.midEarFilter_r.reset;
            end
            
        end
    
    end
    
    methods (Access=protected)
        
        function verifyParameters(pObj)
            
            % TODO: Add more? e.g., what follows
%             if numel(p.pp_refSPLdB)>2
%                 fprintf('More than two refSPLdB values given - only the first two will be used (L/R)');
%             end
%             if pObj.bUnityComp
%                 switch pObj.middleEarModel
%                     case 'jepsen'
%                         pObj.meFilterPeakdB = 55.9986;
%                     case 'lopezpoveda'
%                         pObj.meFilterPeakdB = 66.2888;
%                 end
%             else
%                 pObj.meFilterPeakdB = 0;
%             end
        end
        
    end
    
    % Pre-processor is a multi-channel processor and needs to overload some of the
    % standard Processor methods to function correctly
    methods (Hidden = true)
        function output = instantiateOutput(pObj,dObj)
            %INSTANTIATEOUTPUT  Instantiate the output signal for the pre-processor
            %
            %NB: This method is overloading the parent method to deal with multiple
            %outputs
            
            if dObj.isStereo
                sig_l = feval(pObj.getProcessorInfo.outputType, ...
                            pObj, ...
                            dObj.bufferSize_s, ...
                            'left');
                sig_r = feval(pObj.getProcessorInfo.outputType, ...
                            pObj, ...
                            dObj.bufferSize_s, ...
                            'right');
                dObj.addSignal(sig_l);
                dObj.addSignal(sig_r);
            
                output = {sig_l, sig_r};
                
            else
                sig = feval(pObj.getProcessorInfo.outputType, ...
                            pObj, ...
                            dObj.bufferSize_s, ...
                            pObj.Channel);
            
                dObj.addSignal(sig);
            
                output = {sig};
                
            end
            
        end
        
        function initiateProcessing(pObj)
            %INITIATEPROCESSING    Wrapper calling the processChunk method and routing I/O
            % Because the pre-processor can have two outputs (for stereo signals), it is
            % necessary to overload the parent method here.
            
            
            if size(pObj.Input,2)>1
                [out_l, out_r] = pObj.processChunk( pObj.Input{1,1}.Data('new'),...
                    pObj.Input{1,2}.Data('new'));
            else
                [out_l, out_r] = pObj.processChunk( pObj.Input{1,1}.Data('new'),...
                    []);
            end
            
            pObj.Output{1}.appendChunk(out_l);
            
            if ~isempty(out_r)
                pObj.Output{2}.appendChunk(out_r);
            end
            
        end
        
        function prepareForProcessing(pObj)
            
            fs = pObj.FsHzIn;
            
            % Filter instantiation (if needed)
            if pObj.bRemoveDC
                pObj.dcFilter_l = bwFilter(fs,4,pObj.cutoffHzDC,[],'high');
                pObj.dcFilter_r = bwFilter(fs,4,pObj.cutoffHzDC,[],'high');
            else
                pObj.dcFilter_l = [];
                pObj.dcFilter_r = [];
            end
            
            if pObj.bPreEmphasis
                pObj.preEmphFilter_l = genericFilter([1 -abs(pObj.coefPreEmphasis)],1,fs);
                pObj.preEmphFilter_r = genericFilter([1 -abs(pObj.coefPreEmphasis)],1,fs);
            else
                pObj.preEmphFilter_l = [];
                pObj.preEmphFilter_r = [];
            end
            
            if pObj.bNormalizeRMS
                a = [1 -exp(-1/(pObj.intTimeSecRMS*fs))];
                b = sum(a);
                pObj.agcFilter_l = genericFilter(b,a,fs);
                pObj.agcFilter_r = genericFilter(b,a,fs);
            else
                pObj.agcFilter_l = [];
                pObj.agcFilter_r = [];
            end
            
            if pObj.bMiddleEarFiltering
                switch pObj.middleEarModel
                    case 'jepsen'
                        model = 'jepsenmiddleear';
                    otherwise
                        model = pObj.middleEarModel;
                end
                a = 1;
                b = middleearfilter(fs, model);
                pObj.midEarFilter_l = genericFilter(b,a,fs);
                pObj.midEarFilter_r = genericFilter(b,a,fs);
            else
                pObj.midEarFilter_l = [];
                pObj.midEarFilter_r = [];
            end
            
        end
        
    end
    
    methods (Static)
        
        function dep = getDependency()
            dep = 'input';
        end
        
        function [names, defaultValues, descriptions] = getParameterInfo()
            %getParameterInfo   Returns the parameter names, default values
            %                   and descriptions for that processor
            %
            %USAGE:
            %  [names, defaultValues, description] =  preProc.getParameterInfo;
            %
            %OUTPUT ARGUMENTS:
            %         names : Parameter names
            % defaultValues : Parameter default values
            %  descriptions : Parameter descriptions
            
            
            names = {'pp_bRemoveDC',...
                    'pp_cutoffHzDC',...
                    'pp_bPreEmphasis',...
                    'pp_coefPreEmphasis',...
                    'pp_bNormalizeRMS',...
                    'pp_bBinauralRMS',...
                    'pp_intTimeSecRMS',...
                    'pp_bLevelScaling',...
                    'pp_refSPLdB',...
                    'pp_bMiddleEarFiltering',...
                    'pp_middleEarModel',...
                    'pp_bUnityComp'};
            
            descriptions = {'Flag to activate DC-removal filter',...
                    'Cutoff frequency (Hz) of DC-removal high-pass filter',...
                    'Flag to activate the pre-emphasis high-pass filter',...
                    'Coefficient for pre-emphasis compensation (usually between 0.9 and 1)',...
                    'Flag for activating automatic gain control',...
                    'Flag indicating the use of unified automatic gain control over left and right channel, for preserving channel relative differences.',...
                    'Time constant (s) for automatic gain control',...
                    'Flag to apply level scaling to the given reference',...
                    'Reference dB SPL value to correspond to input signal RMS value of 1',...
                    'Flag to apply middle ear filtering',...
                    'Middle ear filter model (jepsen or lopezpoveda)',...
                    'Compensation to have maximum of unity gain for middle ear filter (automatically true for Gammatone and false for drnl filterbanks)'};
            
            defaultValues = {0,...
                            20,...
                            0,...
                            0.97,...
                            0,...
                            1,...
                            500E-3,...
                            0,...
                            100,...
                            0,...
                            'jepsen',...
                            1};
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'Pre-processor';
            pInfo.label = 'Pre-processing stage';
            pInfo.requestName = 'time';
            pInfo.requestLabel = 'Time domain signal';
            pInfo.outputType = 'TimeDomainSignal';
            pInfo.isBinaural = 2;   % Indicates that the processor can behave as mono or binaural
            
        end
        
    end
    
    % "Getter" methods
    methods
        function value = get.bRemoveDC(pObj)
            value = pObj.parameters.map('pp_bRemoveDC');
        end
        
        function value = get.cutoffHzDC(pObj)
            value = pObj.parameters.map('pp_cutoffHzDC');
        end
        
        function value = get.bPreEmphasis(pObj)
            value = pObj.parameters.map('pp_bPreEmphasis');
        end
        
        function value = get.coefPreEmphasis(pObj)
            value = pObj.parameters.map('pp_coefPreEmphasis');
        end
        
        function value = get.bNormalizeRMS(pObj)
            value = pObj.parameters.map('pp_bNormalizeRMS');
        end
        
        function value = get.intTimeSecRMS(pObj)
            value = pObj.parameters.map('pp_intTimeSecRMS');
        end
        
        function value = get.bLevelScaling(pObj)
            value = pObj.parameters.map('pp_bLevelScaling');
        end
        
        function value = get.refSPLdB(pObj)
            value = pObj.parameters.map('pp_refSPLdB');
        end
        
        function value = get.bMiddleEarFiltering(pObj)
            value = pObj.parameters.map('pp_bMiddleEarFiltering');
        end
        
        function value = get.middleEarModel(pObj)
            value = pObj.parameters.map('pp_middleEarModel');
        end
        
        function value = get.bBinauralRMS(pObj)
            value = pObj.parameters.map('pp_bBinauralRMS');
        end
        
        function value = get.bUnityComp(pObj)
            value = pObj.parameters.map('pp_bUnityComp');
        end
        
        
    end
    
end