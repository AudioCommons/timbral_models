classdef gaborProc < Processor
%GABORPROC Spectro-temporal modulation spectrogram.
%   Models the spectro-temporal receptive field of neurons through two-dimensional Gabor
%   functions [1].
%
%   See also: Processor, modulationProc
%
%   Reference:
%   [1] Schädler, M. R., Meyer, B. T., and Kollmeier, B. (2012), "Spectro-temporal
%       modulation subspace-spanning filter bank features for robust automatic speech
%       recognition," Journal of the Acoustical Society of America 131(5), pp. 4134-4151.
    
% TODO: Limit the parameters of the dependent ratemap to the values this processor is
% "callibrated" with. Move the definition of pObj.nFeat to a specific method.

    properties (Dependent = true)
        maxDynamicRangeDB   % Used to limit the dynamic range of input ratemap
    end
    
    properties (SetAccess = private)
        nFeat               % Number of Gabor features
    end
    
    methods
        function pObj = gaborProc(fs,parObj)
        %gaborProc   Construct a Gabor feature extraction processor
        %
        % USAGE:
        %   pObj = gaborProc(fs, parObj)
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
            pObj = pObj@Processor(fs,fs,'gaborProc',parObj);
            
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
            
            if isempty( in )
                out = [];
                return;
            end
            
            % Maximum ratemap power
            max_pow = max(in(:));
            
            % Minimum ratemap floor to limit dynamic range
            min_pow = db2pow(-(pObj.maxDynamicRangeDB + (0 - pow2db(max_pow))));
            
            % Apply static compression
            in = pow2db(in + min_pow);

            % Compute Gabor features
            gb_feat = gbfb(in.');
            
            % Normalize features
            out = normalizeData(gb_feat','meanvar');
            
        end
        
        function reset(~)
            % Nothing to reset for that processor at the moment..
        end
        
        function output = instantiateOutput(pObj,dObj)
            %INSTANTIATEOUTPUT  Instantiate the output signal for this processor
            %
            %NB: This method is overloaded here from the master Processor class, as
            %feature signals need additional input argument to construct
            
            nChanIn = size(pObj.getDependentParameter('fb_cfHz'),2);
            pObj.nFeat = size(gbfb(ones(nChanIn,1)),1);
            
            featureNames = cell(1,pObj.nFeat);
            for jj = 1:pObj.nFeat
                featureNames{jj} = num2str(jj);
            end
            
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
        
        function verifyParameters(~)
           
            % TODO: This implementation is designed to function only on ratemaps computed
            % with a given window size/overlap. 
            
        end
        
    end
    
    % "Getter" methods
    methods
        function maxDynamicRangeDB = get.maxDynamicRangeDB(pObj)
            maxDynamicRangeDB = pObj.parameters.map('gb_maxDynamicRangeDB');
        end
    end
    
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
            %                           gammatoneProc.getParameterInfo;
            %
            %OUTPUT ARGUMENTS:
            %         names : Parameter names
            % defaultValues : Parameter default values
            %  descriptions : Parameter descriptions
            
            
            names = {'gb_maxDynamicRangeDB'};
            
            descriptions = {'Maximum dynamic range (dB) of input ratemap'};
            
            defaultValues = {80};
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'Gabor features';
            pInfo.label = 'Gabor features extractor';
            pInfo.requestName = 'gabor';
            pInfo.requestLabel = 'Gabor features extraction';
            pInfo.outputType = 'FeatureSignal';
            pInfo.isBinaural = false;
            
        end
        
    end
    
end