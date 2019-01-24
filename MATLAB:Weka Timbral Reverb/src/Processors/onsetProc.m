classdef onsetProc < Processor
%ONSETPROC Onset processor.
%   The onset processor detects the signal onsets by measuring the
%   frame-based increase in the energy of the ratemap representation, and
%   computes the onset strength as a function of time frame and frequency
%   channel [1,2].
%
%   ONSETPROC properties:
%       maxOnsetdB      - Upper limit for onset strength in dB
%
%   See also: Processor, ratemapProc, offsetProc
%
%   Reference:
%   [1] Bregman, A. S. (1990), Auditory scene analysis: The perceptual 
%       organization of sound, the MIT Press, Cambridge, MA, USA.
%   [2] Klapuri, A. (1999), "Sound onset detection by applying psychoacoustic 
%       knowledge," in Proceedings of the IEEE International Conference on 
%       Acoustics, Speech and Signal Processing (ICASSP), pp. 3089-3092.   
    
    properties (Dependent = true)
        maxOnsetdB      % Upper limit for onset value
    end
    
    properties (GetAccess = private)
        buffer          % Buffered last frame of the previous chunk
    end
    
    methods
        function pObj = onsetProc(fs,parObj)
        %onsetProc   Construct an onset detection processor
        %
        % USAGE:
        %   pObj = onsetProc(fs, parObj)
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
            pObj = pObj@Processor(fs,fs,'onsetProc',parObj);
            
            if nargin>0
                % Initialize an empty buffer
                pObj.buffer = [];
            end
            
        end
        
        function out = processChunk(pObj,in)
            %processChunk   Requests the processing for a new chunk
            %
            %USAGE:
            %   out = pObj.processChunk(in)
            %
            %INPUT ARGUMENTS:
            %   in : Input signal (ratemap)
            %
            %OUTPUT ARGUMENTS:
            %  out : Output signal
            
            % Initialize a buffer if empty
            if isempty( in )
                out = [];
                return;
            end
            
            if isempty(pObj.buffer)
                pObj.buffer = 10*log10(in(1,:));
            end
            
            % Concatenate the input with the buffer
            bufIn = cat(1,pObj.buffer,10*log10(in));
            
            % Compute onset
            onset = diff(bufIn);
            
            % Discard offsets and limit onset strength
            out = min(max(onset,0),pObj.maxOnsetdB);
           
            % Update the buffer
            pObj.buffer = 10*log10(in(end,:));
            
        end
        
        function reset(pObj)
            %reset      Resets the internal buffers of the processor
            %
            %USAGE:
            %   pObj.reset()
            %
            %INPUT PARAMETERS:
            %   pObj : Processor instance
            
            % Reset the buffer
            pObj.buffer = [];
            
        end
                
    end
    
    % "Getter" methods
    methods
        function maxOnsetdB = get.maxOnsetdB(pObj)
            maxOnsetdB = pObj.parameters.map('ons_maxOnsetdB');
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
            
            
            names = {'ons_maxOnsetdB'};
            
            descriptions = {'Upper limit for onset value (dB)'};
            
            defaultValues = {30};
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'Onset detection';
            pInfo.label = 'Onset detection';
            pInfo.requestName = 'onsetStrength';
            pInfo.requestLabel = 'Onset strength';
            pInfo.outputType = 'TimeFrequencySignal';
            pInfo.isBinaural = false;
            
        end
        
    end
    
end