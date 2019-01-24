classdef offsetProc < Processor
%OFFSETPROC Onset processor.
%   The offset processor detects the signal offsets by measuring the
%   frame-based decrease in the energy of the ratemap representation, and
%   computes the offset strength as a function of time frame and frequency
%   channel [1].
%
%   OFFSETPROC properties:
%       maxOffsetdB      - Upper limit for offset strength in dB
%
%   See also: Processor, ratemapProc, onsetProc
%
%   Reference:
%   [1] Bregman, A. S. (1990), Auditory scene analysis: The perceptual 
%       organization of sound, the MIT Press, Cambridge, MA, USA.

    properties (Dependent = true)
        maxOffsetdB      % Upper limit for onset value
    end
    
    properties (GetAccess = private)
        buffer          % Buffered last frame of the previous chunk
    end
    
    methods
        function pObj = offsetProc(fs,parObj)
        %offsetProc   Construct an offset detection processor
        %
        % USAGE:
        %   pObj = offsetProc(fs, parObj)
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
            pObj = pObj@Processor(fs,fs,'offsetProc',parObj);
            
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
            if isempty(pObj.buffer)
                pObj.buffer = 10*log10(in(1,:));
            end
            
            % Concatenate the input with the buffer
            bufIn = cat(1,pObj.buffer,10*log10(in));
            
            % Compute offset
            offset = diff(bufIn);
            
            % Discard onsets and limit onset strength
            out = min(abs(min(offset,0)),abs(pObj.maxOffsetdB));
            
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
        function maxOffsetdB = get.maxOffsetdB(pObj)
            maxOffsetdB = pObj.parameters.map('ofs_maxOffsetdB');
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
            
            
            names = {'ofs_maxOffsetdB'};
            
            descriptions = {'Upper limit for offset value (dB)'};
            
            defaultValues = {30};
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'Offset detection';
            pInfo.label = 'Offset detection';
            pInfo.requestName = 'offsetStrength';
            pInfo.requestLabel = 'Offset strength';
            pInfo.outputType = 'TimeFrequencySignal';
            pInfo.isBinaural = false;
            
        end
        
    end
    
    
end