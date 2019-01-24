classdef ildProc < Processor
%ILDPROC Interaural Level Difference processor.
%   This processor estimates the level difference between the left and the
%   right ear signals for individual frequency channels and time frames by
%   comparing the frame-based energy of the left and the right-ear inner 
%   hair cell representations. The value is expressed in dB, and negative
%   values indicate sources on the left-hand side.
%
%   ILDPROC properties:
%       wname       - Window shape descriptor
%       wSizeSec    - Window duration in seconds
%       hSizeSec    - Step size between windows in seconds
%
%   See also: Processor, ihcProc
%
    
    properties (Dependent = true)
        wname       % Window shape descriptor (see window.m)
        wSizeSec    % Window duration in seconds
        hSizeSec    % Step size between windows in seconds
    end
    
    properties (GetAccess = private)
        wSize       % Window duration in samples
        hSize       % Step size between windows in samples
        win         % Window vector
        buffer_l    % Buffered input signals (left ear)
        buffer_r    % Buffered input signals (right ear)
    end
    
    methods
        function pObj = ildProc(fs,parObj)
		%ildProc   Construct an inter-aural level differences extractor processor
        %
        % USAGE:
        %   pObj = ildProc(fs, parObj)
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
            pObj = pObj@Processor(fs,[],'ildProc',parObj);
            
            if nargin>0
                % Initializa the buffers
                pObj.buffer_l = [];
                pObj.buffer_r = [];
            end
            
        end
        
        function out = processChunk(pObj,in_l,in_r)
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
            
            % Compute ILDs:
            
            % Pre-allocate output
            out = zeros(nFrames,nChannels);
            
            % Loop on the time frame
            for ii = 1:nFrames
                % Get start and end indexes for the current frame
                n_start = (ii-1)*pObj.hSize+1;
                n_end = (ii-1)*pObj.hSize+pObj.wSize;
                
                % Loop on the channel
                for jj = 1:nChannels
                    
                    % Energy in the windowed frame for left and right input
                    frame_l = mean(power(pObj.win.*in_l(n_start:n_end,jj),2));
                    frame_r = mean(power(pObj.win.*in_r(n_start:n_end,jj),2));
                    
                    % Compute the ild for that frame
                    out(ii,jj) = 10*log10((frame_r+eps)/(frame_l+eps));
                    
                end
                
            end
            
            % Update the buffer: the input that was not extracted as a
            % frame should be stored
            pObj.buffer_l = in_l(nFrames*pObj.hSize+1:end,:);
            pObj.buffer_r = in_r(nFrames*pObj.hSize+1:end,:);
            
            
        end
        
        function reset(pObj)
             %reset     Resets the internal states of the ILD extractor
             %
             %USAGE
             %      pObj.reset
             %
             %INPUT ARGUMENTS
             %  pObj : ILD extractor processor instance
             
             % Only thing needed to reset is to empty the buffer
             pObj.buffer_l = [];
             pObj.buffer_r = [];
             
        end
        
    end
    
    methods (Hidden = true)
        
        function prepareForProcessing(pObj)
            
            % Compute internal parameters
            pObj.wSize = 2*round(pObj.parameters.map('ild_wSizeSec')*pObj.FsHzIn/2);
            pObj.hSize = round(pObj.parameters.map('ild_hSizeSec')*pObj.FsHzIn);
            pObj.win = window(pObj.parameters.map('ild_wname'),pObj.wSize);
            % Output sampling frequency
            pObj.FsHzOut = 1/(pObj.hSizeSec);
                
        end
        
    end
    
    % "Getter" methods
    methods
        function wSizeSec = get.wSizeSec(pObj)
            wSizeSec = pObj.parameters.map('ild_wSizeSec');
        end
        
        function hSizeSec = get.hSizeSec(pObj)
            hSizeSec = pObj.parameters.map('ild_hSizeSec');
        end
        
        function wname = get.wname(pObj)
            wname = pObj.parameters.map('ild_wname');
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
            
            
            names = {'ild_wname',...
                    'ild_wSizeSec',...
                    'ild_hSizeSec'};
            
            descriptions = {'Window name',...
                    'Window duration (s)',...
                    'Window step size (s)'};
            
            defaultValues = {'hann',...
                            20E-3,...
                            10E-3};
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'ILD Extractor';
            pInfo.label = 'Inter-aural Level Difference Extractor';
            pInfo.requestName = 'ild';
            pInfo.requestLabel = 'Inter-aural level difference';
            pInfo.outputType = 'TimeFrequencySignal';
            pInfo.isBinaural = true;
            
        end
        
    end
        
    
end