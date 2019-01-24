classdef itdProc < Processor
%ITDPROC Interaural Time Difference processor.
%   This processor estimates the time difference between the left and the
%   right ear signals for individual frequency channels and time frames by
%   locating the time lag that corresponds to the most prominent peak in
%   the normalized cross-correlation function. This estimation is further 
%   refined by a parabolic interpolation stage [1].
%
%   ITDPROC properties:
%       frameFsHz       - Sampling frequency of the signal (see below)
%
%   See also: Processor, crosscorrelationProc
%
%   Reference:
%   [1] May, T., van de Par, S., and Kohlrausch, A. (2011), "A probabilistic 
%       model for robust localization based on a binaural auditory front-end,¡± 
%       IEEE Transactions on Audio, Speech, and Language Processing 19(1), 
%       pp. 1?13.

    properties (GetAccess = private)
        frameFsHz   % Sampling frequency of the signal before framing, used 
                    %   for expressing ITDs in seconds
    end
    
    methods
        
        function pObj = itdProc(fs,parObj)
        %itdProc   Construct an inter-aural time difference detection processor
        %
        % USAGE:
        %   pObj = itdProc(fs, parObj)
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
            pObj = pObj@Processor(fs,fs,'itdProc',parObj);
            
        end
        
        function out = processChunk(pObj,in)
            %processChunk   Calls the processing for a new chunk of signal
            %
            %USAGE
            %   out = pObj.processChunk(in)
            %
            %INPUT ARGUMENT
            %  pObj : Processor instance
            %    in : Input (cross-correlation)
            %
            %OUTPUT ARGUMENT
            %   out : Corresponding output
            
            
            % Dimensionality of the input
            [nFrames,nChannels,nLags] = size(in);
            
            % Pre-allocate output
            out = zeros(nFrames,nChannels);
            
            % Create a lag vector
            lags = (0:nLags-1).'-(nLags-1)/2;
            
            % Loop over the time frame
            for ii = 1:nFrames
                
                % Loop over the frequency channel
                for jj = 1:nChannels
                    
                    % Find the peak in the discretized crosscorrelation
                    [c,i] = max(in(ii,jj,:));
                    
                    % Lag of most salient peak
                    lagInt = lags(i);
                    
                    if i>1 && i<nLags
                        % Then interpolate using neighbor points
                        c_l = in(ii,jj,i-1);    % Lower neighbor
                        c_u = in(ii,jj,i+1);    % Upper neighbor
                        
                        % Estimate "true" peak deviation through parabolic
                        % interpolation
                        delta = 0.5*(c_l-c_u)/(c_l-2*c+c_u);
                        
                        % Store estimate
                        out(ii,jj) = (lagInt + delta)/pObj.frameFsHz;
                        
                    else
                        % Do not interpolate if the peak is at a boundary
                        out(ii,jj) = lagInt/pObj.frameFsHz;
                    end
                    
                end
            end
        end
        
        function reset(~)
            % Nothing to reset for this processor
        end
        
    end
    
    methods (Hidden = true)
        
        function prepareForProcessing(pObj)
            
            % Get original signal sampling frequency to express ITDs in ms
            pObj.frameFsHz = pObj.LowerDependencies{1}.FsHzIn;
            
        end
        
    end
    
    methods (Static)
        
        function dep = getDependency()
            dep = 'crosscorrelation';
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
            
            
            names = {};
            
            descriptions = {};
            
            defaultValues = {};
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'ITD Extractor';
            pInfo.label = 'ITD Extractor';
            pInfo.requestName = 'itd';
            pInfo.requestLabel = 'Inter-aural time difference';
            pInfo.outputType = 'TimeFrequencySignal';
            pInfo.isBinaural = true;
            
        end
        
    end
    
end