classdef crosscorrelationProc < Processor
%CROSSCORRELATIONPROC Cross-correlation processor.
%   The cross-correlation between the right and left ears is computed from
%   their respective Inner Hair-Cell representations, in the Fourier 
%   domain for given time frames. This is normalized by the
%   auto-correlation sequence at lag zero, and then evaluated for time lags
%   in a given range, resulting in a three-dimensional function of time
%   frame, frequency channel and lag time.
%
%   CROSSCORRELATIONPROC properties:
%        wname       - Window type 
%        wSizeSec    - Window duration in seconds
%        hSizeSec    - Step size between windows in seconds
%        maxDelaySec - Maximum delay in cross-correlation computation (s)
%
%   See also: Processor, ihcProc
%
    
    properties (Dependent = true)
        wname       % Window shape descriptor (see window.m)
        wSizeSec    % Window duration in seconds
        hSizeSec    % Step size between windows in seconds
        lags        % Vector of lags at which cross-correlation is computed
        maxDelaySec % Maximum delay in cross-correlation computation (s)
    end
    
    properties (GetAccess = private)
        wSize       % Window duration in samples
        hSize       % Step size between windows in samples
        win         % Window vector
        buffer_l    % Buffered input signals (left ear)
        buffer_r    % Buffered input signals (right ear)
        do_mex      % TEMP flag indicating the use of the Tobias' mex code (1)
    end
        
    methods
        
        function pObj = crosscorrelationProc(fs,parObj,do_mex)
		%crosscorrelationProc   Construct a cross-correlation processor
        %
        % USAGE:
        %   pObj = crosscorrelationProc(fs, parObj)
        %   pObj = crosscorrelationProc(fs, parObj, do_mex)
        %
        % INPUT ARGUMENTS:
        %     fs : Input sampling frequency (Hz)
        % parObj : Parameter object instance
        % do_mex : Set to 0 to disable use of pre-compiled mex file in computation
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
            if nargin<3||isempty(do_mex);do_mex = 1; end
            if nargin<2||isempty(parObj); parObj = Parameters; end
            if nargin<1; fs = []; end
            
            % Call superconstructor
            pObj = pObj@Processor(fs,[],'crosscorrelationProc',parObj);
            
            if nargin>0     % Safeguard for Matlab empty calls
            
                pObj.do_mex = do_mex;
                
                % Initialize buffers
                pObj.buffer_l = [];
                pObj.buffer_r = [];
            end 
        end
        
        function out = processChunk(pObj,in_l,in_r)
            %processChunk   Calls the processing for a new chunk of signal
            %
            %USAGE
            %   out = pObj.processChunk(in_l,in_r)
            %
            %INPUT ARGUMENTS
            %  pObj : Processor instance
            %  in_l : Left-ear input (inner hair-cell envelope)
            %  in_r : Right-ear input (inner hair-cell envelope)
            %
            %OUTPUT ARGUMENT
            %   out : Resulting output
            
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
            
            % Determine maximum lag in samples
            maxLag = ceil(pObj.maxDelaySec*pObj.FsHzIn);
            
            % Pre-allocate output
            if ~pObj.do_mex

            out = zeros(nFrames,nChannels,maxLag*2+1);
            % Loop on the time frame
            for ii = 1:nFrames
                % Get start and end indexes for the current frame
                n_start = (ii-1)*pObj.hSize+1;
                n_end = (ii-1)*pObj.hSize+pObj.wSize;
                
                % Loop on the channel
                for jj = 1:nChannels
                    
                    % Extract frame for left and right input
                    frame_l = pObj.win.*in_l(n_start:n_end,jj);
                    frame_r = pObj.win.*in_r(n_start:n_end,jj);
                    
                    % Compute the frames in the Fourier domain
                    X = fft(frame_l,2^nextpow2(2*pObj.wSize-1));
                    Y = fft(frame_r,2^nextpow2(2*pObj.wSize-1));
                    
                    % Compute cross-power spectrum
                    XY = X.*conj(Y);
                    
                    % Back to time domain
                    c = real(ifft(XY));
                    
                    % Adjust to requested maximum lag and move negative
                    % lags upfront
                    if maxLag >= pObj.wSize
                        % Then pad with zeros
                        pad = zeros(maxLag-pObj.wSize+1,1);
                        c = [pad;c(end-pObj.wSize+2:end);c(1:pObj.wSize);pad];
                    else
                        % Else keep lags lower than requested max
                        c = [c(end-maxLag+1:end);c(1:maxLag+1)];
                    end
                    
                    % Normalize with autocorrelation at lag zero and store
                    % output
                    powL = sum(frame_l.^2);
                    powR = sum(frame_r.^2);
                    out(ii,jj,:) = c/sqrt(powL*powR+eps);
                    
                end
                
            end

            else
            out = zeros(max(0,nFrames),nChannels,maxLag*2+1);
                % Use Tobias mex code for framing
                for jj = 1:nChannels
                    % Framing
                    frames_L = frameData(in_l(:,jj),pObj.wSize,pObj.hSize,pObj.win,false);
                    frames_R = frameData(in_r(:,jj),pObj.wSize,pObj.hSize,pObj.win,false);

                    % Cross-correlation analysis
                    output = calcXCorr(frames_L,frames_R,maxLag,'coeff');
                    
                    % Store output
                    out(:,jj,:) = permute(output,[2 3 1]);
                    
                end
            end
            
            % Update the buffer: the input that was not extracted as a
            % frame should be stored
            pObj.buffer_l = in_l(nFrames*pObj.hSize+1:end,:);
            pObj.buffer_r = in_r(nFrames*pObj.hSize+1:end,:);
            
        end
        
        function reset(pObj)
             %reset     Resets the internal states of the processor
             %
             %USAGE
             %      pObj.reset
             %
             %INPUT ARGUMENTS
             %  pObj : Cross-correlation processor instance
             
             % Only thing needed to reset is to empty the buffers
             pObj.buffer_l = [];
             pObj.buffer_r = [];
             
        end
        
    end
    
    methods (Hidden = true)
        
        function prepareForProcessing(pObj)
            
            % Compute internal parameters
            pObj.wSize = 2*round(pObj.parameters.map('cc_wSizeSec')*pObj.FsHzIn/2);
            pObj.hSize = round(pObj.parameters.map('cc_hSizeSec')*pObj.FsHzIn);
            pObj.win = window(pObj.parameters.map('cc_wname'),pObj.wSize);
            % Output sampling frequency
            pObj.FsHzOut = 1/(pObj.hSizeSec);
        
        end
            
    end
    
    % "Getter" methods
    methods
        function wSizeSec = get.wSizeSec(pObj)
            wSizeSec = pObj.parameters.map('cc_wSizeSec');
        end
        
        function hSizeSec = get.hSizeSec(pObj)
            hSizeSec = pObj.parameters.map('cc_hSizeSec');
        end
        
        function wname = get.wname(pObj)
            wname = pObj.parameters.map('cc_wname');
        end
        
        function lags = get.lags(pObj)
            maxLag = ceil(pObj.maxDelaySec*pObj.FsHzIn);
            lags = (-maxLag:maxLag).'./pObj.FsHzIn;
        end
        
        function maxDelaySec = get.maxDelaySec(pObj)
            maxDelaySec = pObj.parameters.map('cc_maxDelaySec');
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
            
            
            names = {'cc_wname',...
                    'cc_wSizeSec',...
                    'cc_hSizeSec',...
                    'cc_maxDelaySec'};
            
            descriptions = {'Window name',...
                    'Window duration (s)',...
                    'Window step size (s)',...
                    'Maximum delay in cross-correlation computation (s)'};
            
            defaultValues = {'hann',...
                            20E-3,...
                            10E-3,...
                            1.1E-3};
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'Cross-correlation';
            pInfo.label = 'Cross-correlation';
            pInfo.requestName = 'crosscorrelation';
            pInfo.requestLabel = 'Crosscorrelation computation';
            pInfo.outputType = 'CorrelationSignal';
            pInfo.isBinaural = true;
            
        end
        
    end
    
    
end