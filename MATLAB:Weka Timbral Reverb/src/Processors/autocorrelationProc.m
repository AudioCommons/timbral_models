classdef autocorrelationProc < Processor
%AUTOCORRELATIONPROC Auto-correlation processor.
%   This processor calculates the auto-correlation in the fast Fourier
%   transform domain for short time frames based on the inner hair-cell
%   representation, towards prediction of perceived pitch [1].
%
%   AUTOCORRELATIONPROC properties:
%        wname       - Window shape descriptor (see window.m)
%        wSizeSec    - Window duration in seconds
%        hSizeSec    - Step size between windows in seconds
%        clipMethod  - Center clipping method ('clc','clp','sgn') [2]
%        alpha       - Threshold coefficient in center clipping
%        K           - Exponent to control the compression [3]
%
%   See also: Processor, ihcProc
%
%   Reference:
%   [1] Meddis, R. and Hewitt, M. J. (1991), "Virtual pitch and phase 
%       sensitivity of a computer model of the auditory periphery. 
%       I: Pitch identification," Journal of the Acoustical Society
%       of America 89(6), pp. 2866?2882.
%   [2] Rabiner, L. R. (1977), "On the use of autocorrelation analysis 
%       for pitch detection," IEEE Transactions on Audio, Speech, and 
%       Language Processing 25(1), pp. 24?33.
%   [3] Tolonen, T. and Karjalainen, M. (2000), "A computationally efficient
%       multipitch analysis model," IEEE Transactions on Audio, Speech, and
%       Language Processing 8(6), pp. 708?716.

    properties (Dependent = true)
        wname       % Window shape descriptor (see window.m)
        wSizeSec    % Window duration in seconds
        hSizeSec    % Step size between windows in seconds
        lags        % Vector of lags
        clipMethod  % Center clipping method ('clc','clp','sgn')
        alpha       % Threshold coefficient in center clipping
        K           % Exponent in auto-correlation
    end
    
    properties (GetAccess = private)
        wSize       % Window duration in samples
        hSize       % Step size between windows in samples
        win         % Window vector
        buffer      % Buffered input signals for framing
        do_mex      % Flag indicating the use of the Tobias' mex code
    end
    
    methods
        function pObj = autocorrelationProc(fs,parObj,do_mex)
        %autocorrelationProc   Construct an autocorrelation processor
        %
        % USAGE:
        %   pObj = autocorrelationProc(fs, parObj)
        %   pObj = autocorrelationProc(fs, parObj, do_mex)
        %
        % INPUT ARGUMENTS:
        %     fs : Input sampling frequency (Hz)
        % parObj : Parameter object instance
        % do_mex : Set to 0 to disable use of pre-compiled mex file in the computation
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
            pObj = pObj@Processor(fs,[],'autocorrelationProc',parObj);
            
            if nargin>0 % Safeguard for Matlab empty calls
                
                pObj.do_mex = do_mex;
                
                % Initialize buffer
                pObj.buffer = [];
            end
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
            %
            %NOTE: This method does not control dimensionality of the
            %provided input. If called outside of a manager instance,
            %validity of the input is the responsibility of the user!
            
           
            % Append the input to existing buffer
            if ~isempty(pObj.buffer)
                in = [pObj.buffer;in];
            end
            
            % Get dimensionality of buffered input
            [nSamples,nChannels] = size(in);
            
            % How many frames are in the buffered input?
            nFrames = floor((nSamples-(pObj.wSize-pObj.hSize))/pObj.hSize);
            
            % Determine maximum lag
            M = pObj.wSize;     % Frame size in samples
            maxLag = M-1;      % Maximum lag in computation

            % Pre-allocate output
            out = zeros(nFrames,nChannels,maxLag);            
            
            if ~pObj.do_mex
                % Loop on the frames
                for ii = 1:nFrames
                    % Get start and end indexes for the current frame
                    n_start = (ii-1)*pObj.hSize+1;
                    n_end = (ii-1)*pObj.hSize+pObj.wSize;

                    % Loop on the channel
                    for jj = 1:nChannels

                        % Extract current frame
                        frame = pObj.win.*in(n_start:n_end,jj);
                        
                        % Perform center clipping
                        frame = applyCenterClipping(frame,pObj.clipMethod,pObj.alpha);
                        
                        % Compute auto-correlation:

                        % Get the frame in the Fourier domain
                        XX = abs(fft(frame,2^nextpow2(2*M-1))).^pObj.K;

                        % Back to time domain
                        x = real(ifft(XX));

                        % Normalize by auto-correlation at lag zero
                        x = x/x(1);

                        % Store results for positive lags only
                        out(ii,jj,:) = x(1:M);

                    end

                end
                
            else
                % Use Tobias previous code
                
                % Loop on the auditory channels
                for jj = 1:nChannels
                    
                    % Framing using mex
                    frames = frameData(in(:,jj),pObj.wSize,pObj.hSize,pObj.win,false);
                    
                    % Perform center clipping
                    frames = applyCenterClipping(frames,pObj.clipMethod,pObj.alpha);
                    
                    % Auto-correlation analysis
                    acf = calcACorr(frames,maxLag,'unbiased',pObj.K);
                    
                    % Normalize by lag zero
                    acf = acf ./ repmat(acf(1,:),[M-1 1]) ;
                    
%                     acf = acorrNorm(frames,maxLag-1,true);
%                     
%                     % ACF of window
%                     acfWin = calcACorr(pObj.win,maxLag,'coeff',pObj.K);
%                     
%                     % Normalize ACF pattern
%                     acf = acf ./ repmat(acfWin + 1E-10,[1 nFrames]);
                                        
                    % Store results for positive lags only
                    out(:,jj,:) = permute(acf,[2 3 1]);
                end
            end

            % Update the buffer: the input that was not extracted as a
            % frame should be stored
            pObj.buffer = in(nFrames*pObj.hSize+1:end,:);
            
        end
        
        function reset(pObj)
            %reset      Resets the auto-correlation processor (cleans the
            %           internal buffer)
            %
            %USAGE
            %   pObj.reset
            %   pObj.reset()
            %
            %INPUT ARGUMENT
            % pObj : Auto-correlation processor instance
            
            % Empty the buffer
            pObj.buffer = [];
        end
        
    end
    
    methods (Hidden = true)
        
        function prepareForProcessing(pObj)
            
            % Compute internal parameters
            pObj.wSize = 2*round(pObj.parameters.map('ac_wSizeSec')*pObj.FsHzIn/2);
            pObj.hSize = round(pObj.parameters.map('ac_hSizeSec')*pObj.FsHzIn);
            pObj.win = window(pObj.parameters.map('ac_wname'),pObj.wSize);
            % Output sampling frequency
            pObj.FsHzOut = 1/(pObj.hSizeSec);
            
        end
        
    end
    
    % "Getter" methods
    methods
        function wSizeSec = get.wSizeSec(pObj)
            wSizeSec = pObj.parameters.map('ac_wSizeSec');
        end
        
        function hSizeSec = get.hSizeSec(pObj)
            hSizeSec = pObj.parameters.map('ac_hSizeSec');
        end
        
        function wname = get.wname(pObj)
            wname = pObj.parameters.map('ac_wname');
        end
        
        function clipMethod = get.clipMethod(pObj)
            clipMethod = pObj.parameters.map('ac_clipMethod');
        end
        
        function alpha = get.alpha(pObj)
            alpha = pObj.parameters.map('ac_clipAlpha');
        end
        
        function K = get.K(pObj)
            K = pObj.parameters.map('ac_K');
        end
        
        function lags = get.lags(pObj)
            lags = ((1:(2 * round(pObj.wSizeSec * pObj.FsHzIn * 0.5)-1))-1)/pObj.FsHzIn;
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
            
            
            names = {'ac_wname',...
                    'ac_wSizeSec',...
                    'ac_hSizeSec',...
                    'ac_clipMethod',...
                    'ac_clipAlpha',...
                    'ac_K'};
            
            descriptions = {'Window name',...
                    'Window duration (s)',...
                    'Window step size (s)',...
                    'Center clipping method (''clc'', ''clp'', or ''sgn'')',...
                    'Threshold in center clipping (between 0 and 1)',...
                    'Exponent in auto-correlation'};
            
            defaultValues = {'hann',...
                            20E-3,...
                            10E-3,...
                            'clp',...
                            0.6,...
                            2};
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'Auto-correlation';
            pInfo.label = 'Auto-correlation';
            pInfo.requestName = 'autocorrelation';
            pInfo.requestLabel = 'Autocorrelation computation';
            pInfo.outputType = 'CorrelationSignal';
            pInfo.isBinaural = false;
            
        end
        
    end
    
end