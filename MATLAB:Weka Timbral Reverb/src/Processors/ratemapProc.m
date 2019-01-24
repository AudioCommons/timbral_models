classdef ratemapProc < Processor
%RATEMAPPROC Ratemap processor.
%   The ratemap represents a map of auditory nerve firing rates [1], computed
%   from the inner hair-cell signal representation for individual frequency 
%   channels. 
%
%   RATEMAPPROC properties:
%        wname       - Window type
%        wSizeSec    - Window duration
%        hSizeSec    - Window step size
%        scaling     - Flag specifying ratemap scaling
%        decaySec    - Signal-smoothing leaky integrator time constant        
%
%   See also: Processor, ihcProc
%
%   Reference:
%   [1] Brown, G. J. and Cooke, M. P. (1994), "Computational auditory scene
%       analysis," ComputerSpeech and Language 8(4), pp. 297?336.

    properties (Dependent = true)
        wname       % Window shape descriptor (see window.m)
        wSizeSec    % Window duration in seconds
        hSizeSec    % Step size between windows in seconds
        scaling     % Flag specifying 'magnitude' or 'power'
        decaySec    % Integration time constant (seconds)
    end
    
    properties (GetAccess = private)
        wSize       % Window duration in samples
        hSize       % Step size between windows in samples
        win         % Window vector
        buffer      % Buffered input signals
        rmFilter    % Leaky integrator filter
        do_mex      % Flag indicating use of mex code for framing
    end
        
    
    methods
        function pObj = ratemapProc(fs,parObj,do_mex)
		%ratemapProc   Construct a rate-map extraction processor
        %
        % USAGE:
        %   pObj = ratemapProc(fs, parObj)
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
            if nargin<3||isempty(do_mex);do_mex = 1; end
            if nargin<2||isempty(parObj); parObj = Parameters; end
            if nargin<1; fs = []; end
            
            % Call superconstructor
            pObj = pObj@Processor(fs,[],'ratemapProc',parObj);
            
            if nargin>0 % Safeguard for Matlab empty calls
                
                pObj.do_mex = do_mex;
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
            
            
            % Filter the input
            in = pObj.rmFilter.filter(in);
            
            % Append filtered input to the buffer
            if ~isempty(pObj.buffer)
                in = [pObj.buffer;in];
            end

            [nSamples,nChannels] = size(in);
            
            % How many frames are in the buffered input?
            nFrames = floor((nSamples-(pObj.wSize-pObj.hSize))/pObj.hSize);

            % Pre-allocate output
            out = zeros(nFrames,nChannels);
            
            if ~pObj.do_mex
                
                % Loop on the time frame
                for ii = 1:nFrames
                    % Get start and end indexes for the current frame
                    n_start = (ii-1)*pObj.hSize+1;
                    n_end = (ii-1)*pObj.hSize+pObj.wSize;

                    % Loop on the channel
                    for jj = 1:nChannels

                        switch pObj.scaling
                            case 'magnitude'
                                % Averaged magnitude in the windowed frame 
                                out(ii,jj) = mean(pObj.win.*in(n_start:n_end,jj));
                            case 'power'
                                % Averaged energy in the windowed frame for left 
                                out(ii,jj) = mean(power(pObj.win.*in(n_start:n_end,jj),2));
                            otherwise
                                error('Incorrect scaling method for ratemap')
                        end
                    end


                end

            else
                
                % Loop over the auditory channels
                for jj = 1:nChannels
                    
                    % Frame the data in that channel
                    frames = frameData(in(:,jj),pObj.wSize,pObj.hSize,pObj.win,false);

                    % modified by I.T.: There has to be a handling of the case of an
                    % empty frames-matrix.
                    if ~isempty( frames )
                        % Average the samples in the frame
                        switch pObj.scaling
                            
                            case 'magnitude'
                                % Average magnitude in the frame
                                out(:,jj) = mean(frames,1);
                                
                            case 'power'
                                % Average energy in the frame
                                out(:,jj) = mean(frames.^2,1);
                                
                        end
                    end
                    
                end
                
            end
            
            % Update the buffer: the input that was not extracted as a
            % frame should be stored
            pObj.buffer = in(nFrames*pObj.hSize+1:end,:);
            
            
        end
            
        function reset(pObj)
            %reset     Resets the internal states of the ratemap extractor
            %
            %USAGE
            %      pObj.reset
            %
            %INPUT ARGUMENTS
            %  pObj : Ratemap processor instance
            
            % Reset the leaky integrators
            if ~isempty(pObj.rmFilter)
                pObj.rmFilter.reset;
            end
            
            % Empty the buffer
            pObj.buffer = [];
            
        end
        
    end
    
    methods (Hidden = true)
        
        function prepareForProcessing(pObj)
            
            % Compute internal parameters
            pObj.wSize = 2*round(pObj.parameters.map('rm_wSizeSec')*pObj.FsHzIn/2);
            pObj.hSize = round(pObj.parameters.map('rm_hSizeSec')*pObj.FsHzIn);
            pObj.win = window(pObj.parameters.map('rm_wname'),pObj.wSize);
            
            % Output sampling frequency
            pObj.FsHzOut = 1/(pObj.hSizeSec);

            % Initialize filter
            pObj.rmFilter = leakyIntegratorFilter(pObj.FsHzIn,pObj.decaySec);
            
        end
        
    end
    
    % "Getter" methods
    methods
        function wSizeSec = get.wSizeSec(pObj)
            wSizeSec = pObj.parameters.map('rm_wSizeSec');
        end
        
        function hSizeSec = get.hSizeSec(pObj)
            hSizeSec = pObj.parameters.map('rm_hSizeSec');
        end
        
        function wname = get.wname(pObj)
            wname = pObj.parameters.map('rm_wname');
        end
        
        function scaling = get.scaling(pObj)
            scaling = pObj.parameters.map('rm_scaling');
        end
        
        function decaySec = get.decaySec(pObj)
            decaySec = pObj.parameters.map('rm_decaySec');
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
            
            
            names = {'rm_wname',...
                    'rm_wSizeSec',...
                    'rm_hSizeSec',...
                    'rm_scaling',...
                    'rm_decaySec'};
            
            descriptions = {'Window name',...
                    'Window duration (s)',...
                    'Window step size (s)',...
                    'Ratemap scaling (''power'' or ''magnitude'')',...
                    'Leaky integrator time constant (s)'};
            
            defaultValues = {'hann',...
                            20E-3,...
                            10E-3,...
                            'power',...
                            8E-3};
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'Ratemap';
            pInfo.label = 'Ratemap';
            pInfo.requestName = 'ratemap';
            pInfo.requestLabel = 'Ratemap extraction';
            pInfo.outputType = 'TimeFrequencySignal';
            pInfo.isBinaural = false;
            
        end
        
    end
    
    methods (Access = private)
        function obj = populateFilters(pObj,nChannels,fs)
            % This function creates an array of filter objects. It returns
            % the array instead of directly setting up the property in pObj
            % as a workaround to a presumable bug
            
            % Preallocate memory
            obj(1,nChannels) = leakyIntegratorFilter(fs,pObj.decaySec);
            
            % Instantiate one filter per channel
            for ii = 1:nChannels-1
                obj(1,ii) = leakyIntegratorFilter(fs,pObj.decaySec);
            end
            
        end
    end
    
end