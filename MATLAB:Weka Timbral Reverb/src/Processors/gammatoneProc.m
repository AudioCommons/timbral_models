classdef gammatoneProc < Processor
%GAMMATONEPROC Gammatone auditory filterbank processor.
%   The Gammatone filterbank models the frequency selectivity of the peripheral auditory
%   system according following [1]. It operates on a time-domain signal and returns a
%   time-frequency representation of the signal. 
%
%   GAMMATONEPROC properties:
%       cfHz       - Channels center frequencies (Hz)
%       nERBs      - Distance between neighboring filters in ERBS (see [1])
%       nGamma     - Gammatone order of the filters (2 or 4)
%       bwERBs     - Bandwidth of the filters in ERBs (see [1])
%       lowFreqHz  - Requested center frequency of lowest channel (Hz)
%       highFreqHz - Requested approximate center frequency of highest channel (Hz)
%       bAlign     - Use phase-aligned filters
%       delaySec   - Time delay in seconds
% 
%   There are three different ways of setting up a vector of channel center frequencies
%   (cfHz) when instantiating this processor:
%       1- By providing the lower and upper center frequencies (lowFreqHz and highFreqHz),
%          and the distance between neighboring filters (nERBs).
%       2- By providing the lower and upper center frequencies (lowFreqHz and highFreqHz),
%          and the number of channels that the representation should have.
%       3- By directly providing a vector of center frequencies (cfHz).
%   In case of conflicting arguments, cfHz is generated from one of the three method above
%   with priority order 3 > 2 > 1.
%
%   See also: Processor, drnlProc
%
%   Reference:
%   [1] Glasberg, B.R. and Moore, B.C.J. (1990), "Derivation of auditory filter shapes
%       from notched-noise data", Hearing Research 47(1-2), pp. 103-138.
    
    properties (Dependent = true)
        cfHz            % Filters center frequencies
        nERBs           % Distance between neighboring filters in ERBs
        nGamma          % Gammatone order of the filters
        bwERBs          % Bandwidth of the filters in ERBs
        lowFreqHz       % Lowest center frequency used at instantiation
        highFreqHz      % Highest center frequency used at instantiation
        bAlign          % Use phase-aligned filters
        delaySec        % Time delay in seconds
    end
    
    properties (GetAccess = private)
        Filters         % Array of filter objects
    end
        
    methods
        function pObj = gammatoneProc(fs,parObj)
        %gammatoneProc   Construct a Gammatone filterbank processor
        %
        % USAGE:
        %   pObj = gammatoneProc(fs, parObj)
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
            
            % TODO: 
            %  - Implement solution to allow for different impulse response
            %    durations for different filters (if necessary)
            %  - Implement bAlign option
            
            if nargin<2||isempty(parObj); parObj = Parameters; end
            if nargin<1; fs = []; end
            
            % Call super-constructor
            pObj = pObj@Processor(fs,fs,'gammatoneProc',parObj);
            
            if nargin>0 && ~isempty(fs)
                
            end
        end
        
        function out = processChunk(pObj,in)
            %processChunk       Passes an input signal through the
            %                   Gammatone filterbank
            %
            %USAGE
            %       out = processChunk(pObj,in)
            %       out = pObj.processChunk(in)
            %
            %INPUT ARGUMENTS
            %      pObj : Gammatone filterbank object
            %        in : One-dimensional array containing the input signal
            %
            %OUTPUT ARGUMENTS
            %       out : Multi-dimensional array containing the filterbank
            %             outputs
            %
            %SEE ALSO:
            %       gammatoneProc.m
            
            % TO DO: Indicate that this function is not buit to deal with
            % multiple channels. Multiple channels should be treated with
            % multiple instances of the filterbank.
            
            % Check inputs
            if min(size(in))>1
                error('The input should be a one-dimensional array')
            end
            
            % Turn input into column vector
            in = in(:);
            
            % Get number of channels
            nFilter = size(pObj.Filters,2);
            
            % Pre-allocate memory
            out = zeros(length(in),nFilter);
            
            % Loop on the filters
            for ii = 1:nFilter
                out(:,ii) = pObj.Filters(ii).filter(in);
            end
            
            % TO DO : IMPLEMENT ALIGNMENT CORRECTION
            
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
            
            nFilter = size(pObj.Filters,2);
            
            % Resetting the internal states of the filters
            for ii = 1:nFilter
                pObj.Filters(ii).reset();
            end
            
        end
        
        function bInBranch = isSuitableForRequest(pObj)
            if strcmp(pObj.parameters.map('fb_type'),'gammatone')
                bInBranch = true;
            else
                bInBranch = false;
            end
        end
        
    end
    
    methods (Access=protected)
        
        function verifyParameters(pObj)
                        
            % Solve the conflicts between center frequencies, number of channels, and
            % distance between channels
            if ~isempty(pObj.parameters.map('fb_cfHz'))
                % Highest priority case: a vector of channels center 
                %   frequencies is provided
                centerFreq = pObj.parameters.map('fb_cfHz');
                
                pObj.parameters.map('fb_lowFreqHz') = centerFreq(1);
                pObj.parameters.map('fb_highFreqHz') = centerFreq(end);
                pObj.parameters.map('fb_nChannels') = numel(centerFreq);
                pObj.parameters.map('fb_nERBs') = 'n/a';
                
                
            elseif ~isempty(pObj.parameters.map('fb_nChannels'))
                % Medium priority: frequency range and number of channels
                %   are provided
               
                % Build a vector of center ERB frequencies
                ERBS = linspace( freq2erb(pObj.parameters.map('fb_lowFreqHz')), ...
                                freq2erb(pObj.parameters.map('fb_highFreqHz')), ...
                                pObj.parameters.map('fb_nChannels') );  
                centerFreq = erb2freq(ERBS);    % Convert to Hz
                
                pObj.parameters.map('fb_nERBs') = (ERBS(end)-ERBS(1)) ...
                                                / pObj.parameters.map('fb_nChannels');
                pObj.parameters.map('fb_cfHz') = centerFreq;
                
                
            else
                % Lowest (default) priority: frequency range and distance 
                %   between channels is provided (or taken by default)
                
                % Build vector of center ERB frequencies
                ERBS = freq2erb(pObj.parameters.map('fb_lowFreqHz')): ...
                                double(pObj.parameters.map('fb_nERBs')): ...
                                            freq2erb(pObj.parameters.map('fb_highFreqHz'));
                centerFreq = erb2freq(ERBS);    % Convert to Hz
                
                pObj.parameters.map('fb_nChannels') = numel(centerFreq);
                pObj.parameters.map('fb_cfHz') = centerFreq;
                
            end
            
            % Calculate filter bandwidth in Hertz
            bwHz = pObj.parameters.fb_bwERBs * (24.7 + 0.108 * centerFreq);
            
            % Map bandwidth to time delay in seconds
            pObj.parameters.map('fb_delaySec') = 3./(2*pi*bwHz);
        end
    end
    
    methods (Hidden = true)
        
        function prepareForProcessing(pObj)
            
            % Instantiate filters
            pObj.Filters = pObj.populateFilters;
            
        end
        
    end
    
    % "Getter" methods
    methods
        function cfHz = get.cfHz(pObj)
            cfHz = pObj.parameters.map('fb_cfHz');
        end
        
        function nERBs = get.nERBs(pObj)
            nERBs = pObj.parameters.map('fb_nERBs');
        end
        
        function nGamma = get.nGamma(pObj)
            nGamma = pObj.parameters.map('fb_nGamma');
        end
        
        function bwERBs = get.bwERBs(pObj)
            bwERBs = pObj.parameters.map('fb_bwERBs');
        end
        
        function lowFreqHz = get.lowFreqHz(pObj)
            lowFreqHz = pObj.parameters.map('fb_lowFreqHz');
        end
        
        function highFreqHz = get.highFreqHz(pObj)
            highFreqHz = pObj.parameters.map('fb_highFreqHz');
        end
        
        function bAlign = get.bAlign(pObj)
            bAlign = pObj.parameters.map('fb_bAlign');
        end
        
        function delaySec = get.delaySec(pObj)
            delaySec = pObj.parameters.map('fb_delaySec');
        end
        
    end
    
    
    methods (Access = private)
        function obj = populateFilters(pObj)
            % This function is a workaround to assign an array of objects
            % as one of the processor's property, should remain private

            fs = pObj.FsHzIn;
            cfHz = pObj.parameters.map('fb_cfHz');
            n = pObj.parameters.map('fb_nGamma');
            bw = pObj.parameters.map('fb_bwERBs');
            bAlign = pObj.parameters.map('fb_bAlign');
            nFilter = numel(cfHz);
            
            % Preallocate memory by instantiating last filter
            obj(1,nFilter) = gammatoneFilter(cfHz(nFilter),fs,n,...
                                        bw,bAlign);
            % Instantiating remaining filters
            for ii = 1:nFilter-1
                obj(1,ii) = gammatoneFilter(cfHz(ii),fs,n,...
                                        bw,bAlign);
            end                        
            
        end
    end

    methods (Static)
        
        function dep = getDependency()
            dep = 'time';
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
            
            names = {'fb_type',...
                    'fb_lowFreqHz',...
                    'fb_highFreqHz',...
                    'fb_nERBs',...
                    'fb_nChannels',...
                    'fb_cfHz',...
                    'fb_nGamma',...
                    'fb_bwERBs',...
                    'fb_bAlign',...
                    'fb_delaySec'};
            
            descriptions = {'Filterbank type (''gammatone'' or ''drnl'')',...
                    'Lowest center frequency (Hz)',...
                    'Highest center frequency (Hz)',...
                    'Distance between neighbor filters in ERBs',...
                    'Number of channels',...
                    'Channels center frequencies (Hz)',...
                    'Gammatone rising slope order',...
                    'Bandwidth of the filters (ERBs)',...
                    'Create phase-aligned filters',...
                    'Time delay in seconds'};
            
            defaultValues = {'gammatone',...
                            80,...
                            8000,...
                            1,...
                            [],...
                            [],...
                            4,...
                            1.018,...
                            false,...
                            []};
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'Gammatone filterbank';
            pInfo.label = 'Gammatone filterbank';
            pInfo.requestName = 'filterbank';
            pInfo.requestLabel = 'Gammatone filterbank output';
            pInfo.outputType = 'TimeFrequencySignal';
            pInfo.isBinaural = false;
            
        end
        
    end
        
end