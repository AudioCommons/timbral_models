classdef filterObj < handle
    % Filter object class inherited from handle master class to allow for
    % calling by reference properties, hence simplifying syntax
    
    properties (GetAccess=public)       % Public properties
        Type                % Descriptor for type of filter (string)
        RealTF = true;      % True if the transfer function is real-valued
        CascadeOrder =1     % Number of times the filter has to be cascaded
    end
   
    properties %(GetAccess=protected)    % Protected properties
        FsHz        % Sampling frequency in Hertz
        b           % Numerator coefficients of transfer function
        a           % Denominator coefficients of transfer function
        States      % Current states of filter
        nChan       % Number of channels
    end
    
    properties (Dependent)              % Dependent properties
        Order       % Filter order
    end
    
    methods
        function [hLin,f] = frequencyResponse(fObj,nfft)
            %calcFilterResponse   Compute frequency response of filter objects.
            %
            %USAGE 
            %   [hlin,f] = calcFilterResponse(fObj,nfft)
            %
            %INPUT ARGUMENTS
            %   fObj : filter object(s) (see genFilterObj)
            %   nfft : resolution of frequency response
            %
            %OUTPUT ARGUMENTS
            %   hlin : Complex frequency response [nfft x 1]
            %      f : Frequency vector           [nfft x 1] 
            
            % CHECK INPUT ARGUMENTS 
            % 
            % 
            % Check for proper input arguments
            if nargin < 1 || nargin > 2
                help(mfilename);
                error('Wrong number of input arguments!')
            end

            % Set default frequency resolution
            if nargin < 2 || isempty(nfft); 
                if ~isempty(fObj.FsHz); 
                    % Sampling frequency is known
                    nfft = 2^nextpow2(fObj.FsHz * 50e-3); 
                else
                    nfft = 512; 
                end
            end


            % *********************  COMPUTE FREQUENCY RESPONSE  *********************

            % Number of filter
%             nFilter = length(fObj);

            % Allocate memory
            hLin = zeros(nfft,1);

            
            if fObj.FsHz
                if fObj.RealTF
                    % Calculate frequency response of ii-th filter
                    [hLin(:),f] = freqz(fObj.b,fObj.a,nfft,fObj.FsHz);
                else
                    % TO DO: Is there a more elegant way to do things here?
                    ir = impz(fObj.b,fObj.a,nfft,fObj.FsHz);
                    [hLin(:),f] = freqz(2*real(ir),1,nfft,fObj.FsHz);
                end
            else
                % Calculate frequency response of ii-th filter
                [hLin(:),f] = freqz(fObj.b,fObj.a,nfft);
            end
        end
        
        function h = plot(fObj,nfft,mode,handle)
            %plot   Plot frequency/phase response of a filter object.
            %
            %USAGE 
            %      h = fObj.plot()
            %      h = fObj.plot(nfft,mode)
            %      h = plot(fObj,nfft,mode)
            %
            %INPUT ARGUMENTS
            %   fObj : filter object 
            %   nfft : resolution of frequency response (default, 
            %          nfft = 512)
            %   mode : string specifying what to plot (default, mode = 'both')
            %          'magnitude' - plot magnitude response
            %              'phase' - plot phase response 
            %               'both' - plot both magnitude and phase response
            %            'impulse' - plot filter's impulse response
            % handle : handle to a previous filter object plot to visualize
            %          alongside other filters
            %
            %OUTPUT ARGUMENTS
            %   h : figure handle
            % 
            % TO DO: Finish use of handle for plotting on existing figure
            
            % CHECK INPUT ARGUMENTS 
            % 
            % 
            % Check for proper input arguments
            if nargin < 1 || nargin > 4
                help(mfilename);
                error('Wrong number of input arguments!')
            end

            % Set default values
            if nargin < 2 || isempty(nfft); nfft = [];     end
            if nargin < 3 || isempty(mode); mode = 'both'; end
            if nargin < 4 || isempty(handle); handle = []; end


            % COMPUTE FREQUENCY RESPONSE  
            % 
            % 
            % Compute filter response
            [hLin,f] = frequencyResponse(fObj,nfft);

            % Magnitude response in dB
            hdB = 20 * log10(abs(hLin));

            % Compute 10th and 90th percentile
            pct    = prctile(hdB,[5 95]);
            yRange = [-50 20];

            % Check if majority of data is within predefined range
            bSetY = pct(1) > yRange(1) && pct(2) < yRange(2);


            % PLOT FREQUENCY RESPONSE
            % 
            % 
            % Select mode
            switch lower(mode)
                case 'magnitude'
                    h = figure;
                    if fObj.FsHz
                        semilogx(f,hdB)
                        xlabel('Frequency (Hz)');
                        xlim([10 f(end)])
                    else
                        plot(f,hdB)
                        xlabel('Normalized frequency (x \pi rad/sample)');
                        xlim([f(1) f(end)])            
                    end
                    ylabel('Magnitude (dB)');
                    grid on;

                    % Show filter label
                    % TO DO: Investigate Type vs. Label
                    if ~isempty(fObj.Type); title(fObj.Type); end

                    if bSetY; ylim(yRange); end

                case 'phase'
                    phi = unwrap(angle(hLin)) * 180/pi;

                    h = figure;

                    if fObj.FsHz
                        semilogx(f,phi)
                        xlabel('Frequency (Hz)');
                        xlim([10 f(end)])
                    else
                        plot(f/pi,phi)
                        xlabel('Normalized frequency (x \pi rad/sample)');
                        xlim([0 0.5])            
                    end
                    ylabel('Phase (degree)');
                    grid on;

                    % Show filter label  % TO DO: Type vs. label
                    if ~isempty(fObj.Type); title(fObj.Type); end


                case 'both'
                    phi = unwrap(angle(hLin)) * 180/pi;

                    if isempty(handle)
                        h = figure;
                    else
                        h = figure(handle);
                        hold on
                    end
                    ax(1) = subplot(211);
                    if fObj.FsHz
                        semilogx(f,hdB)
                        xlabel('Frequency (Hz)');
                        xlim([10 f(end)])
                    else
                        plot(f/pi,hdB)
                        xlabel('Normalized frequency (x \pi rad/sample)');
                        xlim([0 0.5])            
                    end
                    ylabel('Magnitude (dB)');
                    grid on;

                    % Show filter label above the fist subplot % TO DO:
                    % Type vs. label
                    if ~isempty(fObj.Type); title(fObj.Type); end

                    if bSetY; ylim(yRange); end

                    ax(2) = subplot(212);

                    if fObj.FsHz
                        semilogx(f,phi)
                        xlabel('Frequency (Hz)');
                        xlim([10 f(end)])
                    else
                        plot(f/pi,phi)
                        xlabel('Normalized frequency (x \pi rad/sample)');
                        xlim([0 0.5])            
                    end
                    ylabel('Phase (degree)');
                    grid on;

                    linkaxes(ax,'x')
                otherwise
                    error(['Plotting mode "',lower(mode),'" is not supported.'])
            end
        end
        
        function out = filter(fObj,data)
            %filterAudio   Perform digital filtering of data/audio objects.
            %
            %USAGE 
            %   out = fObj.filterAudio(in)
            %
            %INPUT ARGUMENTS
            %     in : single channel audio 
            %   fObj : filter object
            %
            %OUTPUT ARGUMENTS
            %    out : filtered audio object or data matrix
            
            % N.B: As fObj is inherited from the handle master class, the
            % filter properties (e.g., its states), will be updated
            
            % Check for proper input arguments
            if nargin ~= 2
                help(mfilename);
                error('Wrong number of input arguments!')
            end

            % Update size of data
            if isempty(fObj.nChan)
                fObj.nChan = size(data,2);
%                 fObj.reset
            end

            % Check if filter states are initialized
            if isempty(fObj.States) 
                % Initialize filter states
                fObj.reset;
            elseif fObj.Order ~= size(fObj.States,1)
                error(['Dimension mismatch between the filter ',...
                       'order and the filter states.']);
            end
            
            % Processing
            if fObj.CascadeOrder == 1
                % Then apply the filter only once
                [out,fObj.States] = filter(fObj.b,fObj.a,data,fObj.States,1);
                
            else
                % Then cascade the filter
                out = data;
                for ii = 1:fObj.CascadeOrder
                    [out,fObj.States(:,:,ii)] = filter(fObj.b,fObj.a,out,fObj.States(:,:,ii),1);
                end
            end
            
            % Correction for complex-valued transfer function filters
            if ~(fObj.RealTF)
                out = 2*real(out);
            end
            
        end

        function order = get.Order(fObj)
            if isempty(fObj.a)||isempty(fObj.b)
                order = [];
            else
                order = max(length(fObj.b),length(fObj.a))-1;
            end
        end
        
        function reset(fObj,states)
            % This method resets the filter's states to zero or to states if specified
            if nargin<2 || isempty(states)
                states = [];
            end
            if isempty(fObj.Order)
                error('The filter transfer function must have been specified before initializing its states')
            else
                if isempty(states)
                    % Create filter states
                    fObj.States = zeros(fObj.Order,fObj.nChan,fObj.CascadeOrder);
                else
                    fObj.States = states;
                end
            end
        end 
        
        function bInit = isInitialized(fObj)
            %isInitialized      Returns true when filter states are non-empty
            
            if isempty(fObj.States)
                bInit = 0;
            else
                bInit = 1;
            end
        end
        
    end
    
    methods (Access=protected)
        function fObj = populateProperties(fObj,varargin)
            
            % First check on input
            if mod(size(varargin,2),2)||isempty(varargin)
                error('Additional input arguments have to come in pairs of ...,''property name'',value,...')
            end
            
            % List of valid properties % TO DO: should this be hardcoded
            % here?
            validProp = {'Type',...
                         'RealTF',...
                         'FsHz',...
                         'b',...
                         'a',...
                         'States',...
                         'Order'};
                     
            % Loop on the additional arguments
            for ii = 1:2:size(varargin,2)-1
                % Check that provided property name is a string
                if ~ischar(varargin{ii})
                    error('Property names should be given as strings, %s isn''t one!',num2str(varargin{ii}))
                end
                % Check that provided property name is valid
                if ~ismember(varargin{ii},validProp)
                    error('Property name ''%s'' is invalid',varargin{ii})
                end
                % Then add the property value
                fObj.(varargin{ii})=varargin{ii+1};
            end
            
            
        end 
    end
    
    
end