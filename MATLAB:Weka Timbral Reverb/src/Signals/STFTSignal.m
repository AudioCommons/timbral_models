classdef STFTSignal < Signal
%STFTSignalSTFT Signal class for complex two-dimensional, time-frequency representations.
%   see circVBufC and TimeFrequenySignal
    
    properties (SetAccess=protected)
        cfHz        % Center frequencies of the frequency channels
    end
       
    properties (GetAccess = protected)
        scaling
    end
    
    methods
        
        function sObj = STFTSignal(procHandle,bufferSize,channel,data,cfHz)
            %STFTSignal   Class constructor
            %
            %USAGE
            %     sObj = STFTSignal(procHandle,bufferSize,channel,data,freqs)
            %
            %INPUT ARGUMENTS
            % procHandle : Handle to the processor generating this signal as output
            % bufferSize : Size of the ring buffer in s (default: bufferSize = 10)
            %    channel : Flag indicating 'left', 'right', or 'mono' (default:
            %              channel = 'mono')
            %       data : Array of amplitudes to construct an object from existing data
            %
            %OUTPUT ARGUMENT
            %  sObj : Time-frequency signal object inheriting the signal class
            
            if nargin<5; cfHz = procHandle.cfHz; end
            if nargin<4; data = []; end
            if nargin<3||isempty(channel); channel = 'mono'; end
            if nargin<2||isempty(bufferSize); bufferSize = 10; end
            if nargin<1||isempty(procHandle); procHandle = emptyProc; end
            
            sObj = sObj@Signal( procHandle, bufferSize, length(cfHz));
            % buffer will store double reals at this point but have correct dimensions.
            % we will call the overloaded setBufferSize to get a complex buffer
            sObj.setBufferSize(bufferSize);
            
            if nargin>0     % Safeguard for Matlab empty calls
                
                sObj.Dimensions = 'nSamples x nFilters';
                sObj.cfHz = cfHz;
                sObj.setData( data );
                sObj.Channel = channel;
                sObj.scaling = 'magnitude';
                
            end
            
        end
        
        function setBufferSize( sObj, newBufferSize_s )
            %setBufferSize  This method sets the buffer to a new size,
            %               erasing all data previously stored.
            %USAGE:
            %   sObj.setBufferSize( newBufferSize_s )
            %
            %INPUT ARGUMENTS:
            %   newBufferSize_s : new size of the buffer in seconds. The
            %                     dimensionality of the individual elements
            %                     remains the same
            
            bufferSizeSamples = ceil( newBufferSize_s * sObj.FsHz );
            sObj.Buf = circVBufC( bufferSizeSamples, sObj.Buf.matSz );
            sObj.Data = circVBufArrayInterface( sObj.Buf );
        end
        
        function h = plot(sObj,h0,p,varargin)
            %plot       This method plots the magnitude of the data from 
            %           a time-frequency domain signal object
            %
            %USAGE
            %       sObj.plot
            %       sObj.plot(h_prev,p,...)
            %       h = sObj.plot(...)
            %
            %INPUT ARGUMENT
            %  h_prev : Handle to an already existing figure or subplot
            %           where the new plot should be placed
            %       p : Structure of non-default plot parameters (generated
            %           from genParStruct.m)
            %
            %OUTPUT ARGUMENT
            %       h : Handle to the newly created figure
            %
            %OPTIONAL ARGUMENTS
            % 'rangeSec' - Vector of time limits for the plot
            
            if ~isempty(sObj.Data)
            
                % Decide if the plot should be on a linear or dB scale
                switch sObj.Name
                    case {'filterbank','ild','ic','itd','onsetStrength',...
                            'offsetStrength','innerhaircell','adaptation','onsetMap',...
                            'offsetMap','drnl'}
                        do_dB = 0;
                    case {'ratemap', 'stft'}
                        do_dB = 1;
                    otherwise 
                        error('Cannot plot this object')
                end
            
                % Manage plotting parameters
                if nargin < 3 || isempty(p) 
                    % Get default plotting parameters
                    p = Parameters.getPlottingParameters('STFTSignal');
                else
                    defaultPar = Parameters.getPlottingParameters('STFTSignal');
                    defaultPar.replaceParameters(p);
                    p = defaultPar;
                end
                
                % Manage optional arguments
                if nargin>3 && ~isempty(varargin)
                    opt = struct;
                    for ii = 1:2:size(varargin,2)
                        opt.(varargin{ii}) = varargin{ii+1};
                    end
                else
                    opt = [];
                end
                
                if do_dB
                    if strcmp(sObj.scaling,'power')
                        data = 10*log10(abs(sObj.Data(:).'));
                    else
                        % Get the data in dB
                        data = 20*log10(abs(sObj.Data(:).'));
                    end
                else
                    data = sObj.Data(:).';
                end

                % Limit to a certain time period
                if ~isempty(opt) && isfield(opt,'rangeSec')
                    data = data(:,floor(opt.rangeSec(1)*sObj.FsHz):floor(opt.rangeSec(end)*sObj.FsHz));
                    t = opt.rangeSec(1):1/sObj.FsHz:opt.rangeSec(1)+(size(data,2)-1)/sObj.FsHz;
                else
                    t = 0:1/sObj.FsHz:(size(data,2)-1)/sObj.FsHz;
                end

                % Manage handles
                if nargin < 2 || isempty(h0)
                        h = figure;             % Generate a new figure
                    elseif get(h0,'parent')~=0
                        % Then it's a subplot
                        figure(get(h0,'parent')),subplot(h0)
                        h = h0;
                    else
                        figure(h0)
                        h = h0;
                end

                % Managing frequency axis ticks for auditory filterbank
                %
                % Find position of y-axis ticks
                M = size(sObj.cfHz,2);  % Number of channels
                n_points = 500;         % Number of points in the interpolation
                interpolate_ticks = spline(1:M,sObj.cfHz,...
                    linspace(0.5,M+0.5,n_points));
                %
                % Restrain ticks to signal range (+/- a half channel)
                aud_ticks = p.map('aud_ticks');
                aud_ticks=aud_ticks(aud_ticks<=interpolate_ticks(end));
                aud_ticks=aud_ticks(aud_ticks>=interpolate_ticks(1));
                n_ticks = size(aud_ticks,2);        % Number of ticks
                ticks_pos = zeros(size(aud_ticks)); % Tick position
                %
                % Find index for each tick
                for ii = 1:n_ticks
                    jj = find(interpolate_ticks>=aud_ticks(ii),1);
                    ticks_pos(ii) = jj*M/n_points;
                end

                
                % Set the color map
                try
                    colormap(p.map('colormap'))
                catch
                    warning('No colormap %s is available, using ''jet''.',p.map('colormap'))
                    colormap('jet')
                end
                
                % Plot the figure
                switch sObj.Name
                    
                    case {'filterbank','innerhaircell','drnl','adaptation'}
                        waveplot(data(:,1:p.map('wavPlotDS'):end).',t(1:p.map('wavPlotDS'):end),sObj.cfHz,p.map('wavPlotZoom'),1);
                    
                    otherwise
                        imagesc(t,1:M,data)  % Plot the data
                        axis xy              % Use Cartesian coordinates
                        
                        % Set up y-axis
                        set(gca,'YTick',ticks_pos,...
                            'YTickLabel',aud_ticks,'fontsize',p.map('fsize_axes'),...
                            'fontname',p.map('ftype'))
                
                        if p.map('bColorbar')
                            colorbar             % Display a colorbar
                        end
                end

                % Set up a title
                if strcmp(sObj.Name,'ratemap')
                    if strcmp(sObj.scaling,'power')
                        % label = [sObj.Label ' (power)'];
                        label = sObj.Label;
                    else
                        % label = [sObj.Label ' (magnitude)'];
                        label = sObj.Label;
                    end
                elseif strcmp(sObj.Name,'stft')
                    label = [sObj.Label, ' magnitude'];
                else
                    label = sObj.Label;
                end
                if ~strcmp(sObj.Channel,'mono')
                    pTitle = [label ' - ' sObj.Channel];
                else
                    pTitle = label;
                end
                
                % Set up axes labels
                xlabel('Time (s)','fontsize',p.map('fsize_label'),'fontname',p.map('ftype'))
                ylabel('Frequency (Hz)','fontsize',p.map('fsize_label'),'fontname',p.map('ftype'))
                title(pTitle,'fontsize',p.map('fsize_title'),'fontname',p.map('ftype'))

                % Set up plot properties
                

                % Scaling the plot
                switch sObj.Name
                    case {'innerhaircell','ratemap', 'stft'}
                        m = max(data(:));    % Get maximum value for scaling
                        set(gca,'CLim',[m-p.map('dynrange') m])

                    case {'ild','itd'}
                        m = max(abs(data(:)))+eps;
                        set(gca,'CLim',[-m m])

                    case 'ic'
                        set(gca,'CLim',[0 1])

                end
            else
                error('This is an empty signal, cannot be plotted')
            end
                
            
        end
    end
    
    methods (Static)
%        
%         function sObj = construct(fs,bufferSize,name,label,cfHz,channel,data)
%             %construct
%             %
%             %
%             
%             if nargin<7; data = []; end
%             if nargin<6; channel = []; end
%             if nargin<5; cfHz = []; end
%             if nargin<4; label = []; end
%             if nargin<3; name = []; end
%             if nargin<2||isempty(bufferSize); bufferSize = 10; end
%             if nargin<1; fs = []; end
%             
%             % Create a dummy structure with that information to emulate a processor and
%             % correctly call the class constructor
%             dummyStruct = struct;
%             dummyStruct.FsHzOut = fs;
%             dummyStruct.getProcessorInfo.requestName = name;
%             dummyStruct.getProcessorInfo.requestLabel = label;
%             dummyStruct.getDependentParameter = containers.Map;
%             dummyStruct.getDependentParameter('fb_cfHz') = cfHz;
%             
%             % Instantiate the signal
%             sObj = STFTSignal(dummyStruct,bufferSize,channel,data);
%             
%         end
%         
        function [names, defaultValues, descriptions] = getPlottingParameterInfo()
            %GETPLOTTINGPARAMETERINFO   Stores plot parameters that are common to all
            %signals.
            
            
            names = {'colormap',...
                    'bColorbar',...
                    'dynrange',...
                    'aud_ticks',...
                    'wavPlotDS',...
                    'wavPlotZoom',...
                    'binaryMaskColor'};
                 
            descriptions = {'Colormap for time-frequency plots',...
                    'Boolean for displaying colorbar in time-frequency plots',...
                    'Dynamic range for time-frequency plots (dB)',...
                    'Auditory ticks for ERB-based representations',...
                    'Decimation ratio for plotting undecimated wave plot representations',...
                    'Zoom factor in wave plot representations',...
                    'Color for binary mask (in RGB value)'};
                
            defaultValues = {'jet',...
                    1,...
                    80,...
                    [100 250 500 1000 2000 4000 8000 16000 32000],...
                    3,...
                    5,...
                    [0 0 0]};
            
        end
%         
    end

end
