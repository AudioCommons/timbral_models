classdef CorrelationSignal < Signal
%CORRELATIONSIGNAL Signal class for three-dimensional correlation signals.
%   This class collects all signals resulting from a correlation computation on a
%   time-frequency representation in short time windows (e.g., auto-correlation, 
%   cross-correlation). Its data is therefore three dimensional, with first to third
%   dimension respectively related to time, frequency, and lag.
%
%   CORRELATIONSIGNAL properties:
%       cfHz - Center frequencies of the frequency channels (Hz)
%       lags - Lag values (in seconds)
%
% See also Signal, crosscorrelationProc, autocorrelationProc
    
    properties (SetAccess=protected)
        cfHz    % Center frequencies of the frequency channels (Hz)
        lags    % Lag values 
    end
    
    methods
        function sObj = CorrelationSignal(procHandle,bufferSize,channel,data)
            %CorrelationSignal  Class constructor
            %
            %USAGE
            %     sObj = CorrelationSignal(procHandle)
            %     sObj = CorrelationSignal(procHandle,bufferSize,channel,data)
            %
            %INPUT ARGUMENTS
            % procHandle : Handle to the processor generating this signal as output
            % bufferSize : Size of the ring buffer in s (default: bufferSize = 10)
            %    channel : Flag indicating 'left', 'right', or 'mono' (default: 
            %              channel = 'mono')
            %       data : Array of amplitudes to construct an object from existing data
            %
            %OUTPUT ARGUMENT
            %  sObj : Correlation signal object inheriting the signal class
            
            if nargin<4; data = []; end
            if nargin<3||isempty(channel); channel = 'mono'; end
            if nargin<2||isempty(bufferSize); bufferSize = 10; end
            if nargin<1||isempty(procHandle); procHandle = emptyProc; end
            
            % TODO: Might have to change how lags are accessed to prevent error when
            % instantiating empty processor
            numChannel = max(length(procHandle.getDependentParameter('fb_cfHz')),1);
            sObj = sObj@Signal( procHandle, bufferSize, ...
                                            [numChannel length(procHandle.lags)]);
            
            if nargin>0     % Safeguard for Matlab empty calls
                
            % Populate object properties
            sObj.Dimensions = 'nSample x nFilters x nLags';
            sObj.cfHz = procHandle.getDependentParameter('fb_cfHz');
            sObj.lags = procHandle.lags;
            sObj.setData( data );
            sObj.Channel = channel;
                
            end
        end
        
        function h = plot(sObj,h0,p,frameNb,varargin)
            %plot       This method plots the data from a correlation signal object
            %
            %USAGE
            %       sObj.plot
            %       sObj.plot(h_prev,p,frameNb)
            %       h = sObj.plot(...)
            %
            %INPUT ARGUMENT
            %  h_prev : Handle to an already existing figure or subplot
            %           where the new plot should be placed
            %       p : Structure of non-default plot parameters (generated
            %           from genParStruct.m)
            % frameNb : Specify a given time frame to plot. If none specified, the summary
            %           correlation will be plotted as a 2D image.
            %
            %OUTPUT ARGUMENT
            %       h : Handle to the newly created figure
            %
            %OPTIONAL ARGUMENTS
            % 'noTitle' - Flag to avoid displaying a title
            
            % TODO: Add an option to plot the correlation in a given frame instead of
            % the summary
            
            % Manage plotting parameters
            if nargin < 3 || isempty(p) 
                % Get default plotting parameters (same as time-frequency signals)
                    p = Parameters.getPlottingParameters('TimeFrequencySignal');
                else
                    defaultPar = Parameters.getPlottingParameters('TimeFrequencySignal');
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
            
            % Compute the summary correlation
            scorr = squeeze(mean(sObj.Data(:),2));
            
            % Time axis
            t = 0:1/sObj.FsHz:(size(sObj.Data(:),1)-1)/sObj.FsHz;
            
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
            if nargin<4 || isempty(frameNb)
                
                % Set the colormap
                try
                    colormap(p.map('colormap'))
                catch
                    warning('No colormap %s is available, using ''jet''.',p.map('colormap'))
                    colormap('jet')
                end

                % Plot
                imagesc(t,sObj.lags,scorr.');
                axis xy

                if p.map('bColorbar')
                    colorbar
                end

                xlabel('Time (s)','fontsize',p.map('fsize_label'),...
                                  'fontname',p.map('ftype'))
                ylabel('Lag period (s)','fontsize',p.map('fsize_label'),...
                                        'fontname',p.map('ftype'))

                % Set up a title
                if ~strcmp(sObj.Channel,'mono')
                    pTitle = [sObj.Label ' summary -' sObj.Channel];
                else
                    pTitle = [sObj.Label ' summary'];
                end
                
                if ~isempty(opt) && isfield(opt,'noTitle')
                    if opt.noTitle ~= 1
                        title(pTitle,'fontsize',p.map('fsize_title'),...
                                     'fontname',p.map('ftype'))
                    end
                end
            
            else
                
                ax(1) = subplot(4,1,[1:3]);
                waveplot(permute(sObj.Data(frameNb,:,:),[3 1 2]),sObj.lags*1E3,...
                    sObj.cfHz,p.map('wavPlotZoom'),1)
                xlabel('')
                hy = ylabel('Center frequency (Hz)');
%                 hypos = get(hy,'position');
%                 hypos(1) = -2.15;
%                 set(hy,'position',hypos);
                if ~isempty(opt) && isfield(opt,'noTitle')
                    if opt.noTitle ~= 1
                        title(sObj.Label,'fontsize',p.map('fsize_title'),...
                                         'fontname',p.map('ftype'))
                    end
                else
                    title(sObj.Label,'fontsize',p.map('fsize_title'),...
                                     'fontname',p.map('ftype'))
                end
                
                ax(2) = subplot(4,1,4);
                plot(sObj.lags*1E3,mean(permute(sObj.Data(frameNb,:,:),[3 1 2]),3),...
                                        'k','linewidth',1.25)
                grid on
                xlim([sObj.lags(1)*1E3 sObj.lags(end)*1E3])
                ylim([0 1])
                xlabel('Lag period (ms)')
                ylabel('Summary')
                
            end
            
        end
        
    end
    
end