classdef DuetHistogramSignal < Signal
%DUETHISTOGRAMSIGNAL 2d histogram signal
%
%   DuetHistogramSignal properties:
%       nAlpha - number of bins for alpha
%       nDelta - number of bisn for delta
%
% See also Signal, duetProc
    
    properties (SetAccess=protected)
        binsAlpha;
        binsDelta;
        maxAlpha;
        maxDelta;
    end
    
    methods
        function sObj = DuetHistogramSignal(procHandle,bufferSize,channel,data)
            %DuetHistogramSignal  Class constructor
            %
            %USAGE
            %     sObj = DuetHistogramSignal(procHandle)
            %     sObj = DuetHistogramSignal(procHandle,bufferSize,channel,data)
            %
            %INPUT ARGUMENTS
            % procHandle : Handle to the processor generating this signal as output
            % bufferSize : Size of the ring buffer in s (default: bufferSize = 10)
            %    channel : Flag indicating 'left', 'right', or 'mono' (default: 
            %              channel = 'mono')
            %       data : Array of amplitudes to construct an object from existing data
            %
            %OUTPUT ARGUMENT
            %  sObj : Duet histigram signal object inheriting the signal class
            
            if nargin<4; data = []; end
            if nargin<3||isempty(channel); channel = 'mono'; end
            if nargin<2||isempty(bufferSize); bufferSize = 10; end
            if nargin<1||isempty(procHandle); procHandle = emptyProc; end
            
            % super constructor
            sObj = sObj@Signal( procHandle, bufferSize,...
                [procHandle.getDependentParameter('duet_binsAlpha') ...
                 procHandle.getDependentParameter('duet_binsDelta')]);
                        
            if nargin>0
                sObj.Dimensions = 'nSamples x binsAlpha x binsDelta';
                sObj.binsAlpha = procHandle.getDependentParameter('duet_binsAlpha');
                sObj.binsDelta = procHandle.getDependentParameter('duet_binsDelta');
                sObj.maxAlpha = procHandle.getDependentParameter('duet_maxAlpha');
                sObj.maxDelta = procHandle.getDependentParameter('duet_maxDelta');
                sObj.setData( data );
                sObj.Channel = channel;                
            end
        end
        
        function h = plot(sObj,h0,p,frameNb,varargin)
            %plot       This method plots the data from a duet histogram signal object
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
            % frameNb : Specify a given time frame to plot. If none specified, the first
            %           one will be plotted.
            %
            %OUTPUT ARGUMENT
            %       h : Handle to the newly created figure
            %
            %OPTIONAL ARGUMENTS
            % 'noTitle' - Flag to avoid displaying a title
            
            % TODO: 
            
            % Manage plotting parameters
            if nargin < 3 || isempty(p)
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
            
            % axis labeling
            alphaAxisVal = linspace(-sObj.maxAlpha, sObj.maxAlpha, sObj.binsAlpha);
            deltaAxisVal = linspace(-sObj.maxDelta, sObj.maxDelta, sObj.binsDelta);
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
            
            
            % Set the colormap
            try
                colormap(p.map('colormap'))
            catch
                warning('No colormap %s is available, using ''jet''.',p.map('colormap'))
                colormap('jet')
            end
            
            % Plot
            if nargin<4 || isempty(frameNb)
                frameNb = 1;
            end
            plot_data = squeeze(sObj.Data(frameNb,:,:));
            mesh(deltaAxisVal, alphaAxisVal, plot_data);
            
            if p.map('bColorbar')
                colorbar
            end
            
            xlabel('delta (relative phase shift)', 'fontsize', p.map('fsize_label'),...
                'fontname',p.map('ftype'));
            ylabel('alpha (symetric attenuation)', 'fontsize', p.map('fsize_label'),...
                'fontname',p.map('ftype'));
            
            % Set up a title
            if ~strcmp(sObj.Channel,'mono')
                pTitle = [sObj.Label '-' sObj.Channel];
            else
                pTitle = sObj.Label;
            end
            pTitle = [pTitle ' - t=' frameNb];
            
            if ~isempty(opt) && isfield(opt,'noTitle')
                if opt.noTitle ~= 1
                    title(pTitle,'fontsize',p.map('fsize_title'),...
                        'fontname',p.map('ftype'))
                end
            end
            
        end
        
    end
    
end
