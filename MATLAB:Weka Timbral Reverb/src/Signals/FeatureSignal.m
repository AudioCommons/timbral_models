classdef FeatureSignal < Signal
%FEATURESIGNAL Signal class for groups of related time-domain signals ('features')
%   This class is used for grouping related individual time-domain features into a unique
%   signal (e.g., spectral features). Data is stored in a matrix where first dimension
%   corresponds to time and the second dimension stores different features. The features
%   are labeled by the fList property of the signal.
%
%   FEATURESIGNAL properties:
%       fList - Cell array of features names, ordered as features are arranged in the data
%
% See also Signal, spectralFeaturesProc, pitchProc, gaborProc
    
    properties (SetAccess=protected)
        fList       % Ordered list of the features (cell array of strings)
    end
    
    methods
        
        function sObj = FeatureSignal(procHandle,bufferSize,channel,data,fList)
            %FeatureSignal  Class constructor
            %
            %USAGE
            %     sObj = FeatureSignal(procHandle)
            %     sObj = FeatureSignal(procHandle,bufferSize,channel,data,fList)
            %
            %INPUT ARGUMENTS
            % procHandle : Handle to the processor generating this signal as output
            % bufferSize : Size of the ring buffer in s (default: bufferSize = 10)
            %    channel : Flag indicating 'left', 'right', or 'mono' (default: 
            %              channel = 'mono')
            %       data : Array of amplitudes to construct an object from existing data
            %      fList : Cell array of feature names
            %
            %OUTPUT ARGUMENT
            %  sObj : Feature signal object inheriting the signal class
            
            if nargin<5; fList = []; end
            if nargin<4; data = []; end
            if nargin<3||isempty(channel); channel = 'mono'; end
            if nargin<2||isempty(bufferSize); bufferSize = 10; end
            if nargin<1||isempty(procHandle); procHandle = emptyProc; end
            
            sObj = sObj@Signal( procHandle, bufferSize, size(fList,2));
            
            if nargin>0     % Failsafe for Matlab empty calls
                
                sObj.Dimensions = ['nSamples x ' num2str(size(fList,2)) 'features'];
                sObj.setData( data );
                sObj.Channel = channel;
                sObj.fList = fList;
                
            end
            
        end
        
        function h = plot(sObj,h0,feature,varargin)
            %plot   Plots the requested spectral features
            %
            %USAGE:
            %     sObj.plot
            % h = sObj.plot(mObj,h0,feature)
            %
            %INPUT ARGUMENTS:
            %    sObj : Spectral features signal instance
            %      h0 : Handle to already existing figure or subplot
            % feature : Name of a specific feature to plot
            %
            %OUTPUT ARGUMENTS:
            %       h : Handle to the figure
            %
            %OPTIONAL ARGUMENTS:
            % Keyvalues:
            % 'overlay'  - Handle to a signal object to plot together with the feature
            % 'pitchRange' - Vector of valid pitch range
            % 'confThresh' - Confidence threshold in percent of the maximum
            % 'lagDomain'  - True for plotting 1/pitch (i.e., in the lag domain)
            % 'noSubplots' - True to obtain each feature in its own figure 
            
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
            
            % Manage parameters
            % TODO: Do we want to use common plot parameters (e.g., fontsize)
            
            % Manage optional arguments
            if nargin>3 && ~isempty(varargin)
                opt = struct;
                for ii = 1:2:size(varargin,2)
                    opt.(varargin{ii}) = varargin{ii+1};
                end
            else
                opt = [];
            end
            
            if nargin<3||isempty(feature)
                feature = sObj.fList;
            end
            
            if ~iscell(feature)
                feature = {feature};
            end
            
            if ~strcmp(sObj.Name,'gabor')
            
                % Number of subplots
                nFeatures = size(feature,2);

                % Time axis
                tSec = 0:1/sObj.FsHz:(size(sObj.Data(:,:),1)-1)/sObj.FsHz;

                % Plots
                for ii = 1 : nFeatures

                    % Find the feature
                    jj = find(ismember(sObj.fList,feature{ii}),1);

                    if ~isempty(jj)

                        % Create a subplot if more than one representation is needed
                        if nFeatures > 1
                            if ~isempty(opt) && isfield(opt,'noSubPlots')
                                if opt.noSubPlots == 1
                                    if ii>1
                                        h(ii) = figure;
                                    end
                                    ax(ii) = subplot(1,1,1);
                                else
                                    nSubplots = ceil(sqrt(nFeatures));
                                    ax(ii) = subplot(nSubplots,nSubplots,ii);
                                end
                            else
                                nSubplots = ceil(sqrt(nFeatures));
                                ax(ii) = subplot(nSubplots,nSubplots,ii);
                            end
                        end

                        % Raw plot
                        hp = plot(tSec,sObj.Data(:,jj));

                        % Some feature dependent styling below
                        switch feature{ii}
                            % Spectral features...
                            case {'variation' 'brightness' 'flatness' 'entropy' 'rolloff' 'spread' 'centroid'}

                                if ~isempty(opt) && isfield(opt,'overlay')
                                    hold on;
                                    imagesc(tSec,(0.5:size(opt.overlay.cfHz,2)-0.5)/size(opt.overlay.cfHz,2),10*log10(opt.overlay.Data(:)'));axis xy;
                                end

                                %Repeat plot to be on top
                                hp = plot(tSec,sObj.Data(:,jj));
                                box on

                                % Linestyle
                                set(hp,'LineStyle','-','LineWidth',1,'Color','k')

                                xlim([tSec(1) tSec(end)])
                                xlabel('Time (s)')
                                ylabel('Normalized frequency')
                                ylim([0 1])
                                title(['Spectral ',sObj.fList{ii}])

                            case {'irregularity' 'hfc' 'skewness' 'kurtosis' 'flux' 'decrease' 'crest'}

                                % Linestyle
                                set(hp,'LineStyle','-','LineWidth',1,'Color','k')

                                xlim([tSec(1) tSec(end)])

                                xlabel('Time (s)')
                                ylabel('Feature magnitude')
                                title(['Spectral ',sObj.fList{ii}])


                            % Pitch features
                            case 'pitch'

                                if ~isempty(opt) && isfield(opt,'lagDomain')
                                    if opt.lagDomain
                                        % Plot in terms of lag period
                                        set(hp,'YData',1./get(hp,'YData'));
                                    end
                                end

                                % Linestyle
                                set(hp,'marker','o','markerfacecolor','k','color','k','linestyle','none')

                                xlabel('Time (s)')
                                ylabel('Frequency (Hz)')
                                title('Estimated pitch contour')

                                if ~isempty(opt) && isfield(opt,'pitchRange')
                                    ylim(opt.pitchRange)
                                end

                            case 'rawPitch'

                                if ~isempty(opt) && isfield(opt,'lagDomain')
                                    if opt.lagDomain
                                        % Plot in terms of lag period
                                        set(hp,'YData',1./get(hp,'YData'));
                                    end
                                end

                                % Linestyle
                                set(hp,'marker','x','markerfacecolor','k','color','k',...
                                    'linestyle','none','markersize',8,'linewidth',2)

                                % Valid pitch indication
                                if ~isempty(opt) && isfield(opt,'pitchRange')
                                    rangeLags = 1./opt.pitchRange;
                                    plot([tSec(1) tSec(end)],[rangeLags(1) rangeLags(1)],'w--','linewidth',2)
                                    plot([tSec(1) tSec(end)],[rangeLags(2) rangeLags(2)],'w--','linewidth',2)
                                end

                            case 'confidence'

                                set(hp,'LineStyle','-','color','k','linewidth',1.25)

                                % Plot the maximum
                                [maxVal,maxIdx] = max(sObj.Data(:,jj));
                                hold on
                                plot(tSec(maxIdx),maxVal,'rx','linewidth',2,'markersize',12);

                                % And the threshold if available
                                if ~isempty(opt) && isfield(opt,'confThres')
                                    plot([tSec(1) tSec(end)],[opt.confThres opt.confThres],'--k','linewidth',2);
                                    legend({'SACF magnitude' 'global maximum' 'confidence threshold'},'location','southeast');
                                else
                                    legend({'SACF magnitude' 'global maximum'},'location','southeast');
                                end

                                xlabel('Time (s)')
                                ylabel('Magnitude')
                                ylim([0 1.1*maxVal])
                                title('Confidence measure')
                        end


                    else
                        warning('There is no feature names %s in the signal',feature{ii})
                    end




                end
                if nFeatures > 1
                    linkaxes(ax,'x');
                end
                set(gca,'xLim',[0 tSec(end)])
    %             set(h,'units','normalized','outerposition',[0 0 1 1])
            else
                
                % Special plot for Gabor features
                
                % Time axis
                tSec = 0:1/sObj.FsHz:(size(sObj.Data(:,:),1)-1)/sObj.FsHz;
                
                imagesc(tSec,1:size(sObj.Data(:),2),sObj.Data(:).');axis xy;
                xlim([tSec(1) tSec(end)])
                colorbar;
                xlabel('Time (s)')
                ylabel('# feature dimensions')
                title('Gabor features')
                
                
            end
        end
        
    end
    
end