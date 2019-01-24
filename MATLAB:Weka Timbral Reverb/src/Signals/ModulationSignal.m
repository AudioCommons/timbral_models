classdef ModulationSignal < Signal
%MODULATIONSIGNAL Signal class for three-dimensional amplitude modulation signals
%   This class collects three-dimensional outputs of a modulation filterbank, respectively
%   time, audio frequency, and modulation frequency.
%
%   MODULATIONSIGNAL properties:
%       cfHz    - Center frequencies of audio frequency channels (Hz)
%       modCfHz - Center frequencies of modulation frequency channels (Hz)
%
% See also Signal, modulationProc
    
    properties (SetAccess=protected)
        cfHz    % Audio center frequencies of the channels (Hz)
        modCfHz % Modulation center frequencies (Hz)
    end
    
    methods
        
        function sObj = ModulationSignal(procHandle,bufferSize,channel,data)
            %ModulationSignal  Class constructor
            %
            %USAGE
            %     sObj = ModulationSignal(procHandle)
            %     sObj = ModulationSignal(procHandle,bufferSize,channel,data)
            %
            %INPUT ARGUMENTS
            % procHandle : Handle to the processor generating this signal as output
            % bufferSize : Size of the ring buffer in s (default: bufferSize = 10)
            %    channel : Flag indicating 'left', 'right', or 'mono' (default: 
            %              channel = 'mono')
            %       data : Array of amplitudes to construct an object from existing data
            %
            %OUTPUT ARGUMENT
            %  sObj : Modulation signal object inheriting the signal class
            
            if nargin<4; data = []; end
            if nargin<3||isempty(channel); channel = 'mono'; end
            if nargin<2||isempty(bufferSize); bufferSize = 10; end
            if nargin<1||isempty(procHandle); procHandle = emptyProc; end
            
            cfHz = procHandle.getDependentParameter('fb_cfHz'); %#ok<PROP>
            if isprop(procHandle,'modCfHz')
                modCfHz = procHandle.modCfHz; %#ok<PROP>
            else
                modCfHz = []; %#ok<PROP>
            end
            
            sObj = sObj@Signal( procHandle, bufferSize, ...
                                [length(cfHz) length(modCfHz)]); %#ok<PROP>
            
            if nargin>0     % Safeguard for Matlab empty calls
                
            % Data dimensionality check: if there is a mismatch between
            % provided center frequencies and data dimension, try to
            % transpose the data to fix the problem. If impossible, the
            % provided center frequencies are discarded, data is not
            % transposed, and a warning is issued.
            s = size(data);             % Data dimensions
            n_f = size(cfHz(:),1);      %#ok<PROP> Number of audio frequency bins 
            n_fm = size(modCfHz(:),1);  %#ok<PROP> Number of modulation freq. bins
            
            if size(s,2)==2
                % Add a third dimension to avoid generating an
                % error
                s = [s 0];
            end
            
            % Check for a mismatch
            if ismember(n_f*n_fm,s) && s(3)==0
                % Then the provided data contained interleaved channels,
                % need to sort them out
                
                if s(1) == n_f*n_fm
                    % Transpose the data to have time first
                    data = data.';
                end
                
                % Reshape the data
                data = permute(reshape(data,[size(data,1) n_fm n_f]),[1 3 2]);
                
            elseif ((n_f == s(2))&&(n_fm == s(3)))||(max(s)==0)
                % Then everything should be fine, carry on
                
            elseif (~ismember(n_f,s)||~ismember(n_fm,s))&&(max(s)~=0)
                % Then provided center frequencies do not match the data
                
                warning('Provided vectors of center frequencies do not match the data, ignoring them')
                cfHz = [];      %#ok<PROP>
                modCfHz = [];   %#ok<PROP>
                
            else
                % Dimensions matches, but a permutation is required
                
                % Find the number of time bins
                n_t = s(s~=n_f);
                n_t = n_t(n_t~=n_fm);
                
                % Permutation vector
                p = [find(n_t==s) find(n_f==s) find(n_fm==s)];
                
                % Request the permutation
                data = permute(data,p);
                
                % Issue a warning
                warning('Data was permuted to match the provided vectors of center frequencies')
                
            end
            
            % NB: This functionality is implemented to allow defining such
            % signals from a 2D matrix with interleaved audio and
            % modulation frequencies, but does not (cannot) check that
            % proper input was provided in this case!
            
            % Populate object properties
            sObj.Dimensions = 'nSample x nFilters x nModulationFilters';
            sObj.cfHz = cfHz(:)';       %#ok<PROP>
            sObj.setData(data);
            sObj.Channel = channel;
            sObj.modCfHz = modCfHz(:)'; %#ok<PROP>
                
            end
            
            
        end
        
        function h = plot(sObj,h0,p)
            %plot   Plot a modulation signal
            %
            %USAGE:
            %      sObj.plot
            %      sObj.plot(...)
            %  h = sObj.plot(...)
            %
            %INPUT ARGUMENTS
            %  sObj : Signal instance
            %
            %OUTPUT ARGUMENTS
            %     h : Handle to the figure
            %
            %TODO: - Clean up and use plot properties. 
            %      - Add other plot options (given audio frequency?)
            
            % Manage plotting parameters
            if nargin < 3 || isempty(p) 
                % Get default plotting parameters
                p = Parameters.getPlottingParameters('TimeFrequencySignal');
            else
                defaultPar = Parameters.getPlottingParameters('TimeFrequencySignal');
                defaultPar.replaceParameters(p);
                p = defaultPar;
            end
            
            % Extract the data from the buffer
            data = sObj.Data(:);
            s = size(data);
            
            % Limit dynamic range of AMS representation
            maxDynamicRangedB = p.map('dynrange');

            % Reshape data to incorporate borders
            rsAMS = permute(data,[3 2 1]);
            rsAMS = reshape(20*log10(abs(rsAMS)),[s(3) s(2) s(1)]);
            
            maxValdB = max(rsAMS(:));
            
            % Minimum ratemap floor to limit dynamic range
            minValdB = -(maxDynamicRangedB + (0 - maxValdB));

            rangedB = [round(minValdB/5)*5 round(maxValdB/5)*5];

            rsAMS(end+1,:,:) = NaN;  % insert borders
            
            % Reshape for plotting (interleaved audio and modulation
            % frequencies)
            rsAMS = reshape(rsAMS,[(s(3)+1)*s(2) s(1)]);
            
            % Generate a time vector
            timeSec = 0:1/sObj.FsHz:(s(1)-1)/sObj.FsHz;
                
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

            im = imagesc(timeSec,1:s(2)*(1+s(3)),rsAMS);
            
            % Make borders appear black
            set(gca,'Color','k')                % Black background
            set(im,'AlphaData',~isnan(rsAMS))   % Transparent borders
            
            % Set color map
            try
                colormap(p.map('colormap'))
            catch
                warning('No colormap %s is available, using ''jet''.',p.map('colormap'))
                colormap('jet')
            end
            
            if p.map('bColorbar')
                colorbar;
            end

            set(gca,'CLim',rangedB)
            
            axis xy
            
            if ~strcmp(sObj.Channel,'mono')
                pTitle = [sObj.Label ' - ' sObj.Channel];
            else
                pTitle = sObj.Label;
            end
            
            % Set up title and labels
            title(pTitle,'fontsize',p.map('fsize_title'),'fontname',p.map('ftype'))
            xlabel('Time (s)','fontsize',p.map('fsize_label'),'fontname',p.map('ftype'))
            ylabel('Center frequency (Hz)','fontsize',p.map('fsize_label'),...
                   'fontname',p.map('ftype'))
            
            nYLabels = 5;
            
            yPosInt = round(linspace(1,s(2),nYLabels));
            
            % Center yPos at mod filter frequencies
            yPos = (yPosInt - 1) * (s(3)+1) + round((s(3)+1)/2);
            
            % Find the spacing for the y-axis which evenly divides the y-axis
            set(gca,'ytick',yPos);
            set(gca,'yticklabel',round(sObj.cfHz(yPosInt)));

            for ii = 1:s(2)
                text([timeSec(3) timeSec(3)], ...
                    floor(s(3)/2)+[((ii-1)*(s(3)+1)) ((ii-1)*(s(3)+1))], ...
                    num2str(ii),'verticalalignment','middle','fontsize',8);
            end
            
            % 
            
        end
        
    end
    
    
end