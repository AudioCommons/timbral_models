classdef TimeDomainSignal < Signal
%TIMEDOMAINSIGNAL Signal class for simple, one dimensional, time-domain signals.
%
% See also Signal, preProc
    
    properties
        % Only inherited properties
    end
    
    methods
        function sObj = TimeDomainSignal(procHandle,bufferSize,channel,data)
            %TimeDomainSignal       Constructor for the "time domain signal"
            %                       children signal class
            %
            %USAGE
            %     sObj = TimeDomainSignal(procHandle)
            %     sObj = TimeDomainSignal(procHandle,bufferSize,channel,data)
            %
            %INPUT ARGUMENTS
            % procHandle : Handle to the processor generating this signal as output
            % bufferSize : Size of the ring buffer in s (default: bufferSize = 10)
            %    channel : Flag indicating 'left', 'right', or 'mono' (default: 
            %              channel = 'mono')
            %       data : Vector of amplitudes to construct an object from existing data
            %
            %OUTPUT ARGUMENT
            %  sObj : Time domain signal object inheriting the signal class
            
            % Check input arguments
            if nargin<4; data = []; end
            if nargin<3||isempty(channel); channel = 'mono'; end
            if nargin<2||isempty(bufferSize); bufferSize = 10; end
            if nargin<1||isempty(procHandle); procHandle = emptyProc; end
            
            sObj = sObj@Signal( procHandle, bufferSize, 1 );
            
            if nargin>0  % Failproof for Matlab empty calls
            
            % Check dimensionality of data if it was provided
            if ~isempty(data) && min(size(data))>1
                error(['The data used to instantiate this object should be a' ...
                    'single vector of amplitude values'])
            end

            % Format data to a column vector
            data = data(:);
            
            % Populate object properties
            sObj.Dimensions = 'nSamples x 1';
            sObj.setData( data );
            sObj.Channel = channel;
            
            end
        end
       
        function h = plot(sObj,h0,p,varargin)
            %plot       This method plots the data from a time domain
            %           signal object
            %
            %USAGE
            %       sObj.plot
            %       sObj.plot(h_prev,p)
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
                
                % Manage plot parameters
                if nargin < 3 || isempty(p) 
                    % Get default plotting parameters
                    p = Parameters.getPlottingParameters('TimeDomainSignal');
                else
%                     p.fs = sObj.FsHz;   % Add the sampling frequency to satisfy parseParameters
%                     p = parseParameters(p);
                    defaultPar = Parameters.getPlottingParameters('TimeDomainSignal');
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
                
                % Limit to a certain time period
                if ~isempty(opt) && isfield(opt,'rangeSec')
                    data = sObj.Data(floor(opt.rangeSec(1)*sObj.FsHz):floor(opt.rangeSec(end)*sObj.FsHz));
                    t = opt.rangeSec(1):1/sObj.FsHz:opt.rangeSec(1)+(size(data,1)-1)/sObj.FsHz;
                else
                    data = sObj.Data(:);
                    t = 0:1/sObj.FsHz:(size(data,1)-1)/sObj.FsHz;
                end
                

                % Set up a title (include channel in the title)
                if ~strcmp(sObj.Channel,'mono')
                    pTitle = [sObj.Label ' - ' sObj.Channel];
                else
                    pTitle = sObj.Label;
                end
                
                % Plot
                plot(t,data,'color',p.map('color'),'linewidth',p.map('linewidth_s'))
                xlabel('Time (s)','fontsize',p.map('fsize_label'),'fontname',p.map('ftype'))
                ylabel('Amplitude','fontsize',p.map('fsize_label'),'fontname',p.map('ftype'))
                title(pTitle,'fontsize',p.map('fsize_title'),'fontname',p.map('ftype'))
                set(gca,'fontsize',p.map('fsize_axes'),'fontname',p.map('ftype'))
                
                box on 
                
                % Center the waveform
                m = max(abs(data));
                if m == 0, m = 1; end;
                set(gca,'XLim',[t(1) t(end)],'YLim',[-1.1*m 1.1*m])
            
            else
                warning('This is an empty signal, cannot be plotted')
            end
                
            
        end
        
        function play(sObj)
            %play       Playback the audio from a time domain signal
            %
            %USAGE
            %   sObj.play()
            %
            %INPUT ARGUMENTS
            %   sObj : Time domain signal object
            
            sound(sObj.Data(:),sObj.FsHz)
            
        end
        
    end
    
    methods (Static)
        
        function sObj = construct(fs,bufferSize,name,label,channel,data)
            %construct    Constructs a time domain signal given its properties.
            %             To be used to generate a standalone signal, i.e., one that is
            %             not given as the output of a processor.
            %
            %USAGE:
            %  sObj = TimeDomainSignal.construct(fs)
            %  sObj = TimeDomainSignal.construct(fs, bufferSize, name,  label, ...
            %                                    channel, data)
            %
            %INPUT ARGUMENTS:
            %         fs : Sampling frequency (Hz)
            % bufferSize : Ring buffer duration in seconds
            %       name : Name for the signal
            %      label : Label for the signal (used in e.g., plotting)
            %    channel : Channel tag ('left', 'right', or 'mono')
            %       data : Initial data to fill the buffer with
            %
            % NB: If arguments are missing, default values from the class constructor are
            %     used
            %
            %OUTPUT ARGUMENT:
            %   sObj : Time domain signal instance
            
            if nargin<6; data = []; end
            if nargin<5; channel = []; end
            if nargin<4; label = []; end
            if nargin<3; name = []; end
            if nargin<2||isempty(bufferSize); bufferSize = 10; end
            if nargin<1; fs = []; end
            
            % Create a structure that contains all the info needed with correct formating
            dummyStruct = struct;
            dummyStruct.FsHzOut = fs;
            dummyStruct.getProcessorInfo.requestName = name;
            dummyStruct.getProcessorInfo.requestLabel = label;
            
            % Instantiate the signal
            sObj = TimeDomainSignal(dummyStruct,bufferSize,channel,data);
            
        end
        
        function [names, defaultValues, descriptions] = getPlottingParameterInfo()
            %GETPLOTTINGPARAMETERINFO   Stores plot parameters that are common to all
            %signals.
            
            
            names = {};
                 
            descriptions = {};
                
            defaultValues = {};
            
        end
        
        
    end
    
end