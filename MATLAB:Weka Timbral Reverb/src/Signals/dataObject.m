classdef dataObject < dynamicprops
%DATAOBJECT: Signal container class for the auditory front-end (AFE) framework.
%   A data object is necessary for running the AFE framework. It contains as individual
%   properties all signals: the input signal, the output signal(s) requested by the user, 
%   but also all the intermediary representations necessary to compute the request.
%   All signals are individual objects inheriting the Signal class.
%
%   DATAOBJECT properties:
%       bufferSize_s - Global buffer size for all signals in seconds.
%       isStereo     - Flag indicating a binaural signal.
%       'signalname' - All the contained signals are defined as dynamic properties. The
%                      property name is taken from the signal object name property.
%                      Multiple signals with same names are arranged in a cell array, with
%                      columns spanning left/right channels.
%
%   DATAOBJECT methods:
%       dataObject          - Constructor for the class.
%       addsignal           - Adds a signal to the data object.
%       clearData           - Clears buffered data in all contained signals.
%       getParameterSummary - Returns the parameters used for computing all signals.
%       play                - Plays back the audio from the input signal.
%       plot                - Plots the input signal waveform.
%                      
%
% See also Signal, signals (folder)
    
    
    properties
        % bufferSize_s - Global buffer size for all signals (seconds). Due to the
        % compatibility with chunk-based processing, all signals need to be buffered to a
        % certain duration as the duration of the input signal is unknown and assumed
        % infinite.
        % See also circVBuf, circVBufArrayInterface
        bufferSize_s;
        % isStereo - Flag indicating if the data structure will be based on a stereo 
        % (binaural) signal. This is necessary in online scenario, when an empty
        % dataObject is initialized and there is no explicit information about the number
        % of channels in the input signal.
        isStereo;
    end
    
    
    methods
        function dObj = dataObject(s,fs,bufferSize_s,channelNumber)
            %dataObject     Constructs a data object from an (optional) input signal
            %
            %USAGE
            %       dObj = dataObject(s,fs)
            %       dObj = dataObject(s,fs,bufferSize)
            %       dObj = dataObject([],fs,[],channelNb)
            %
            %INPUT ARGUMENTS
            %          s : Initial time-domain signal
            %         fs : Sampling frequency
            % bufferSize : length of the signal buffer in seconds (default = 10)
            %  channelNb : Number of channels, 1 (mono-default) or 2 (stereo) in the 
            %              original signal. Used in chunk-based scenario when the 
            %              dataObject is initialized with an empty signal (third example)
            %
            %OUTPUT ARGUMENTS
            %       dObj : Data object
            
            if nargin==1
                error('The sampling frequency needs to be provided') 
            end
            if nargin==0
                s = [];
                fs = [];
            end
            if nargin < 3
                bufferSize_s = 10;
            end
            dObj.bufferSize_s = bufferSize_s;
            
            % Check number of channels of the provided signal
            if size(s,2)==2 
                % Set to stereo when provided signal is
                channelNumber = 2;
            elseif nargin<4||isempty(channelNumber)
                % Set to default if channel number not provided
                channelNumber = 1;
            end
            
            if channelNumber>2 || channelNumber<1
                error('Provided number of channel is incorrect, should be 1 (mono) or 2 (stereo).')
            end
            
            % Set the is_stereo property
            dObj.isStereo = channelNumber-1;
            
            % Populate the signal property
            if dObj.isStereo
                if ~isempty(s)
                    sig_l = TimeDomainSignal.construct(fs, dObj.bufferSize_s, ...
                            'input', 'Ear Signal', 'left', s(:,1));
                    sig_r = TimeDomainSignal.construct(fs, dObj.bufferSize_s, ...
                            'input', 'Ear Signal', 'right', s(:,2));
                else
                    sig_l = TimeDomainSignal.construct(fs, dObj.bufferSize_s, ...
                            'input', 'Ear Signal', 'left', []);
                    sig_r = TimeDomainSignal.construct(fs, dObj.bufferSize_s, ...
                            'input', 'Ear Signal', 'right', []);
                end
                dObj.addSignal(sig_l);
                dObj.addSignal(sig_r);
            else
                if ~isempty(s)
%                     sig = TimeDomainSignal(fs,dObj.bufferSize_s,'input','Ear signal (mono)',s);
                    sig = TimeDomainSignal.construct(fs, dObj.bufferSize_s, ...
                            'input', 'Ear signal', 'mono', s);
                else
                    sig = TimeDomainSignal.construct(fs, dObj.bufferSize_s, ...
                            'input', 'Ear signal', 'mono', []);
                end
                dObj.addSignal(sig);
            end          
        end
        
        function addSignal(dObj,sObj)
            %addSignal  Incorporates an additional signal object to a data object
            %
            %USAGE
            %     dObj.addSignal(sObj)
            %
            %INPUT ARGUMENTS
            %      dObj : Data object to add the signal to
            %      sObj : Signal object to add
            %
            %N.B. This method uses dynamic property names. The data object dObj will
            %     contain the signal sObj as a new property, named after sObj.Name. If
            %     such a property existed beforehand, the new signal will be incorporated
            %     in a cell array under that property name.
            
            % TODO: Left/right channel handling works here only because
            % left channel is always instantiated first. Might want to
            % check if that is a limitation.
            
            % Which channel (left, right, mono) is this signal
            if strcmp(sObj.Channel,'right')
                jj = 2;
            else
                jj = 1;
            end
            
            % Check if a signal with this name already exist
            if isprop(dObj,sObj.Name)
                ii = size(dObj.(sObj.Name),1)+2-jj;
                dObj.(sObj.Name){ii,jj} = sObj;
            else
                dObj.addprop(sObj.Name);
                dObj.(sObj.Name) = {sObj};
            end
            
        end
        
        function clearData(dObj,bClearSignal)
            %clearData  Clears data of all signals in the data structure
            %
            %USAGE:
            %   dObj.clearData
            %   
            %N.B. Use dObj.clearData(0) to clear all signals BUT the input signal.
            
            if nargin<2 || isempty(bClearSignal)
                bClearSignal = 1;
            end
            
            % Get a list of the signals contained in the data object
            sig_list = fieldnames(dObj);
            
            % Remove the "isStereo" and "bufferSize_s" properties from the list
            sig_list = setdiff(sig_list,{'isStereo' 'bufferSize_s'});
            
            % Remove the signal from the list if needed
            if ~bClearSignal
                sig_list = setdiff(sig_list,{'input'});
            end
                
            % Loop over all the signals
            for ii = 1:size(sig_list,1)
                
                % There should always be a left or mono channel
                dObj.(sig_list{ii}){1}.clearData;
                
                % Check if there is a right channels
                if size(dObj.(sig_list{ii}),2)>1
                    
                    % It could still be empty (e.g. for "mix" signals)
                    if isa(dObj.(sig_list{ii}){2},'Signal')
                        dObj.(sig_list{ii}){2}.clearData;
                    end
                    
                end
                
            end
           
            
        end
        
        function p = getParameterSummary(dObj,mObj)
            %getParameterSummary  Returns a structure parameters used for computing each 
            %                     signal in the data object.
            %
            %USAGE:
            %   p = dObj.getParameterSummary(mObj)
            %
            %INPUT ARGUMENTS: 
            %   dObj : Data object instance
            %   mObj : Manager instance associated to the data
            %
            %OUTPUT ARGUMENTS:
            %      p : Structure of used parameter values
            
            % TODO: Update to refactored code version
            
            % Get a list of instantiated signals
            prop_list = properties(dObj);
            sig_list = setdiff(prop_list,{'isStereo' 'bufferSize_s'});
            
            % Initialize the output
            p = struct;
            
            % Loop on each signal
            for ii = 1:size(sig_list,1)
                
                % Test if multiple representations exist 
                if size(dObj.(sig_list{ii}),1)>1
                    % There are multiple representations with this name
                    
                    % Use a cell array
                    p.(sig_list{ii}) = cell(size(dObj.(sig_list{ii}),1),1);
                    
                    % Get the parameters
                    for jj = 1:size(dObj.(sig_list{ii}),1)
                        p.(sig_list{ii}){jj} = dObj.(sig_list{ii}){jj,1}.getParameters(mObj);
                    end
                    
                else
                    % There is only one such representation
                    p.(sig_list{ii}) = dObj.(sig_list{ii}){1,1}.getParameters(mObj);
                end
                    
                
            end
            
            
        end
        
        function play(dObj,bPreProcessed)
            %play   Playback the audio from the input signal contained in the data object
            %
            %USAGE
            %   dObj.play
            %
            %INPUT ARGUMENTS
            %   dObj : Data object
            %
            %
            %N.B. Use dObj.play(1) to play the pre-processed input signal, if available
            
            if nargin<2 || isempty(bPreProcessed) || ~bPreProcessed
                if ~isprop(dObj,'input')||isempty(dObj.input)||...
                        isempty(dObj.input{1}.Data)
                    warning('There is no audio in the data object to playback')
                else
                    if size(dObj.input,2)==1
                        % Then mono playback
                        soundsc(dObj.input{1}.Data(:),dObj.input{1}.FsHz)
                    else
                        % Stereo playback
                        temp_snd = [dObj.input{1}.Data(:) dObj.input{2}.Data(:)];
                        soundsc(temp_snd,dObj.input{1}.FsHz)
                    end
                end
            else
                if ~isprop(dObj,'time')||isempty(dObj.time)||...
                        isempty(dObj.time{1}.Data)
                    warning('There is no audio in the data object to playback')
                else
                    if size(dObj.time,2)==1
                        % Then mono playback
                        soundsc(dObj.time{1}.Data(:),dObj.time{1}.FsHz)
                    else
                        % Stereo playback
                        temp_snd = [dObj.time{1}.Data(:) dObj.time{2}.Data(:)];
                        soundsc(temp_snd,dObj.time{1}.FsHz)
                    end
                end
            end
        end
        
        function h = plot(dObj,h0,p,varargin)
            %plot       Plots the time-domain waveforms of the signal after pre-processing
            %
            %USAGE
            %       dObj.plot
            %       dObj.plot(h_prev,p,...)
            %       h = dObj.plot(...)
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
            %OPTIONAL PARAMETERS - (key,value) pair
            % 'bGray'         - Set to 1 for gray shades plot
            % 'rangeSec'      - Vector of time limits for the plot 
            % 'bSignal'       - Set to 1 to plot the input signal before pre-processing
            % 'decimateRatio' - Integer ratio to decimate the plot (for smaller files)
            
            % Manage plotting parameters
            if nargin < 3 || isempty(p) 
                % Get default plotting parameters
                p = Parameters.getPlottingParameters();
            else
                defaultPar = Parameters.getPlottingParameters();
                defaultPar.replaceParameters(p);
                p = defaultPar;
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
            
            % Manage optional arguments
            if nargin>3 && ~isempty(varargin)
                opt = struct;
                for ii = 1:2:size(varargin,2)
                    opt.(varargin{ii}) = varargin{ii+1};
                end
            else
                opt = [];
            end
            
            % Colors
            if ~isempty(opt) && isfield(opt,'bGray')
                if opt.bGray
                    colors = {[0 0 0] [.5 .5 .5]};
                else
                    colors = p.map('colors');
                end
            else
                colors = p.map('colors');
            end
            
            % Plot before/after pre-processing
            if ~isempty(opt) && isfield(opt,'bSignal')
                if opt.bSignal
                    sig = dObj.input;
                else
                    sig = dObj.time;
                end                
            else
                sig = dObj.time;
            end
            
            % Decimate the waveforms for "lighter" plots
            if ~isempty(opt) && isfield(opt,'decimateRatio')
                dcR = opt.decimateRatio;
            else
                dcR = 1;
            end
            
            % Limit to a certain time period
            if ~isempty(opt) && isfield(opt,'rangeSec')
                data = sig{1}.Data(floor(opt.rangeSec(1)*sig{1}.FsHz):dcR:floor(opt.rangeSec(end)*sig{1}.FsHz));
                if size(sig,2)>1
                    data = [data sig{2}.Data(floor(opt.rangeSec(1)*sig{1}.FsHz):dcR:floor(opt.rangeSec(end)*sig{1}.FsHz))];
                end
                t = opt.rangeSec(1):dcR/sig{1}.FsHz:opt.rangeSec(1)+dcR*(size(data,1)-1)/sig{1}.FsHz;
            else
                data = sig{1}.Data(1:dcR:end);
                if size(sig,2)>1
                    data = [data sig{2}.Data(1:dcR:end)];
                end
                t = 0:dcR/sig{1}.FsHz:dcR*(size(data,1)-1)/sig{1}.FsHz;
            end
            
            % Time vector
%             t = 0:1/sig{1}.FsHz:(size(sig{1}.Data(:),1)-1)/sig{1}.FsHz;
            
            % Plot channel with highest energy (or mono channel) first
            if size(sig,2)==1 || norm(data(:,1),2) > norm(data(:,2),2)
                plot(t,data(:,1),'linewidth',p.map('linewidth_s'),'color',colors{1});
                
                if size(sig,2)>1
                    hold on
                    plot(t,data(:,2),'linewidth',p.map('linewidth_s'),'color',colors{2});
                    legend([sig{1}.Channel ' ear'],[sig{2}.Channel ' ear'],'location','NorthEast')
                end
            else
                plot(t,data(:,2),'linewidth',p.map('linewidth_s'),'color',colors{1});
                hold on
                plot(t,data(:,1),'linewidth',p.map('linewidth_s'),'color',colors{2});
                legend([sig{2}.Channel ' ear'],[sig{1}.Channel ' ear'],'location','NorthEast')
            end
            xlabel('Time (s)','fontsize',p.map('fsize_axes'))
            ylabel('Amplitude','fontsize',p.map('fsize_axes'))
            xlim([t(1) t(end)])
            title('Time domain signals','fontsize',p.map('fsize_title'))
            set(gca,'fontname',p.map('ftype'),'fontsize',p.map('fsize_label'))
            
        end
        
    end
    
    methods (Hidden = true)
       
        function replaceInputSignal(dObj,s)
            %REPLACEINPUTSIGNAL    Clears out all data and appends a new input signal
            %
            %USAGE:
            %  dObj.replaceInputSignal(s)
            %
            %INPUT ARGUMENTS:
            %  dObj : Data object instance
            %     s : New input signal
            
            if size(s,2)-1 ~= dObj.isStereo
                warning(['Cannot replace input signal: number of channels in data ' ...
                    'object and new input signal do not match. Consider instanciating ' ...
                    'a new data object.'])
                return
            end
            
            % Clear all data
            dObj.clearData;
            
            % Append new input
            dObj.input{1}.appendChunk(s(:,1));
            
            if size(s,2) == 2
                % Stereo signal
                dObj.input{2}.appendChunk(s(:,2));
            end
            
        end
        
    end
    
    
    
end
