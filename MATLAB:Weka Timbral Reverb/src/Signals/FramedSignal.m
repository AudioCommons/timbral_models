classdef FramedSignal < Signal
    
    properties
        frameFsHz      % Sampling frequency inside  
    end
    
    methods
    
        function sObj = FramedSignal(fs,bufferSize_s,frameSize,frameFs,name,label,channel)
            %FramedSignal   Constructor for the framed signal class
            %
            %USAGE:
            %   sObj = FramedSignal(fs,frameSize,frameFs)
            %   sObj = FramedSignal(fs,frameSize,frameFs,name,label,channel)
            %
            %INPUT ARGUMENTS
            %        fs : Sampling frequency (inverse of frame step-size)
            % frameSize : Frame size in samples
            %   frameFs : Sampling frequency inside a frame (Hz)
            %      name : Name tag of the signal
            %     label : Label for the signal
            %   channel : 'left', 'right', or 'mono' (default)
            %
            %OUTPUT ARGUMENTS
            %    sObj : Signal instance
            
            sObj = sObj@Signal(fs,bufferSize_s,frameSize);
            
            if nargin>0
                
            if nargin<7||isempty(channel);channel='mono';end
            if nargin<5||isempty(name);name='framedSignal';end
            if nargin<6||isempty(label);label=name;end
            
            if nargin<4
                error('Sampling frequencies are needed to instantiate a framed signal')
            end
            
            % Populate signal properties
            populateProperties(sObj,'Label',label,'Name',name,...
                'Dimensions','nFrames x frameSize');
            sObj.frameFsHz = frameFs;
            sObj.Channel = channel;
                
            end
            
        end
        
        function  h = plot(sObj,h0)
            % TO DO: Implement (if that is ever needed)
            h = [];
            warning('Framed signal have no automated plotting routine yet.')
        end
        
    end
    
end