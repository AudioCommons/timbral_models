classdef framingProc < Processor
    
    properties
        wname       % Window shape descriptor (see window.m)
        wSize       % Window duration in samples
        hSize       % Step size between windows in samples
    end
    
    properties (GetAccess = private)
        win
        buffer
    end
    
    methods
        
        function pObj = framingProc(fs,wname,wSize,hSize)
            %framingProc        Constructs a framing processor
            %
            %USAGE:
            %   pObj = framingProc(wname,wSize,hSize)
            %
            %INPUT PARAMETERS:
            %  wname : Window shape descripto (see window.m)
            %  wSize : Window size in samples
            %  hSize : Step size between windows in samples
            %
            %OUTPUT PARAMETERS:
            %   pObj : Processor instance
            
            if nargin>0     % Failsafe for Matlab empty calls
                
            if nargin<3
                error('Not enough input arguments')
            end
            
            pObj.Type = 'Frame extractor';
            pObj.FsHzIn = fs;
            pObj.FsHzOut = fs/hSize;
            
            pObj.wname = wname;
            pObj.wSize = wSize;
            pObj.hSize = hSize;
            
            pObj.win = window(wname,wSize);
            pObj.buffer = [];

            end

            % Hide the processor (for the moment)
            pObj.bHidden = 1;
        
        end
        
        function out = processChunk(pObj,in)
            %processChunk       Apply the processor to a new chunk of input
            %                   signal
            %
            %USAGE
            %   out = pObj.processChunk(in)
            %
            %INPUT ARGUMENT
            %    in : New chunk of input data
            %
            %OUTPUT ARGUMENT
            %   out : Corresponding output
            
            % Input should be a vector
            if size(in,2)>1
                error('Input should be a column vector')
            end
            
            % Append input to previous buffer
            if ~isempty(pObj.buffer)
                in = cat(1,pObj.buffer,in);
            end
            
            % Frame the data
            out = frameData(in,pObj.wSize,pObj.hSize,pObj.win,false).';
            
            % Update the buffer
            nFrames = size(out,1);      % Number of frames extracted
            pObj.buffer = in(nFrames*pObj.hSize+1:end);
              
        end
        
        function reset(pObj)
            %reset     Resets the internal states of the ratemap extractor
            %
            %USAGE
            %      pObj.reset
            %
            %INPUT ARGUMENTS
            %  pObj : Processor instance
            
            % Empty the buffer
            pObj.buffer = [];
            
        end
        
        function hp = hasParameters(pObj,p)
            %hasParameters  This method compares the parameters of the
            %               processor with the parameters given as input
            %
            %USAGE
            %    hp = pObj.hasParameters(p)
            %
            %INPUT ARGUMENTS
            %  pObj : Processor instance
            %     p : Structure containing parameters to test
            
            p_list_proc = {'wname','wSize','hSize'};
            p_list_par = {'fr_wname','fr_wSize','fr_hSize'};
            
            % Initialization of a parameters difference vector
            delta = zeros(size(p_list_proc,2),1);
            
            % Loop on the list of parameters
            for ii = 1:size(p_list_proc,2)
                try
                    if ischar(pObj.(p_list_proc{ii}))
                        delta(ii) = ~strcmp(pObj.(p_list_proc{ii}),p.(p_list_par{ii}));
                    else
                        delta(ii) = abs(pObj.(p_list_proc{ii}) - p.(p_list_par{ii}));
                    end
                    
                catch err
                    % Warning: something is missing
                    warning('Parameter %s is missing in input p.',p_list_par{ii})
                    delta(ii) = 1;
                end
            end
            
            % Check if delta is a vector of zeros
            if max(delta)>0
                hp = false;
            else
                hp = true;
            end
            
            
            
        end
    end
    
end