classdef ccFeatureProc < Processor
    
    properties
        factor      % Factor by which the cross-correlation is downsampled
    end
    
    methods
    
        function pObj = ccFeatureProc(fs,factor)
            %ccFeatureProc      Instantiate a cross-correlation feature
            %processor
            %
            %USAGE:
            %   pObj = ccFeatureProc(fs,factor)
            %
            %INPUT ARGUMENTS
            %     fs : Sampling frequency
            % factor : Downsampling factor 
            %
            %OUTPUT ARGUMENTS
            %   pObj : Processor instance
            
            if nargin>0     % Safeguard for Matlab empty calls
                
            if nargin<2
                error('Not enough input arguments')
            end
            
            if mod(factor,1)~=0||factor<=0
                error('Provided factor should be an positive integer')
            end
            
            % Populate processor properties
            pObj.Type = 'Cross-correlation feature';
            pObj.FsHzIn = fs;
            pObj.FsHzOut = fs;
            pObj.factor = factor;

            end

            % Hide the processor (for the moment)
            pObj.bHidden = 1;
            
        end
        
        function out = processChunk(pObj,in)
            %processChunk   Calls the processing for a new chunk of signal
            %
            %USAGE
            %   out = pObj.processChunk(in)
            %
            %INPUT ARGUMENTS
            %  pObj : Processor instance
            %    in : Input signal
            %
            %OUTPUT ARGUMENTS
            %   out : Output signal
            
            % Get the dimensions of the input
            dim = size(in);
            
            % We need to downsample the 3rd dimension of the input (lags)
            % but still include the origin (zero-lag) in the output
            
            
            % 1- Find the index of the origin. dim(3) is the number of lags
            % and should be an odd integer. The origin is the middle point
            % of this vector.
            origin = (dim(3)+1)/2;
            
            % 2- Figure out a vector of index of the points to retain, that
            % should include the origin
            
            % New number of lags on each side of the origin
            n_lags = floor((origin-1)/pObj.factor); 
            
            % Vector of index to retain:
            index = origin-n_lags*pObj.factor:pObj.factor:origin+n_lags*pObj.factor;
            
            % Downsample the input
            out = in(:,:,index);
            
        end
        
        function reset(pObj)
            % Nothing to reset for this processor, just implementing an
            % abstract class
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
            
            if pObj.factor == p.ccf_factor
                hp = 1;
            else
                hp = 0;
            end
            
        end
        
    end
    
end