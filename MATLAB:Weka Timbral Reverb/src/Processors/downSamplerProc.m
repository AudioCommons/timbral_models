classdef downSamplerProc < Processor
    
    properties
        WorkingDim      % Index of the dimension to be downsampled
        Ratio           % Downsampling ratio (before/after)
    end
    
    properties (GetAccess = private)
        bStart          % Starting index for next block (acts as a buffer)
    end
    
    methods
        
        function pObj = downSamplerProc(fs,ratio,dim)
            %downSampler    Instantiate a downsampling processor
            %
            %USAGE:
            %  pObj = downSampler(fs,ratio)
            %  pObj = downSampler(fs,ratio,dim)
            %
            %INPUT PARAMETERS
            %    fs : Sampling frequency
            % ratio : Downsampling ratio
            %   dim : Dimension to work upon (default: 1)
            %
            %OUTPUT PARAMETERS
            %  pObj : Processor instance
            
            % Hide this processor from the framework
            pObj.bHidden = 1;
            
            if nargin>0     % Safeguard for Matlab empty calls
                
            % Checking input parameters
            if nargin<3||isempty(dim)
                dim = 1;
            end
            
            if nargin<2
                error(['Downsampler processor needs sampling frequency '...
                    'and downsampling ratio to be instantiated.'])
            end
            
            % The down-sampling ratio should be a positive integer
            if mod(ratio,1)~=0
                error('Downsampling ratio should be a positive integer')
            end
            
            % Populate properties
            pObj.Type = 'Down-sampler';
            pObj.FsHzIn = fs;
            if dim == 1
                % Then time is down-sampled
                pObj.FsHzOut = fs/ratio;
            else
                % Then another dimension is down-sampled. Sampling
                % frequency is unchanged
                pObj.FsHzOut = fs;
            end
            pObj.WorkingDim = dim;
            pObj.Ratio = ratio;
                
            % Initialize block shift:
            % This processor uses the downsample.m script from Matlab. It
            % simply picks every N sample, where N is the downsampling
            % ratio, starting always with the first sample. Hence there is
            % no need to store a buffer, but only the index (or shift of
            % index) at which the next block should pick the first sample
            % of the downsampled signal, to be consistent with the previous
            % block.
            pObj.bStart = 1;   
            
            end
            
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
            %   out : Output (down-sampled) signal
            
            % Append provided input to the buffer
%             if ~isempty(pObj.buffer)
%                 
%                 % This operation is dimension-dependent!
%                 in = cat(pObj.WorkingDim,pObj.buffer,in);
%                 
%             end
            
            % Do nothing and return a warning if input has fewer
            % dimensions than pObj.WorkingDim
            if size(size(in),2)<pObj.WorkingDim
                warning(['Working dimension of the downsampler is not a'...
                    ' valid dimension of the input']) 
            else
            
            % Get buffered input dimensions
            dim = size(in);
            
            % We want to move the dimension along which the downsampling
            % occurs to be the first dimension if needed:
            if pObj.WorkingDim~=1
                %  - Generate a permutation order vector
                order = 1:size(dim,2);
                order(1) = pObj.WorkingDim;
                order(pObj.WorkingDim) = 1;

                %  - Permute the dimensions
                in = permute(in,order);
                
                %  - Update the dimension vector
                dim = size(in);
            end
            
            % Matlab's downsample.m function can operate on max 2
            % dimensions. Need to reshape the input if it is 3D or more
            if size(dim,2)>2
                % Then additional dimensions are squeezed into one
                in = reshape(in,[dim(1) prod(dim(2:end))]);
                
                % Perform resampling colum-wise, starting from the block
                % index
                out = downsample(in(pObj.bStart:end,:),pObj.Ratio);
                
                % Re-organize data
                out = reshape(out,[size(out,1) dim(2:end)]);
            else
                if size(dim,2)==1
                    out = downsample(in(pObj.bStart:end),pObj.Ratio);
                else
                    % Perform resampling colum-wise
                    out = downsample(in(pObj.bStart:end,:),pObj.Ratio);
                end
            end
            
            % Reorganize the dimensions of the output if they were modified
            if pObj.WorkingDim~=1
                out = permute(out,order);
            end
            
            % Length of the actual input
            L = size(in,pObj.WorkingDim)-pObj.bStart+1;
            
            % Update the starting index for next block
            pObj.bStart = mod(pObj.Ratio-rem(L,pObj.Ratio),pObj.Ratio)+1;
            
            end
            
        end
            
        function reset(pObj)
            %reset     Resets the internal states of the processor
             %
             %USAGE
             %      pObj.reset
             %
             %INPUT ARGUMENTS
             %  pObj : Down-sampler processor instance
             
             % Empty the buffer
             pObj.buffer = [];
            
        end
        
        
    end
   
    
    
end