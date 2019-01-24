classdef onsetMapProc < Processor
%TRANSIENTMAPPROC Binary onset and offset maps processor.    
%   Based on the transient strength which is derived from the corresponding 
%   onset strength and offset strength processor, a binary decision about 
%   transient activity is formed, where only the most salient information
%   is retained. This can be used to group the acoustic input according to 
%   individual auditory events [1].
%
%   TRANSIENTMAPPROC properties:
%        minStrengthdB   - Minimum transient strength for mapping
%        minSpread       - Minimum spread of the transient (number of time-frequency units)
%        fuseWithinSec   - Time constant below which transients are fused
%        minValuedB      - Lower limit for the input representation below which transient are discarded
%
%   See also: Processor, ratemapProc, onsetProc, offsetProc
%
%   Reference:
%   [1] Turgeon, M., Bregman, A. S., and Ahad, P. A. (2002), "Rhythmic 
%       masking release: Contribution of cues for perceptual organization
%       to the cross-spectral fusion of concurrent narrow-band noises," 
%       Journal of the Acoustical Society of America 111(4), pp. 1819?1831.

    properties (Dependent = true)
        minStrengthdB   % Minimum transient strength for mapping
        minSpread       % Minimum spread of the transient (number of frequency channels)
        fuseWithinSec   % Events within that period (in sec) are fused together
        minValuedB      % Lower limit for the input representation below which transient are discarded
    end
    
    properties (GetAccess = private)
        buffer          % Last fuseWithinSec seconds of input chunk are stored there
        fuseWithinSamples
    end
    
    methods
        function pObj = onsetMapProc(fs,parObj)
        %onsetMapProc   Construct an onset mapping processor
        %
        % USAGE:
        %   pObj = onsetMapProc(fs, parObj)
        %
        % INPUT ARGUMENTS:
        %     fs : Input sampling frequency (Hz)
        % parObj : Parameter object instance
        %
        % OUTPUT ARGUMENTS:
        %   pObj : Processor instance
        %
        % NOTE: Parameter object instance, parObj, can be generated using genParStruct.m
        % User-controllable parameters for this processor and their default values can be
        % found by browsing the script parameterHelper.m
        %
        % See also: genParStruct, parameterHelper, Processor
    
            % Checking input parameter
            if nargin<2||isempty(parObj); parObj = Parameters; end
            if nargin<1; fs = []; end
            
            % Call superconstructor
            pObj = pObj@Processor(fs,fs,'onsetMapProc',parObj);
            
            if nargin > 0
                pObj.buffer = [];
            end
            
        end
    
        function out = processChunk(pObj,in)
            %processChunk       Apply the processor to a new chunk of input signal
            %
            %USAGE
            %   out = pObj.processChunk(in)
            %
            %INPUT ARGUMENT
            %    in : New chunk of input data
            %
            %OUTPUT ARGUMENT
            %   out : Corresponding output
    
            % Append the buffer
            in = [pObj.buffer; in];
            
            % Store the last fuseWithinSamples of the input in the buffer
            Lbuf = size(pObj.buffer,1);
            L = size(in,1);
            
            if L <= pObj.fuseWithinSamples
                % Then the buffered input goes back into the buffer and there is no
                % processing for this chunk
                pObj.buffer = in;
                out = [];
            else
                
                % Discard transients if the representation is below a threshold
                if ~isempty(pObj.minValuedB)
                    try
                        rmap = pObj.LowerDependencies{1}.LowerDependencies{1}.Output{1}.Data(end-L+1:end);   
                        bSet2zero = 10*log10(rmap) < pObj.minValuedB;
                        in(bSet2zero) = 0;
                    
                    catch
                        warning('Could not access the ratemap representation from which transients were derived. Skipping transient discarding ...')
                    end
                        
                end
                
                % This "valid" input is then processed
                out = detectOnsetsOffsets(in,1/pObj.FsHzIn,pObj.minStrengthdB,...
                            pObj.minSpread,pObj.fuseWithinSec);
                
                % The output for the first Lbuf samples was already provided in the
                % previous chunk
                out = out(Lbuf+1:end,:);
                
                % Update the buffer
                pObj.buffer = in(end-pObj.fuseWithinSamples+1:end,:);
                
            end
            
            
        end
    
        function reset(pObj)
            pObj.buffer = [];
        end
        
        function output = instantiateOutput(pObj,dObj)
            %INSTANTIATEOUTPUT  Instantiate the output signal for this processor
            %
            %NB: This method is overloaded here from the master Processor class, as
            %binary mask signals need additional input argument to construct
            
            maskedSignal = pObj.LowerDependencies{1}.LowerDependencies{1}.Output{1};
            
            sig = feval(pObj.getProcessorInfo.outputType, ...
                        pObj, ...
                        dObj.bufferSize_s, ...
                        pObj.Channel,...
                        [],...
                        maskedSignal);
            
            dObj.addSignal(sig);
            
            output = {sig};
            
        end
        
    end
    
    methods (Hidden = true)
        
        function prepareForProcessing(pObj)
            
            % Compute internal parameter
            pObj.fuseWithinSamples = ceil(pObj.fuseWithinSec*pObj.FsHzIn);
            
        end
        
    end
    
    % "Getter" methods
    methods
        
        function minStrengthdB = get.minStrengthdB(pObj)
            minStrengthdB = pObj.parameters.map('trm_minStrengthdB');            
        end
        
        function minSpread = get.minSpread(pObj)
            minSpread = pObj.parameters.map('trm_minSpread');
        end
        
        function fuseWithinSec = get.fuseWithinSec(pObj)
            fuseWithinSec = pObj.parameters.map('trm_fuseWithinSec');            
        end
        
        function minValuedB = get.minValuedB(pObj)
            minValuedB = pObj.parameters.map('trm_minValuedB');            
        end
        
    end
    
    methods (Static)
        
        function dep = getDependency()
            dep = 'onsetStrength';
        end
        
        function [names, defaultValues, descriptions] = getParameterInfo()
            %getParameterInfo   Returns the parameter names, default values
            %                   and descriptions for that processor
            %
            %USAGE:
            %  [names, defaultValues, description] =  ...
            %                           gammatoneProc.getParameterInfo;
            %
            %OUTPUT ARGUMENTS:
            %         names : Parameter names
            % defaultValues : Parameter default values
            %  descriptions : Parameter descriptions
            
            
            names = {'trm_minStrengthdB',...
                     'trm_minSpread',...
                     'trm_fuseWithinSec',...
                     'trm_minValuedB'};
            
            descriptions = {'Minimum transient strength for mapping',...
                            'Minimum spread of the transient over frequency channels',...
                            'Events within that period (in sec) are fused together',...
                            ['Lower limit of the original representation, below which' ...
                            ' its transient will not be considered']};
            
            defaultValues = {3,...
                             5,...
                             30E-3,...
                             []};
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'Onset mapping';
            pInfo.label = 'Onset mapping';
            pInfo.requestName = 'onsetMap';
            pInfo.requestLabel = 'Onset map';
            pInfo.outputType = 'BinaryMask';
            pInfo.isBinaural = false;
            
        end
        
    end
    
end