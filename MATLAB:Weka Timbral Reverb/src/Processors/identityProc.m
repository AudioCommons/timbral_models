classdef identityProc < Processor
%IDENTITYPROC A processor that copies the input directly to the output.
%   Used mainly for place holder, e.g., for instantiating arrays of processors in the 
%   manager. 
    
    methods
        function pObj = identityProc(fs,parObj)
            %identityProc   This processor does nothing! Used as a placeholder in the AFE.
            
            if nargin<1; fs = []; end
            if nargin<2||isempty(parObj); parObj = Parameters; end
            
            pObj = pObj@Processor(fs,fs,'identityProc',parObj);
            
            % Hide the processor from the list of processors
            pObj.bHidden = 1;
        end
        
        function out = processChunk(~,in)
            out = in;
        end
        
        function reset(~)
            % Nothing to reset here
        end
                
    end
    
    methods (Static)
        
        function dep = getDependency()
            dep = 'n-a';
        end
        
        function [names, defaultValues, description] = getParameterInfo()
            % No parameters here
            names = [];
            defaultValues = [];
            description = [];
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'Identity processor';
            pInfo.label = 'Identity processor';
            pInfo.requestName = 'n-a';
            pInfo.requestLabel = 'n-a';
            pInfo.outputType = 'n-a';
            
        end
        
    end
    
end