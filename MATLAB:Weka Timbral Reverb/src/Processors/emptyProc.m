classdef emptyProc < Processor
%emptyPROC A processor that is used only for place holder and will return an error when
%trying to processo anything
    
    methods
        function pObj = emptyProc(fs,parObj)
            %emptyProc   This processor does nothing! Used as a placeholder in the AFE.
            
            if nargin<1; fs = []; end
            if nargin<2||isempty(parObj); parObj = Parameters; end
            
            pObj = pObj@Processor(fs,fs,'emptyProc',parObj);

            % Hide the processor from the list of processors
            pObj.bHidden = 1;
        end
        
        function out = processChunk(~,~) %#ok<STOUT>
            % Throw an error here
            error('This is an empty processor and it cannot perform any processing')
        end
        
        function reset(pObj) %#ok<MANU>
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
            
            pInfo.name = 'Empty processor';
            pInfo.label = 'Empty processor';
            pInfo.requestName = 'n-a';
            pInfo.requestLabel = 'n-a';
            pInfo.outputType = 'n-a';
            
        end
        
    end
    
end