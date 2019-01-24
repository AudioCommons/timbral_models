classdef bwFilter < filterObj
    
    methods
        function obj = bwFilter(fs,order,cutOffHz,cascade,type)
            if nargin > 0
                if nargin < 5||isempty(type);type = 'low';end
                if nargin < 4||isempty(cascade);cascade = 1;end
                if nargin < 3||isempty(cutOffHz); cutOffHz = 1000; end % Dau1996 model
                if nargin < 2||isempty(order); order = 2; end % Dau1996
                
                % Generate filter coefficients
                [b,a]=butter(order,cutOffHz/(0.5*fs),type);
                
                % Populate filter properties
                obj = populateProperties(obj,'Type', 'Butterworth low-pass filter','FsHz',fs,'b',b,'a',a);
            
                obj.CascadeOrder = cascade;
               
                
            end
        end
        
        
    end
    
end