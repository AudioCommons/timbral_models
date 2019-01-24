classdef genericFilter < filterObj
    
    methods 
        function obj = genericFilter(b,a,fs,cascade)
            % genericFilter Design a "generic" filter object from
            %               its transfer function coefficients
            %
            %USAGE
            %       F = genericFilter(b,a)
            %       F = genericFilter(b,a,fs,struct,cascade)
            %
            %INPUT ARGUMENTS
            %       b : filter coefficient numerator
            %       a : filter coefficient denominator
            %      fs : sampling frequency on which the filter operates
            %  struct : implementation of filter structure (default:
            %           'Direct-Form II Transposed')
            % cascade : Cascading order of the filter (default : 1)
            %
            %OUTPUT ARGUMENT
            %       F : filter object
            
            if nargin>0
                % CHECK INPUT ARGUMENTS
                if nargin<2
                    error('Provide both numerator and denominator filter coefficients')
                end
                if nargin<3 || isempty(fs); fs=1; end
                if nargin<4 || isempty(cascade); cascade = 1; end

                % POPULATE THE FILTER OBJECT PROPERTIES
                obj = populateProperties(obj,'Type','Generic Filter','FsHz',fs,...
                    'b',b,'a',a);
                
                obj.CascadeOrder = cascade;
                obj.RealTF = true;
            end
        end
    end
end