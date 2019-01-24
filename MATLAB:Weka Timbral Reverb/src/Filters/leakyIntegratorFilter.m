classdef leakyIntegratorFilter < filterObj
    
    properties
        decaySec    % Integration time constant in seconds
    end
    
    methods
        function obj = leakyIntegratorFilter(fs,decaySec,cascade)
            %leakyIntegrator   Design leaky integration object filter.
            %
            %USAGE
            %          F = leakyIntegrator(fs)
            %          F = leakyIntegrator(fs,decaySec,cascade)
            %
            %INPUT ARGUMENTS
            %         fs : sampling frequency in Hz
            %   decaySec : leaky integration time constant in seconds
            %              (default, decaySec = 8E-3)
            %    cascade : filter cascading order (default : 1)
            %
            %OUTPUT ARGUMENTS
            %          F : filter object
            %EXAMPLE
            %   genFilterLeakyIntegrator(44.1E3);

            %   Developed with Matlab 7.5.0.342 (R2007b). Please send bug reports to:
            %   
            %   Author  :  Tobias May, © 2007-2008 
            %              TUe Eindhoven and Philips Research  
            %              t.may@tue.nl      tobias.may@philips.com
            %
            %   History :   
            %   v.0.1   2008/07/29
            %   v.0.2   2009/05/17 added 3Dim input support
            %   ***********************************************************************
            
            % CHECK INPUT ARGUMENTS 
            % 
            %
            if nargin>0 % Need to allow the constructor to be called without arguments
                
                if nargin > 3
                    help(mfilename);
                    error('Wrong number of input arguments!')
                end

                % Set default value
                if nargin < 3 || isempty(cascade); cascade = 1; end
                if nargin < 2 || isempty(decaySec); decaySec = 8E-3; end


                % COMPUTE FILTER COEFFICIENTS
                % 
                % 
                % Filter decay
                intDecay = exp(-(1/(fs * decaySec)));

                % Integration gain
                intGain = 1 - intDecay;

                % Set up standard filter object properties
                obj = populateProperties(obj,'Type','Leaky Integrator','FsHz',fs,...
                    'b',intGain,'a',[1 -intDecay]);
                obj.CascadeOrder = cascade;

                % Set up specific properties
                obj.decaySec = decaySec;
            end            
        end
    end
    
    
end