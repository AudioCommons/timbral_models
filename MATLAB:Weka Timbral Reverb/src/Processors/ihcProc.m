classdef ihcProc < Processor
    
     properties (Dependent = true)
         method        % Label for the IHC model used
     end
     
     properties %(GetAccess = private)
         IHCFilter     % Filter involved in the extraction, if any
     end
     
    methods
        function pObj = ihcProc(fs,parObj)
		%ihcProc   Construct an inner hair-cell envelope extractor processor
        %
        % USAGE:
        %   pObj = ihcProc(fs, parObj)
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
             
             
             if nargin<2||isempty(parObj); parObj = Parameters; end
             if nargin<1; fs = []; end
             
             % Call super-constructor
             pObj = pObj@Processor(fs,fs,'ihcProc',parObj);
             
        end
         
         function out = processChunk(pObj,in)
                        
            % Carry out the processing for the chosen IHC method
            switch pObj.method
                case 'none'
                    out = in;

                case 'halfwave'
                    % Half-wave rectification
                    out = max(in,0);

                case 'fullwave'
                    % Full-wave rectification
                    out = abs(in);

                case 'square'
                    out = abs(in).^2;

                case 'hilbert'
                    out = abs(hilbert(in));

                case 'joergensen'
                    out = pObj.IHCFilter.filter(abs(hilbert(in)));

                case 'dau'
                    out = pObj.IHCFilter.filter(max(in,0));

                case 'breebart'
                    out = pObj.IHCFilter.filter(max(in,0));

                case 'bernstein'
                    env = max(abs(hilbert(in)).^(-.77).*in,0).^2;
                    out = pObj.IHCFilter.filter(env);

                otherwise
                    error('%s: Method ''%s'' is not supported!',upper(mfilename),pObj.IHCMethod)
            end
            
         end
         
         function reset(pObj)
             %reset     Resets the internal states of the IHC envelope
             %          extractor
             %
             %USAGE
             %      pObj.reset
             %
             %INPUT ARGUMENTS
             %  pObj : Inner haircell envelope extractor processor instance
             
             % A reset is needed only if the extractor involves filters
             if ~isempty(pObj.IHCFilter)
                 pObj.IHCFilter.reset
             end
         end
         
     end
      
     methods (Access=protected)
         
         function verifyParameters(pObj)
             
             % Check that the IHC method name is valid
             validMeth = {'none',...
                         'halfwave',...
                         'fullwave',...
                         'square',...
                         'hilbert',...
                         'joergensen',...
                         'dau',...
                         'breebart',...
                         'bernstein'};
             
             if ~ismember(pObj.parameters.map('ihc_method'),validMeth)
                 [~,defaultMethod] = ihcProc.getParameterInfo;
                 warning(['''%s'' is an invalid name for envelope extraction method. '...
                          'Setting it to the default value, ''%s'''],...
                          pObj.parameters.map('ihc_method'),defaultMethod)
                 pObj.parameters.map('ihc_method') = defaultMethod;
             end
             
         end
         
     end
     
     methods (Hidden = true)
         
         function update(pObj,~,~)
            % Overloading the update method for that processor, for testing purposes
            %disp(['I am ' pObj.Type ', and I will not display the same message as the ' ...
            %   'others because my update method is overloaded!'])
            
            % Update means reset
            pObj.reset;
            
            % Notify possible listeners that there was a modification
            notify(pObj,'hasChanged');
            
         end
        
         function prepareForProcessing(pObj)
             
             % Instantiate the filters
             pObj.populateFilters;
             
         end
         
     end
     
     methods (Static)
        function dep = getDependency()
            dep = 'filterbank';
        end
        
        function [names, defaultValues, descriptions] = getParameterInfo()
            %getParameterInfo   Returns the parameter names, default values
            %                   and descriptions for that processor
            %
            %USAGE:
            %  [names, defaultValues, description] =  ihcProc.getParameterInfo;
            %
            %OUTPUT ARGUMENTS:
            %         names : Parameter names
            % defaultValues : Parameter default values
            %  descriptions : Parameter descriptions
            
            
            names = {'ihc_method'};
            
            descriptions = {['Inner hair-cell envelope extraction method (''none'', ' ...
                            '''halfwave'', ''fullwave'', ''square'', ''hilbert'', '...
                            '''joergensen'', ''dau'', ''breebart'', ''berstein'')']};
            
            defaultValues = {'dau'};
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'IHC envelope';
            pInfo.label = 'Inner hair-cell envelope extraction';
            pInfo.requestName = 'innerhaircell';
            pInfo.requestLabel = 'Inner hair-cell envelope';
            pInfo.outputType = 'TimeFrequencySignal';
            pInfo.isBinaural = 0;
            
        end
        
     end
    
     % "Getter" methods
     methods
         function method = get.method(pObj)
             method = pObj.parameters.map('ihc_method');
         end
     end
     
     methods (Access = private)
         function populateFilters(pObj)
             
             % Using a try/catch to allow instantiation of dummy ihcProc
             try
                 % Instantiate a low-pass filter if needed
                 switch pObj.method

                     case 'joergensen'
                         % First order butterworth filter @ 150Hz
                         pObj.IHCFilter = bwFilter(pObj.FsHzIn,1,150);

                     case 'dau'
                         % Second order butterworth filter @ 1000Hz
                         pObj.IHCFilter = bwFilter(pObj.FsHzIn,2,1000);

                     case 'breebart'
                         % First order butterworth filter @ 2000Hz cascaded 5
                         % times
                         pObj.IHCFilter = bwFilter(pObj.FsHzIn,1,2000,5,[]);

                     case 'bernstein'
                         % Second order butterworth filter @ 425Hz
                         pObj.IHCFilter = bwFilter(pObj.FsHzIn,2,425);

                     otherwise
                         pObj.IHCFilter = [];
                 end
             catch
                 pObj.IHCFilter = [];
             end
         end
     end
     
end