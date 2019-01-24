classdef templateProc < Processor
%TEMPLATEPROC Blank template for creating new processors.
%
% Copy this file and rename it to the name of your processor, also changing the class name
% above and in the constructor. A quick search and replace for "templateProc" should allow
% not to miss any!
%
% Populate the fields below, with the help of the tutorial available at:
%   http://twoears.aipa.tu-berlin.de/doc/afe/add-processor.html
%
% Remember to remove such "tutorial" comments and document/comment your implementation
% Some h1 lines have been pre-populated, remember to update them with more accurate
% information on your processor. You can use the class definition of any of the existing
% processors as an inspiration or an coding style example.
%
% Some methods are or not to be implemented depending on your processor's behavior. They
% have been commented out here, and indications as well as example processors to refer to
% have been provided.
%
% NOTE: Some methods have never been overridden as it was not necessary for the current
% version of the AFE. The option of overriding them is available as it paves the road for
% using more complex processors. However, such features have not been tested so there
% could be a risk that further work is needed at a lower programming level.

%% Properties

    properties (Dependent = true)
        % List here the external properties of your processor, i.e., those associated with
        % a user-controllable parameter
        % The values of these properties are not stored here, but instead read via a
        % "getter" method implemented below
        par1            % Add a brief description here
        par2            
        % ...
        
    end
    
    properties (GetAccess = private)
        % List here internal properties of your processor that do not need public access
        % (e.g., internal buffers, filters,...)
        
    end
        
    % You can add more properties with different attributes if needed
    
%% Methods
    methods
        function pObj = templateProc(fs,parObj)
        %templateProc   Construct a processor
        %
        % USAGE:
        %   pObj = templateProc(fs, parObj)
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
            pObj = pObj@Processor(fs,fs,'templateProc',parObj);
            
            if nargin>0 && ~isempty(fs)
                % Add here any additional initialization steps which need only be called
                % once at instantiation. Additional initialization steps that should be
                % performed again after receiving feedback are implemented below
            end
        end
        
        function out = processChunk(pObj,in)
            %processChunk       Processing method of your processor
            %
            %USAGE
            %       out = processChunk(pObj,in)
            %       out = pObj.processChunk(in)
            %
            %INPUT ARGUMENTS
            %      pObj : templateProc processor object
            %        in : X-dimensional array containing the input signal
            %
            %OUTPUT ARGUMENTS
            %       out : X-dimensional array containing the output
            %
            %SEE ALSO:
            %       templateProc.m
            
            % Perform your processing here...
            
        end
        
        function reset(pObj)
            %reset          Order the processor to reset its internal
            %               states, e.g., when some critical parameters in
            %               the processing have been changed
            %USAGE
            %       pObj.reset()
            %       reset(pObj)
            %
            %INPUT ARGUMENT
            %       pObj : Processor object
            
            % Reset the internal states of your processor here, if any. 
            % Otherwise, leave empty.
            
        end
        
        % Uncomment and implement the following method only if your processor is one of
        % many alternative for a given processing step
        % See, e.g., gammatoneProc.m and drnlProc.m
        
%         function bInBranch = isSuitableForRequest(pObj)
%             % Look into pObj.parameters to determine if this processor is the right one
%             % considering the parameters entered by the user, return 1 if it is, or 0 if
%             % it is not.
%         end
        
    end
    
    methods (Access=protected)
        
        % Uncomment and implement the following method only if you should perform a check
        % on user-provided parameters or solve potential conflicts between parameters.
        % See, e.g., gammatoneProc.m or ihcProc.m
        
%         function verifyParameters(pObj)
%                         
%             % Solve any potential conflicts and/or control the validity of parameters
%             % here.
%             % Invalid parameter values are usually replaced with their corresponding
%             % default value.
%             
%         end
        
    end
    
    % "Overridden" methods
    methods (Hidden = true)
        % NOTE: All the methods in this block are for processors with "special" behavior.
        % Examples were provided when available.
        % The methods below are approximately ordered from the most likely to the least 
        % likely to be overridden.
        
        
        %% Uncomment and implement only if some internal parameters or properties of your
        % processor are subject to change following user feedback
        % All initialization steps that are potentially part of a "re-initialization"
        % following a change of parameter should go in here. If an internal parameter as a
        % value dependent on other processors this one depends on, it should be
        % initialized here, using e.g., the getDependentProperty or getDependentParameter
        % method of your processor.
        % See, e.g., ratemapProc.m
        
%         function prepareForProcessing(pObj)
%             
%             % Finalize the initialization of your processor. 
%             
%         end
        

        %% Uncomment and modify only if your processor generates multiple outputs 
        % (including left and right channels of a same type, e.g., pre-processor) or if
        % the signal type it generates needs additional arguments in its constructor
        % (e.g., featureSignal)
        % See, e.g., preProc.m and spectralFeaturesProc.m

%         function output = instantiateOutput(pObj,dObj)
%             %INSTANTIATEOUTPUT  Instantiate the output signal for this processor
%             %
%             % USAGE:
%             %  sObj = pObj.instantiateOutput(dObj)
%             %
%             % INPUT ARGUMENTS:
%             %  pObj : Processor instance
%             %  dObj : Data object instance
%             %
%             % OUTPUT ARGUMENTS:
%             %  sObj : Handle to the instantiated output signal 
%             %
%             % NB: This method should be overridden in children processor where output 
%             % differs from standard (e.g., multiple output or additional inputs to the 
%             % signal constructor).
%             
%             sig = feval(pObj.getProcessorInfo.outputType, ...
%                         pObj, ...
%                         dObj.bufferSize_s, ...
%                         pObj.Channel);
%             
%             dObj.addSignal(sig);
%             
%             output = {sig};
%             
%         end

        %% Uncomment and modify only if your processor is not a single or binaural input 
        % and single output processor.
        % See, e.g., preProc.m

%         function initiateProcessing(pObj)
%             %INITIATEPROCESSING    Wrapper calling the processChunk method and routing I/O
%             % Main purpose is to allow overloading of I/O routing in processors with
%             % "unusual" number of input/outputs.
%             %
%             % Two cases considered here, monaural and binaural processors producing single
%             % outputs. In other cases, the method should be overloaded in the particular
%             % processor.
%             %
%             % USAGE:
%             %   pObj.initiateProcessing
%             %
%             % INPUT ARGUMENT:
%             %   pObj : Processor instance
%             
%             if size(pObj.Input,1) > 1 || numel(pObj.Output) > 1
%                 % Then it is a multiple-input processor, return an error
%                 error(['Cannot initiate the processing for this ' ...
%                     'processor. Consider overloading this method in the children ' ... 
%                     'class definition.'])
%             
%             elseif size(pObj.Input,2) == 1
%                 % Then monaural processor
%                 pObj.Output{1}.appendChunk( ...
%                     pObj.processChunk( pObj.Input{1}.Data('new')));
%                 
%             elseif size(pObj.Input,2) == 2
%                 % Then binaural processor
%                 pObj.Output{1}.appendChunk( ...
%                     pObj.processChunk( ...
%                         pObj.Input{1}.Data('new'), pObj.Input{2}.Data('new')));
%                 
%                 
%             else
%                 % TODO: Remove after testing
%                 error('Something is wrong with the inputs of this processor, investigate.')
%             end
%             
%         end


        %% Uncomment and modify only if your processor depends on multiple type of other
        % processors
        % NOTE: There are no such processor in the AFE at the moment, so you're on your
        % own!
        
%         function addUpperDependencies(pObj,dependentProcs)
%             %ADDUPPERDEPENDENCIES   Populate link to higher processors relying on this one
%             %
%             % USAGE:
%             %  pObj.addUpperDependencies(dependentProcs)
%             %
%             % INPUT ARGUMENTS:
%             %           pObj : Processor instance
%             % dependentProcs : Cell array of (correctly ordered) processors dependent on
%             %                  this processor
%             
%             pObj.UpperDependencies = [pObj.UpperDependencies dependentProcs];
%             
%         end


        %% Uncomment and modify only if your processor depends on multiple type of other
        % processors
        % NOTE: There are no such processor in the AFE at the moment, so you're on your
        % own!

%         function removeUpperDependency(pObj,upperProc)
%             %REMOVEUPPERDEPENDENCY   Removes a processor from the upper dependency list
%             %
%             % USAGE:
%             %    pObj.removeUpperDependency(upperProc)
%             %
%             % INPUT ARGUMENTS:
%             %       pObj : Processor instance
%             %  upperProc : Processor depending on pObj to remove from the dependency list
%             
%             for ii = 1:size(pObj.UpperDependencies,2)
%                 if pObj.UpperDependencies{ii} == upperProc
%                     pObj.UpperDependencies{ii} = [];
%                 end
%             end
%             
%             % Clean up empty elements
%             pObj.UpperDependencies = ...
%                 pObj.UpperDependencies( ~cellfun( @isempty, pObj.UpperDependencies ));
%             
%         end
%         


        %% Uncomment and implement only if your processor has multiple inputs of different
        % types (i.e., not just 'left' and 'right' channel of a same type)
        % NOTE: There are no such processor in the AFE at the moment, so you're on your
        % own!

%         function addInput(pObj,dependency)
%             %ADDINPUT   Adds a signal to the list of input of the processor
%             % Will consider that the input signal is the output of the dependency. If
%             % the dependency has a left- and right-channel output, will pick the suitable
%             % one.
%             % Three cases are implemented here:
%             % 1- single dependency with single output
%             % 2- two dependencies (left and right channels) with each single output
%             % 3- single dependency with left and right outputs (e.g., pre-processor)
%             %
%             % Should the input attribution works in any other way for a given processor,
%             % this method should be overloaded for that specific children processor.
%             %
%             
%             
%         end


        %% Uncomment and implement only if your processors returns multiple outputs of
        % different types (i.e., not just 'left' and 'right' channels of a same type)

%         function addOutput(pObj,sObj)
%             %ADDOUTPUT  Adds a signal to the list of output of the processor
%             %
%             % USAGE:
%             %   pObj.addOutput(sObj)
%             %
%             % INPUT ARGUMENTS:
%             %   pObj : Processor instance
%             %   sObj : Single instance or cell array of left and right channel instances
%             %          of signal objects
%             
%             
%         end

    end
    
    %% "Getter" methods
    methods
        % All properties listed as "Dependent" in the beginning of the class definition
        % should have a corresponding "getter" method here.
        
        function par1Val = get.par1(pObj)
            % For example, par1 corresponds to the user-controlled parameter 'XX_para1'
            par1Val = pObj.parameters.map('XX_para1');
        end
        
        function par2Val = get.par2(pObj)
            % Another example, par2 is just used to simplify the code and avoid having a
            % long expression re-occuring throughout the code.
            par2Val = pObj.par1 + sqrt(pObj.par1)*exp(pObj.par1); % Complete dummy value!
        end
        
        % ...
        
    end
    
    %% Static methods
    % These methods store all hard-coded information regarding your processor, remember to
    % update them with actual information!
    methods (Static)
        
        function dep = getDependency()
            dep = 'theSignalNameITakeAsInput';
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
            
            
            names = {'XX_para1',...
                    };
            
            descriptions = {'A first parameter that does something (unit)',...
                    };
            
            defaultValues = {'Some default value of any type',...
                            };
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'Short name';
            pInfo.label = 'Longer, more explicit name';
            pInfo.requestName = 'myNewRequest';
            pInfo.requestLabel = 'A description of my new request';
            pInfo.outputType = 'ClassNameOfOutputSignal';
            pInfo.isBinaural = false; % 0 or 1 (or 2 if it can do both, e.g., preProc.m)
            
        end
        
    end
        
end