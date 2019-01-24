classdef manager < handle
%MANAGER Processor managing class for the auditory front-end (AFE) framework. A manager 
%   object controls the processing of the AFE framework. It is responsible for 
%   instantiating the required processors as well as correctly routing their respective 
%   inputs/outputs, given a request from the user. In addition, the manager methods allow 
%   the user to request a new representation or ask for the processing to be performed. 
%   Hence, the manager object represents the core of the AFE framework. 
%
%   MANAGER properties:
%       Processors - Cell array of processor objects.
%       InputList  - Handles to the input of each processors.
%       OutputList - Handles to the output of each processors.
%       Data       - Handle to the data object containing all computed signals.
%
%   MANAGER methods:
%       manager       - Constructor for the class. Requires a dataObject instance.
%       addProcessor  - Request a new auditory representation to extract.
%       processSignal - Requests the (offline) processing of an input signal.
%       processChunk  - Requests the (online) processing for a new chunk of signal..
%       hasProcessor  - Test if a given processor is already instantiated.
%       reset         - Resets internal states of all processors.
%       
%   See also dataObject, requestList, parameterHelper
%
% Disclamer: Known limitations that will be addressed in future releases
%   - When a processor becomes obsolete, its instance is not cleared from memory
%   - Few processors are not fully compatible with chunk-based processing (will return
%     erroneous representations in the vicinity of chunk boundaries):
%       * IHC methods involving Hilbert envelope extraction ('hilbert', 'joergensen', 
%         and 'bernstein')
%       * Spectro-temporal modulation extraction
%   - Pitch estimation might (though highely unlikely) be misestimated at chunk boundaries
%     in chunk-based processing scenarios
    
    
    properties (SetAccess = protected)
        % Processors - Cell array of processor objects. First column of the array contains
        % processors in charge of the left (or single) channel, second column of the
        % right channel. Different lines in the array are for different processor
        % instances.
        Processors = {};     
        
        % InputList - Cell array of handles to the input signal of each processors. A
        % signal at a given position in the array is the input to the processor stored
        % at the same position in the Processors property.
        InputList       
        
        % OutputList - Cell array of handles to the output signal of each processors. A
        % signal at a given position in the array is the output from the processor stored
        % at the same position in the Processors property.
        OutputList      
        
        % Data - Handle to the data object associated with this instance of the manager.
        Data
        
    end
    
    properties (GetAccess = protected)
        use_mex         % Flag for using mex files to speed up computation 
                        % when available
        Map             % Vector mapping the processing order to the 
                        % processors order. Allows for avoiding to reorder
                        % the processors array when new processors are
                        % added.
    end
    
    
    methods
        function mObj = manager(data,request,p,use_mex)
            %manager    Constructs a manager object
            %
            %USAGE
            %     mObj = manager(data)
            %     mObj = manager(data,request)
            %     mObj = manager(data,request,p)
            %
            %INPUT ARGUMENTS
            %     data : Handle of an existing data structure
            %  request : Single request as a string (e.g., 'ild'), OR cell array of
            %            requested signals, cues or features.
            %        p : Single parameter structure, if all requests share the same 
            %            parameters, OR cell array of individual parameter structures 
            %            corresponding to each request.
            %
            %OUTPUT ARGUMENTS
            %     mObj : Manager instance
            %
            %EXAMPLE USE (given an instance of dataObject, dObj)
            %- 'Empty' manager:
            %   mObj = manager(dObj)
            %- Single request, default parameters:
            %   mObj = manager(dObj,'autocorrelation')
            %- Multiple request with same parameters
            %   mObj = manager(dObj,{'ild','itd'},genParStruct('fb_nChannels',16))
            %  
            %
            %SEE ALSO: dataObject requestList genParStruct
            
            if nargin>0     % Failproof for Matlab empty calls
            
            % Input check
            if nargin<4||isempty(use_mex);use_mex=1;end
            if nargin<3||isempty(p);p=[];end
            if nargin<2
                request = [];
            end
            if nargin<1
                error(['Too few arguments, the manager is built upon '...
                    'an existing data Object'])
            end
            
            % Add use_mex property for the manager
            mObj.use_mex = use_mex;
            
            % Add pointer to the data structure
            mObj.Data = data;
            
            % Instantiate the requested processors
            if ~isempty(request)
                if iscell(request) && numel(request) == 1
                    % Then we have a one request with multiple parameters
                    if iscell(p)
                        %... with individual parameters
                        for ii = 1:size(p,2)
                            mObj.addProcessor(request,p{ii});
                        end
                    else
                        mObj.addProcessor(request,p);
                    end
                elseif iscell(request)
                    % Then we have a multiple request...
                    if iscell(p)
                        %... with individual parameters
                        if size(request,2)~=size(p,2)
                            error('Number of requests and number of provided parameters do not match')
                        else
                            for ii = 1:size(request,2)
                                mObj.addProcessor(request{ii},p{ii});
                            end
                        end
                    else
                        %... all with the same set of parameters
                        for ii = 1:size(request,2)
                            mObj.addProcessor(request{ii},p);
                        end
                    end
                elseif iscell(p)
                    % Then it is a same request but with multiple parameters
                    for ii = 1:size(p,2)
                        mObj.addProcessor(request,p{ii});
                    end
                else
                    % Then it is a single request
                     mObj.addProcessor(request,p);
                end
            end
            end
        end
        
        function processSignal(mObj,s)
            %processSignal      Requests a manager object to extract the requested 
            %                   features for the complete signal in mObj.Data.input
            %
            %USAGE
            %    mObj.processSignal()
            %
            %INPUT ARGUMENT
            %   mObj : Manager object
            %
            %NB: As opposed to the method processChunk, this method will
            %reset the internal states of the processors prior to
            %processing, assuming a completely new signal.
            %
            %SEE ALSO: processChunk
            
            if nargin == 2 && ~isempty(s)
                % Then use the input signal provided as input
                mObj.Data.replaceInputSignal(s);
            end
            
            % Check that there is an available signal
            if isempty(mObj.Data.input)
                warning('No signal available for processing')
            else            
                % Reset the processors internal states
                mObj.reset;
                
                % Number of processors
                n_proc = size(mObj.Processors,1);

                % Loop on each processor
                for ii = 1:n_proc
                    % Get index of current processor
                    jj = mObj.Map(ii);

                    mObj.Processors{jj,1}.initiateProcessing;
                    
                    if size(mObj.Processors,2) == 2 && ~isempty(mObj.Processors{jj,2})
                        mObj.Processors{jj,2}.initiateProcessing;
                    end
                end
            end
        end
        
        function processChunk(mObj,sig_chunk,do_append)
            %processChunk   Update the signal with a new chunk of data and calls the 
            %               processing chain for this new chunk.
            %
            %USAGE
            %   mObj.processChunk(sig_chunk)
            %   mObj.processChunk(sig_chunk,append)
            %
            %INPUT ARGUMENTS
            %      mObj : Manager object
            % sig_chunk : New signal chunk
            %    append : Flag indicating if the newly generated output
            %             should be appended (append = 1) to previous
            %             output or should overwrite it (append = 0,
            %             default)
            %
            %NB: Even if the previous output is overwritten, the
            %processChunk method allows for chunk-based processing by keeping
            %track of the processors' internal states between chunks.
            %
            %SEE ALSO: processSignal
            
            if nargin<3||isempty(do_append);do_append = 0;end
            
            % Check that the signal chunk has correct number of channels
            if size(sig_chunk,2) ~= mObj.Data.isStereo+1
                % TO DO: Change that to a warning and handle appropriately
                error(['The dimensionality of the provided signal chunk'...
                    'is incompatible with previous chunks'])
            end
            
            % Delete previous output if necessary
            if ~do_append
                mObj.Data.clearData;
            end
            
            
            % Append the signal chunk
            if mObj.Data.isStereo
               mObj.Data.input{1}.appendChunk(sig_chunk(:,1));
               mObj.Data.input{2}.appendChunk(sig_chunk(:,2));
            else            
               mObj.Data.input{1}.appendChunk(sig_chunk);
            end
            
            % Number of processors
            n_proc = size(mObj.Processors,1);
            
            % Loop on each processor
            for ii = 1:n_proc
                % Get index of current processor
                jj = mObj.Map(ii);

                mObj.Processors{jj,1}.initiateProcessing;

                if size(mObj.Processors,2) == 2 && ~isempty(mObj.Processors{jj,2})
                    mObj.Processors{jj,2}.initiateProcessing;
                end
            end
            
            % Loop on each processor
%             for ii = 1:n_proc
%                 % Get index of current processor
%                 jj = mObj.Map(ii);
%                 
%                 if ~mObj.Processors{jj,1}.isBinaural
%                     % Apply processing for left channel (or mono if
%                     % interaural cue/feature):
% 
%                     % Getting input signal handle (for code readability)
%                     in = mObj.InputList{jj,1};
% 
%                     % Perform the processing
%                     out = mObj.Processors{jj,1}.processChunk(in.Data('new'));
% 
%                     % Store the result
%                     mObj.OutputList{jj,1}.appendChunk(out);
% 
%                     % Apply similarly for right channel if binaural cue/feature
%                     if mObj.Data.isStereo && ~isempty(mObj.Processors{jj,2})
%                         in = mObj.InputList{jj,2};
%                         out = mObj.Processors{jj,2}.processChunk(in.Data('new'));
%                         mObj.OutputList{jj,2}.appendChunk(out);
%                     end
%                     
%                 else
%                     % Inputs from left AND right channels are needed at
%                     % once
%                     
%                     % Getting input signal handles for both channels
%                     in_l = mObj.InputList{jj,1};
%                     
%                     if ~mObj.Processors{jj,1}.hasTwoOutputs
%                         
%                         in_r = mObj.InputList{jj,2};
%                         
%                         % Perform the processing
%                         out = mObj.Processors{jj,1}.processChunk(...
%                             in_l.Data('new'),...
%                             in_r.Data('new'));
% 
%                         % Store the result
%                         mObj.OutputList{jj,1}.appendChunk(out);
%                     else
%                         
%                         if size(mObj.InputList,2)>1
%                             in_r = mObj.InputList{jj,2};
%                             
%                             % Perform the processing
%                             [out_l, out_r] = mObj.Processors{jj,1}.processChunk(...
%                                 in_l.Data('new'),...
%                                 in_r.Data('new'));
%                         else
%                             % Perform the processing
%                             [out_l, out_r] = mObj.Processors{jj,1}.processChunk(...
%                                 in_l.Data('new'));
%                         end
% 
%                         % Store the result
%                         mObj.OutputList{jj,1}.appendChunk(out_l);
%                         
%                         if ~isempty(out_r)
%                             mObj.OutputList{jj,2}.appendChunk(out_r);
%                         end
%                     end
%                 end
                
%                 % Getting input signal handle (for code readability)
%                 in = mObj.InputList{jj};
%                 
%                 % Perform the processing
%                 out = mObj.Processors{jj}.processChunk(in.Data('new'));
%                 
%                 % Store the result
%                 mObj.OutputList{jj}.appendChunk(out);
                
            
        end
        
        function hProc = hasProcessor(mObj,name,p,channel)
            %hasProcessor    Determines if a processor with a given set of parameters
            %                (including those of its dependencies) is already instantiated
            %
            %USAGE
            %   hProc = mObj.hasProcessor(name,p)
            %   hProc = mObj.hasProcessor(name,p,channel)
            %
            %INPUT ARGUMENTS
            %    mObj : Instance of manager object
            %    name : Name of processor
            %       p : Complete structure of parameters for that processor
            % channel : Channel the sought processor should be acting on
            %           ('left', 'right', or 'mono'). If unspecified, any
            %           processor with matching parameter will be returned.
            %
            %OUTPUT ARGUMENT
            %   hProc : Handle to an existing processor, if any, 0 else
            
            %TODO: Will need maintenance when introducing processors with multiple lower
            %dependencies
            
            ch_name = {'left','right','mono'};
            
            if nargin<4 %|| isempty(channel)
                channel = ch_name;
            elseif ~ismember(channel,ch_name)
                error('Invalid tag for channel name. Valid tags are as follow: %s',strjoin(ch_name))
            end
            
            if ~iscell(channel)
                channel = {channel};
            end
            
            % Initialize the output
            hProc = 0;
            
            % Look into corresponding ear depending on channel request.
            % Left and mono are always in the first column of the
            % processors cell array, right in the second.
            if strcmp(channel,'right')
                earIndex = 2;
            else
                earIndex = 1;
            end
            
            % Loop over the processors to find the ones with suitable name
            for ii = 1:size(mObj.Processors,1)
                
                % Get a handle to that processor, for readability in the
                % following
                proc = mObj.Processors{ii,earIndex};
                
                % Is the current processor one of the sought type?
                if isa(proc,name) && ismember(proc.Channel,channel)
                    
                    % Does it have the requested parameters?
                    if proc.hasParameters(p)
                        
                        % Then it is a suitable candidate, we should
                        % investigate its dependencies
                        while true
                            
                            if isempty(proc.LowerDependencies)
                                % Then we reached the end of the dependency
                                % list without finding a mismatch in
                                % parameters. The original processor is a
                                % solution:
                                hProc = mObj.Processors{ii,earIndex};
                                return
                            end
                            
                            % Set current processor to proc dependency
                            proc = proc.LowerDependencies{1};
                            
                            % Does the dependency also have requested
                            % parameters? If not, break of the while loop
                            if ~proc.hasParameters(p)
                                break
                            end
                            
                        end
                        
                        
                    end
                    
                end
                
                % If not, move along in the loop
                
            end
            
        end
        
        function [out,varargout] = addProcessor(mObj,request,p)
            %addProcessor   Add new processor(s) needed to compute a user request.
            %               Optionally returns a handle to the corresponding output signal
            %
            %USAGE:
            %           mObj.addProcessor(request,p)
            %    sOut = mObj.addProcessor(...)
            %
            %INPUT ARGUMENTS
            %    mObj : Manager instance
            % request : Requested signal (string)
            %       p : Structure of non-default parameters
            %
            %OUTPUT ARGUMENTS
            %    sOut : Handle to the requested signal
            %
            %EXAMPLE USE
            %- Single request, default parameters:
            %   sOut = mObj.addProcessor('autocorrelation');
            %- Multiple request with same non-default parameters
            %   [sOut1,sOut2] = manager({'ild','itd'}, genParStruct('fb_nChannels',16));
           
%             if nargin<3 || isempty(p)
%                 % Initialize parameter structure
%                 p = struct;
%             end
            
            if nargin<3; p = []; end
                

            % Deal with multiple requests via pseudo-recursion
            if iscell(request) || iscell(p)
                
                if iscell(request) && ~iscell(p)
                    % All the requests have the same parameters, replicate
                    % them
                    p = repmat({p},size(request));
                elseif ~iscell(request) && iscell(p)
                    % One request with different parameters, replicate the request
                    request = repmat({request},size(p));
                end
                
                if size(p,2)~=size(request,2)
                    error(['Provided number of parameter structures'...
                        ' does not match the number of requests made'])
                end
                
                % Call addProcessor method for each individual request
                varargout = cell(1,size(request,2)-1);
                out = mObj.addProcessor(request{1},p{1});
                for ii = 2:size(request,2)
                    varargout{ii-1} = mObj.addProcessor(request{ii},p{ii});
                end
                return
                
            end
            
%             if ~isfield(p,'fs')
%                 % Add sampling frequency to the parameter structure
%                 p.fs = mObj.Data.input{1}.FsHz;
%             end

            fs = mObj.Data.input{1}.FsHz;
            
            % Find most suitable initial processor for that request
            [initProc,dep_list] = mObj.findInitProc(request,p);
            
            % Replace the initProc with dummy processor(s) if empty
            if isempty(initProc)
                if mObj.Data.isStereo
                    initProc = {identityProc(fs), identityProc(fs)};
                    initProc{1}.Output = mObj.Data.input(1);
                    initProc{2}.Output = mObj.Data.input(2);
                else
                    initProc = {identityProc(fs)};
                    initProc{1}.Output = mObj.Data.input;
                end
            end
            
            % Algorithm should proceed further even if the requested
            % processor already exists
            if isempty(dep_list)
                proceed = 1;
            end
            
            % The processing order is the reversed list of dependencies
            dep_list = fliplr(dep_list);
 
            % Former and new number of processors
            n_proc = size(mObj.Processors,1);
            n_new_proc = size(dep_list,2);
            
            % Preallocation
            if isempty(mObj.Processors)
                if mObj.Data.isStereo
                    n_chan = 2;
                else
                    n_chan = 1;
                end
                mObj.Processors = cell(n_new_proc,n_chan);   
            end
            
            
            % Initialize pointer to dependency 
            dependency = initProc;
            
            % Processors instantiation and data object property population
            for ii = n_proc+1:n_proc+n_new_proc   
                
                proceed = 1;     % Initialize a flag to identify invalid requests (binaural representation requested on a mono signal)
                
                % Get the name of the processor to instantiate
                procName = Processor.findProcessorFromRequest(dep_list{ii-n_proc},p);
                
                
                % Check if one or two processors should be instantiated (mono or stereo)
                procInfo = feval([procName '.getProcessorInfo']);
                if size(dependency,2) == 2 && procInfo.isBinaural == 0
                    % Instantiate two processors, each having a single-channel dependency
                    newProc_l = mObj.addSingleProcessor(procName, p, dependency(1), 1, ...
                                                        ii, 'stereo');
                    newProc_r = mObj.addSingleProcessor(procName, p, dependency(2), 2, ...
                                                        ii, 'stereo');
                    dependency = {newProc_l, newProc_r};
                elseif numel(dependency) == 1 && size(dependency{1}.Output,2) == 2 ...
                                              && ~procInfo.isBinaural
                    % Instantiate two processors, each having the same multi-channel
                    % dependency
                    newProc_l = mObj.addSingleProcessor(procName, p, dependency, 1, ...
                                                        ii, 'stereo');
                    newProc_r = mObj.addSingleProcessor(procName, p, dependency, 2, ...
                                                        ii, 'stereo');
                    dependency = {newProc_l, newProc_r};
                elseif numel(dependency) == 1 && size(dependency{1}.Output,2) == 2 ...
                                              && procInfo.isBinaural
                    % Instantiate a single processor
                    newProc = mObj.addSingleProcessor(procName, p, dependency, 1, ...
                                                         ii, 'mono');
                    dependency = {newProc};
                else
                    if procInfo.isBinaural == 1 && ~mObj.Data.isStereo
                        warning(['Cannot instantiate a binaural processor with a '...
                            'mono input signal!'])
                        proceed = 0;
                    else
                        % Instantiate a single processor having a single dependency
                        newProc = mObj.addSingleProcessor(procName, p, dependency, 1, ii,...
                                                            'mono');
                        dependency = {newProc};
                    end
                end
                
                
                %% Old code again, commented for folding
                % Instantiate processors
%                 mObj.Processors{ii,1} = feval(procName, dep_proc{1}.FsHzOut, p);
%                 if mObj.Data.isStereo && ~mObj.Processors{ii,1}.isBinaural
%                     mObj.Processors{ii,2} = feval(procName, dep_proc{1}.FsHzOut, p);
%                 end
%                 
%                 mObj.findInitProc(mObj.Processors{ii,1}.getProcessorInfo.requestName,p) == dep_proc{1}
%                 
%                 % Link to dependencies
%                 if mObj.Processors{ii,1}.isBinaural
%                     mObj.Processors{ii,1}.Dependencies = dep_proc;
%                 else
%                     mObj.Processors{ii,1}.Dependencies = dep_proc(1);
%                     if mObj.Data.isStereo
%                         mObj.Processors{ii,2}.Dependencies = dep_proc(2);
%                     end
%                 end
%                 
%                 % Instantiate output signal
%                 sig = {feval(mObj.Processors{ii,1}.getProcessorInfo.outputType,...
%                             mObj.Processors{ii,1},...
%                             mObj.Data.bufferSize_s,...
%                             'mono')};
%                 if (mObj.Data.isStereo && ~mObj.Processors{ii,1}.isBinaural) || ...
%                         mObj.Processors{ii,1}.hasTwoOutputs
%                     sig = [sig feval(mObj.Processors{ii,1}.getProcessorInfo.outputType,...
%                                     mObj.Processors{ii,1},...
%                                     mObj.Data.bufferSize_s,...
%                                     'mono')];
%                 end
%                         
%                         
%                       
%                 % Add signal to the Data object
%                 mObj.Data.addSignal(sig);
                
%                 if ~isempty(mObj.Processors{ii})
%                 
%                     % Add input/output pointers, dependencies, and update dependencies.
%                     % Three possible scenarios:
% 
%                     if mObj.Processors{ii}.isBinaural
% 
%                         if ~mObj.Processors{ii}.hasTwoOutputs
%                             % 1-Then there are two inputs (left&right) and one output
%                             mObj.InputList{ii,1} = dep_sig_l;
%                             mObj.InputList{ii,2} = dep_sig_r;
%                             mObj.OutputList{ii,1} = sig;
%                             mObj.OutputList{ii,2} = [];
% 
%                             mObj.Processors{ii}.Input{1} = dep_sig_l;
%                             mObj.Processors{ii}.Input{2} = dep_sig_r;
%                             mObj.Processors{ii}.Output = sig;
% 
%                             mObj.Processors{ii,1}.Dependencies = {dep_proc_l,dep_proc_r};
%                             dep_sig = sig;
%                             dep_proc = mObj.Processors{ii};
%                         else
%                             if exist('sig','var')&&strcmp(sig.Channel,'mono')
%                                 % 1bis - Two inputs and two outputs
%                                 mObj.InputList{ii,1} = dep_sig;
%                                 mObj.OutputList{ii,1} = sig;
% 
%                                 mObj.Processors{ii}.Input = dep_sig;
%                                 mObj.Processors{ii}.Output = sig;
% 
%                                 mObj.Processors{ii,1}.Dependencies = {dep_proc};
%                                 dep_sig = sig;
%                                 dep_proc = mObj.Processors{ii};
%                             else
%                                 % 1bis - Two inputs and two outputs
%                                 mObj.InputList{ii,1} = dep_sig_l;
%                                 mObj.InputList{ii,2} = dep_sig_r;
%                                 mObj.OutputList{ii,1} = sig_l;
%                                 mObj.OutputList{ii,2} = sig_r;
% 
%                                 mObj.Processors{ii}.Input{1} = dep_sig_l;
%                                 mObj.Processors{ii}.Input{2} = dep_sig_r;
% %                                 mObj.Processors{ii}.Output{1} = sig_l;
% %                                 mObj.Processors{ii}.Output{2} = sig_r;
%                                 mObj.Processors{ii}.Output = sig_l;
%                                 
% 
% 
%                                 mObj.Processors{ii,1}.Dependencies = {dep_proc_l,dep_proc_r};
%                                 dep_sig_l = sig_l;
%                                 dep_sig_r = sig_r;
%                                 dep_proc_l = mObj.Processors{ii};
%                                 dep_proc_r = mObj.Processors{ii};
%                             end
%                         end
%                     elseif exist('sig','var')&&strcmp(sig.Channel,'mono') && proceed
% 
%                         % 2-Then there is a single input and single output
%                         mObj.InputList{ii,1} = dep_sig;
%                         mObj.OutputList{ii,1} = sig;
% 
%                         mObj.Processors{ii}.Input = dep_sig;
%                         mObj.Processors{ii}.Output = sig;
% 
%                         mObj.Processors{ii}.Dependencies = {dep_proc};
%                         dep_sig = sig;
%                         dep_proc = mObj.Processors{ii};
% 
%                     elseif ~proceed
% 
%                         % Do nothing, this request is invalid and should be
%                         % skipped
% 
%                     else
% 
%                         % 3-Else there are two inputs and two outputs
%                         mObj.InputList{ii,1} = dep_sig_l;
%                         mObj.InputList{ii,2} = dep_sig_r;
%                         mObj.OutputList{ii,1} = sig_l;
%                         mObj.OutputList{ii,2} = sig_r;
% 
%                         mObj.Processors{ii,1}.Input = dep_sig_l;
%                         mObj.Processors{ii,2}.Input = dep_sig_r;
%                         mObj.Processors{ii,1}.Output = sig_l;
%                         mObj.Processors{ii,2}.Output = sig_r;
% 
%                         mObj.Processors{ii,1}.Dependencies = {dep_proc_l};
%                         mObj.Processors{ii,2}.Dependencies = {dep_proc_r};
%                         dep_sig_l = sig_l;
%                         dep_sig_r = sig_r;
%                         dep_proc_l = mObj.Processors{ii,1};
%                         dep_proc_r = mObj.Processors{ii,2};
% 
%                     end
%                     
%                 else
%                     % Then the processor was not instantiated as the
%                     % request was invalid, exit the for loop
%                     break
%                 end
% 
%                 
%                 % Clear temporary handles to ensure no inconsistencies 
%                 clear sig sig_l sig_r
                
%% Resume
            end
            
            % The mapping at this point is linear
            mObj.Map(n_proc+1:n_proc+n_new_proc) = n_proc+1:n_proc+n_new_proc;
            
            % Provide the user with a pointer to the requested signal
            if nargout>0 && proceed
                if ~isempty(dep_list)
                    if size(mObj.Processors,2)==2
                        if isempty(mObj.Processors{n_proc+n_new_proc,2})
%                             out{1} = mObj.Processors{n_proc+n_new_proc,1}.Output{1};
                            % Modify to handle cases where one binaural
                            % processor has multiple outputs (precedence)
                            out = mObj.Processors{n_proc+n_new_proc,1}.Output;
                        else
                            out{1,1} = mObj.Processors{n_proc+n_new_proc,1}.Output{1};
                            out{1,2} = mObj.Processors{n_proc+n_new_proc,2}.Output{1};
                        end
                    else
                        out{1} = mObj.Processors{n_proc+n_new_proc,1}.Output{1};
                    end
                else
                    % Else no new processor was added as the requested one
                    % already existed
                    if size(initProc,2)==2
                        out{1,1} = initProc{1}.Output{1};
                        out{1,2} = initProc{2}.Output{1};
                    elseif size(initProc{1}.Output,2) == 2
                        out{1,1} = initProc{1}.Output{1};
                        out{1,2} = initProc{1}.Output{2};
                    else
                        out{1} = initProc{1}.Output{1};
                    end
                end
            elseif ~proceed
                % The request was invalid, return a empty handle
                out = [];
                
                % And remove the processors added by mistake
                if ~isempty(mObj.Processors{n_proc+1})
                    mObj.Processors{n_proc+1}.remove;
                    mObj.cleanup;
                end
            end
            
        end
        
        function cleanup(mObj)
            %CLEANUP  Clears the list of processors from handles to deleted processors
            
            %N.B.: We cannot use cellfun here as some elements of the .Processors array
            %are empty (e.g. when using binaural processors)
            
            % Loop through all elements to remove invalid handles
            for ii = 1:numel(mObj.Processors)
                if ~isempty(mObj.Processors{ii}) && ~isvalid(mObj.Processors{ii})
                    mObj.Processors{ii} = [];
                end
            end
            
            % Removes whole lines of empty elements from the list
            mObj.Processors( all( cellfun( @isempty, mObj.Processors), 2), : ) = [];
            
        end
        
        function reset(mObj)
            %reset  Resets the internal states of all instantiated processors
            %
            %USAGE:
            %  mObj.reset
            %
            %INPUT ARGUMENTS
            %  mObj : Manager instance
            
            % Is the manager working on a binaural signal?
            if size(mObj.Processors,2)==2
                
                % Then loop over the processors
                for ii = 1:size(mObj.Processors,1)
                   
                    % There should always be a processor for left/mono
                    mObj.Processors{ii,1}.reset;
                    
                    % Though there might not be a right-channel processor
                    if isa(mObj.Processors{ii,2},'Processor')
                        mObj.Processors{ii,2}.reset;
                    end
                        
                end
                
            else
            
                % Loop over the processors
                for ii = 1:size(mObj.Processors,1)
                    
                    mObj.Processors{ii,1}.reset;
                        
                end
            end
            
        end
        
    end
    
    methods (Access = protected)
       
        function [hProc,list] = findInitProc(mObj,request,p)
            %findInitProc   Find an initial compatible processor for a new
            %               request
            %
            %USAGE:
            %         hProc = mObj.findInitProc(request,p)
            %  [hProc,list] = mObj.findInitProc(request,p)
            %
            %INPUT PARAMETERS
            %    mObj : Manager instance
            % request : Requested signal name
            %       p : Parameter structure associated to the request
            %
            %OUTPUT PARAMETERS
            %   hProc : Handle to the highest processor in the processing 
            %           chain that is compatible with the provided
            %           parameters. In case two instances exist for the
            %           processor for a stereo signal, hProc is a cell
            %           array of the form {'leftEarProc','rightEarProc'}
            %    list : List of signal names that need to be computed,
            %           starting from the output of hProc, to obtain the
            %           request
        
            % Input parameter checking
%             if nargin<3 || isempty(p)
%                 % Initialize parameter structure
%                 p = struct;
%             end
%             if ~isfield(p,'fs')
%                 % Add sampling frequency to the parameter structure
%                 p.fs = mObj.Data.input{1}.FsHz;
%             end
%             % Add default values for parameters not explicitly defined in p
%             p = parseParameters(p);
%         
            % Try/Catch to check that the request is valid
%             try
%                 getDependencies(request);
%             catch err
%                 % Buid a list of available signals for display
%                 list = getDependencies('available');
%                 str = [];
%                 for ii = 1:size(list,2)-1
%                     str = [str list{ii} ', ']; %#ok
%                 end
%                 % Return the list
%                 error(['The requested signal, %s is unknown. '...
%                     'Valid names are as follows: %s'],request,str)
%             end
            
            % Get the full list of dependencies corresponding to the request
%             if ~strcmp(request,'time')
%                 dep_list = [request getDependencies(request)];
%             else
%                 % Time is a special case as it is listed as its own dependency
%                 dep_list = getDependencies(request);
%             end
            dep_list = [request ...
                Processor.getDependencyList( ...
                Processor.findProcessorFromRequest(request,p), p)];
            
            
            % Initialization of while loop
            ii = 1;
%             dep = signal2procName(dep_list{ii},p);
            dep = Processor.findProcessorFromRequest(dep_list{ii},p);
            hProc = mObj.hasProcessor(dep,p);
            list = {};
            
            % Looping until we find a suitable processor in the list of
            % dependency
            while hProc == 0 && ii<size(dep_list,2)
                
                % Then we will need to re-compute that signal
                list = [list dep_list{ii}]; %#ok
                
                % Move on to next level of dependency
                ii = ii + 1;
%                 dep = signal2procName(dep_list{ii},p);
                dep = Processor.findProcessorFromRequest(dep_list{ii},p);
                hProc = mObj.hasProcessor(dep,p);
                
            end
            
            if hProc == 0
                % Then all the signals need recomputation, including time
                list = [list dep_list{end}];
                
                % Return a empty handle
                hProc = [];
            end
            
            % If the processor found operates on the left channel of a stereo
            % signal, we need to find its twin processor in charge of the
            % right channel
            if ~isempty(hProc) && numel(hProc.Output) == 1 && ...
                    strcmp(hProc.Output{1}.Channel,'left')
                
                % Then repeat the same loop, but specifying the "other"
                % channel
                Channel = 'right';
                
                % Initialization of while loop
                ii = 1;
                dep = Processor.findProcessorFromRequest(dep_list{ii},p);
%                 dep = signal2procName(dep_list{ii},p);
                hProc2 = mObj.hasProcessor(dep,p,Channel);
                list = {};

                % Looping until we find a suitable processor in the list of
                % dependency
                while hProc2 == 0 && ii<size(dep_list,2)

                    % Then we will need to re-compute that signal
                    list = [list dep_list{ii}];     %#ok

                    % Move on to next level of dependency
                    ii = ii + 1;
%                     dep = signal2procName(dep_list{ii},p);
                    dep = Processor.findProcessorFromRequest(dep_list{ii},p);
                    hProc2 = mObj.hasProcessor(dep,p,Channel);

                end
                
                % Quick check that both found processor have the same task
                % (else there was probably an issue somewhere in channel
                % attribution)
%                 if ~strcmp(class(hProc),class(hProc2))
%                     error('Found different processors for left and right channels.')
%                 end
                
                % Put results in a cell array
                hProc = {hProc hProc2};
                
            else
                if ~isempty(hProc)
                    hProc = {hProc};
                end
            end
            
        end
        
        function newProcessor = addSingleProcessor(mObj,procName,parameters, ...
                                                dependencies,channelNb,index,channelTag)
            %addSingleProcessor     Instantiates a new processor and integrates it to the
            %                       manager instance
            %
            % Note about channelNb: 1 for left or mono, 2 for right. 
            %                       channelTag is 'mono' or 'stereo'
            %
            % The following steps are carried out:
            %   - Instantiate a processor, add a pointer to it in mObj.Processors
            %   - Generate a mutual link to its dependency
            %   - Instantiate a new output signal (possibly multiple)
            %   - Link it/them as output(s) of the processor
            %   - Provide link to the input signal(s)
            
            
            % TODO: test if index is necessary (to make use of preallocation), remove else
            if nargin<5||isempty(index)
                index = size(mObj.Processors,1)+1;
            end
            
            % Instantiate processor and add it to the list
            newProcessor = feval(procName, dependencies{1}.FsHzOut, parameters);
            mObj.Processors{index,channelNb} = newProcessor;
            
            % Labeling channels
            % TODO: Could be more flexible, to allow e.g., multi-channel processors
            if strcmp(channelTag,'mono')
                newProcessor.Channel = 'mono';
            else
                if channelNb == 1
                    newProcessor.Channel = 'left';
                else
                    newProcessor.Channel = 'right';
                end
            end
            
            % Mutual link to dependencies, unless it has none
            % IVO Comment: Maybe get the linked list of processors outside of the
            % processors, e.g. in a separate object/class
            if ~strcmp(newProcessor.getDependency,'input')
                newProcessor.addLowerDependencies(dependencies);
            end
            
            % Finalize processor initialization
            newProcessor.prepareForProcessing;
            
            % Instantiate and integrate new output signal
            output = newProcessor.instantiateOutput(mObj.Data);
            newProcessor.addOutput(output);
            
            % Integrate input signal pointer
            newProcessor.addInput(dependencies);
            
        end
            
            
            
        
    end
    
end