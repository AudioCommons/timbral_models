classdef Signal < matlab.mixin.Copyable
%SIGNAL Superclass for all signals involved in the auditory front-end (AFE) framework.
%   All signals involved in the AFE are inheriting this class, which defines properties
%   and methods shared among all signals.
%
%   SIGNAL properties:
%       Label      - Short and plain description of the signal (used e.g., as plot titles)
%       Name       - Single-word nametag for the signal
%       Dimensions - String describing signal's dimensions
%       FsHz       - Sampling frequency of the signal (in Hz)
%       Channel    - String indicating which audio channel the signal corresponds to
%
%   SIGNAL abstract method:
%       plot - Plots the signal. Should be implemented by any children class.
%
%   SIGNAL methods:
%       Signal         - Super-constructor for all signal objects
%       appendChunk    - Append a chunk of data to a signal buffer
%       setData        - Initialize the signal buffer
%       clearData      - Clears a signal buffer
%       getSignalBlock - Returns a data block from signal buffer
%       findProcessor  - Finds the processor that computed the signal
%       getParameters  - Returns the parameters used to compute the signal
%
% See also signals (folder), circVBuf, circVBufArrayInterface

    properties (SetAccess=protected)
        Label           % Used to label the signal (e.g., in plots)
        Name            % Used as an instance name in Matlab
        Dimensions      % String describing the dimensions of the signal
        FsHz            % Sampling frequency
        Channel         % Flag keeping track of the channel: 'mono', 'left'
                        %   or 'right'
    end
    
    properties (SetAccess = protected)
        Data;
    end

    properties (Access = protected)
        Buf;
    end

    methods (Abstract = true)
        h = plot(sObj,h_prev)
            % This method should create a figure displaying the data in the
            % signal object. As data can be of different nature for
            % different signal children objects, this method remains
            % abstract here.
               
    end
    
    methods
        
        function sObj = Signal( procHandle, bufferSize_s, bufferElemSize )
            %Signal     Super-constructor for the signal class
            %
            %USAGE:
            %   sObj = Signal(procHandle,bufferSize_s,bufferElemSize)
            %
            %INPUT ARGUMENTS:
            %     procHandle : Handle to the processor generating this signal
            %   bufferSize_s : Buffer duration in s
            % bufferElemSize : Additional dimensions of the buffer
            %                  [dim2,dim3,...]
            %
            %OUTPUT ARGUMENT:
            %           sObj : Signal instance
            
            % Set up sampling frequency
            sObj.FsHz = procHandle.FsHzOut;
            
            % Get the buffer size in samples
            bufferSizeSamples = ceil( bufferSize_s * sObj.FsHz );
            
            % Instantiate a buffer, and an array interface
            sObj.Buf = circVBuf( bufferSizeSamples, bufferElemSize );
            sObj.Data = circVBufArrayInterface( sObj.Buf );
            
            % Populate name and label
            sObj.Name = procHandle.getProcessorInfo.requestName;
            sObj.Label = procHandle.getProcessorInfo.requestLabel;
        end

        function setBufferSize( sObj, newBufferSize_s )
            %setBufferSize  This method sets the buffer to a new size,
            %               erasing all data previously stored.
            %USAGE:
            %   sObj.setBufferSize( newBufferSize_s )
            %
            %INPUT ARGUMENTS:
            %   newBufferSize_s : new size of the buffer in seconds. The
            %                     dimensionality of the individual elements
            %                     remains the same
            
            bufferSizeSamples = ceil( newBufferSize_s * sObj.FsHz );
            sObj.Buf = circVBuf( bufferSizeSamples, sObj.Buf.matSz );
            sObj.Data = circVBufArrayInterface( sObj.Buf );
        end
        
        function appendChunk(sObj,data)
            %appendChunk   This method appends the chunk in data to the signal object
            %
            %USAGE
            %   sObj.appendChunk(data)
            %
            %INPUTS
            %   sObj : Signal object
            %   data : New chunk (re. time) to be appended to existing data
            %
            %DISCLAIMER:
            %This method is likely called numerous times hence there is no
            %checking on the provided data chunk (e.g., regarding its
            %dimensionality) to improve efficiency. Although for practical
            %reason it remains public, it should not be called "manually"
            %by a user.
            
            sObj.Buf.append( data );
            
        end

        function setData(sObj, data)
            %setData    Initialize the cyclic buffer using provided data
            %
            %USAGE:
            %   sObj.setData(data)
            %
            %INPUT ARGUMENTS:
            %   sObj : Signal instance
            %   data : Block of data to initialize the signal with
            
            % First: clear the buffer
            sObj.Buf.clear();

            % Then append the provided data
            if ~isempty(data), sObj.Buf.append( data ); end;

        end
        
        function clearData(sObj)
            %clearData  Clears the data in a signal object without changing
            %           its other properties
            %
            %USAGE
            %   sObj.clearData
            
            sObj.Buf.clear();
        end
        
        function newSobj = cutSignalCopy( sObj, blocksize_s, backOffset_s )
            %cutSignalCopy  This method copies the Signal object into a new
            %               instance, cutting out the specified data block.
            %
            %USAGE:
            %   cutSignalCopy = sObj.cutSignalCopy( blockSize_s, backOffset_s )
            %
            %INPUT ARGUMENTS:
            %  blocksize_s : Length of the required data block in seconds
            % backOffset_s : Offset from the end of the signal to the 
            %                requested block's end in seconds (default: 0s)
            newSobj = sObj.copy();
            newSobj.setBufferSize( ceil( blocksize_s ) );
            dataBlock = sObj.getSignalBlock( blocksize_s, backOffset_s );
            newSobj.setData( dataBlock );
        end
        
        function dataBlockResampled = resampleToFsHz( sObj, dataBlock, srcFsHz )
            %resampleToFsHz  This method resamples a data block
            %               to match the signal's sampling rate.
            %               Resampling is perforrmed via linear
            %               interpolation without decimation.
            %
            %USAGE:
            %   resample = sObj.resampleToFsHz( dataBlock, srcFsHz )
            %
            %INPUT ARGUMENTS:
            %  dataBlock : The data block to resample with signal's FsHz
            %  srcFsHz : the original sampling rate of the data block
            %  
            [rows, ~] = size(dataBlock);
            x = 0 : 1 / srcFsHz : (rows-1) / srcFsHz;
            xq = 0 : size( sObj.Data, 1 ) - 1;
            xq = xq * 1 / sObj.FsHz;
            x = x + (max( xq ) - max( x ));
            dataBlockResampled = interp1( x, dataBlock, xq, 'linear', 'extrap' );
        end
                
        function newSobj = maskSignalCopy( sObj, mask, freq_src, maskHopSize )
            %maskSignalCopy  mask a copy of a Signal's dataBlock
            %
            %USAGE:
            %   newSobj = sObj.maskSignalCopy( mask, maskHopSize )
            %
            %INPUT ARGUMENTS:
            %  mask : A 2-d mask to apply on the Signal's dataBlock
            %  maskHopSize : the mask's sampling rate (tmeporal resolution)
            %  
            if nargin < 3
                % assume both are sampled at equal rates
                maskHopSize = 1./sObj.FsHz;
            end
            newSobj = sObj.copy();
            dataBlock = newSobj.Data(:,:,:,:,:);
            [rd, ~, dd] = size(dataBlock);
            if ((1/maskHopSize) ~= sObj.FsHz) || (size( mask, 1 ) ~= rd)
                mask = sObj.resampleToFsHz( mask, 1./maskHopSize );
            end
            [rm, ~, dm] = size(mask);
            if rm > rd
                % crop
                mask = mask(end+1-max(1, rd):end, :);
            elseif rm < rd
                error('mask too short for dataBlock');
            end
            if any(size( freq_src ) ~= size( sObj.cfHz )) || any(freq_src ~= sObj.cfHz)
                mask = interp1( freq_src, mask', sObj.cfHz, 'linear', 'extrap' )';
            end
            if dd > dm
                mask = repmat(mask, [1, 1, dd]);
            elseif dd < dm
                error('cannot mask dataBlock dimensions < mask');
            end
            dataBlock = dataBlock .* mask;
            if isa(newSobj.Buf, 'circVBuf')
                newSobj.setData( dataBlock );
            else
                newSobj.Data = dataBlock;
            end
        end
        
        function reduceBufferToArray( sObj )
            %reduceBufferToArray    This method converts the
            %                       buffer+interface combination into a
            %                       native matlab array. The object should
            %                       not be used for adding new data any
            %                       more.
            %
            %USAGE:
            %   sObj.reduceBufferToArray
            %
            data = sObj.Data(:);
            delete( sObj.Data );
            sObj.Data = data;
            delete( sObj.Buf );
        end
        
        function newSobj = cutSignalCopyReducedToArray( sObj, blocksize_s, backOffset_s )
            %cutSignalCopyReducedToArray  This method copies the Signal object into a new
            %                             instance, cutting out the specified data block
            %                             and reducing it to a matlab array.
            %
            %USAGE:
            %   cutSignalCopy = sObj.cutSignalCopyReducedToArray( blockSize_s, backOffset_s )
            %
            %INPUT ARGUMENTS:
            %  blocksize_s : Length of the required data block in seconds
            % backOffset_s : Offset from the end of the signal to the 
            %                requested block's end in seconds (default: 0s)
            Signal.doShallowCopy( true, true );
            newSobj = sObj.copy();
            Signal.doShallowCopy( true, false );
            newSobj.Buf = [];
            newSobj.Data = sObj.getSignalBlock( blocksize_s, backOffset_s );
        end
        
        function sb = getSignalBlock(sObj,blocksize_s,backOffset_s,padFront)
            %getSignalBlock   Returns this Signal object's signal data
            %truncated to the last blocksize_s seconds. In case of too
            %little data, the block can get filled with zeros from beginning.
            %
            %USAGE:
            %   sb = sObj.getSignalBlock(blocksize_s)
            %   sb = sObj.getSignalBlock(blocksize_s,backOffset_s)
            %
            %INPUT ARGUMENTS:
            %         sObj : Signal instance
            %  blocksize_s : Length of the required data block in seconds
            % backOffset_s : Offset from the end of the signal to the 
            %                requested block's end in seconds (default: 0s)
            %     padFront : whether to pad with zeros or not. Default: true
            %
            %OUTPUT ARGUMENTS:
            %    sb : signal data block
            
            % Get the block duration in samples
            blocksize_samples = ceil( sObj.FsHz * blocksize_s );
            
            % Set default value for backOffset_s
            if nargin < 3, backOffset_s = 0; end;
            
            % Get the offset in samples...
            % floor it because the last few milliseconds usually are lost in the current
            % processchunk()
            % subtract 1 because the last frame is the last full window (before: shifts)
            offset_samples = max( 0, floor( sObj.FsHz * backOffset_s ) - 1 );
            
            % ... with a warning if the requested signal is "too old"
            if offset_samples >= size(sObj.Data,1)
                warning( ['You are requesting a block that is not in the ',...
                    'buffer anymore.'] );
            end
            
            % Figure out the starting index in the buffer
            blockStart = max( 1, size( sObj.Data, 1 ) - ...
                blocksize_samples - offset_samples + 1 );
            
            % Extract the data block
            sb = sObj.Data(blockStart:end-offset_samples,:,:,:,:,:,:);
            
            % Zero-pad the data if not enough samples are available
            if nargin < 4, padFront = true; end
            if padFront && (size( sb, 1 ) < blocksize_samples)
                sb = [zeros( blocksize_samples - size(sb,1), size(sb,2), size(sb,3) ); sb];
            end
            
        end
            
        function pObj = findProcessor(sObj,mObj)
            %findProcessor   Returns a handle to the processor instance
            %that generated this signal.
            %
            %USAGE:
            %    pObj = sObj.findProcessor(mObj)
            %
            %INPUT ARGUMENTS:
            %    sObj : Signal instance
            %    mObj : Manager instance containing the sought processor
            %
            %OUTPUT ARGUMENTS:
            %    pObj : Handle to the processor instance which computed the
            %           signal (empty if none found)
            
            % NB: This brute force approach could be made more subtle by
            % looking into the type of signal sObj is. However it shouldn't
            % be necessary.
            
            % Number of instantiated processors
            n = numel(mObj.Processors);
            
            % Initialize output
            pObj = [];
            
            % Loop over all of them
            for ii = 1:n
                % Check that it is actually a processors
                if isa(mObj.Processors{ii},'Processor')
                    % Check if it outputs the signal of interest
                    for jj = 1:size(mObj.Processors{ii}.Output,2)
                        if sObj == mObj.Processors{ii}.Output{jj}
                            pObj = mObj.Processors{ii};
                        end
                    end
                end
            end 
        end
        
        function parStruct = getParameters(sObj,mObj)
            %getParameters  This methods returns a list of parameter
            %values used to compute a given signal.
            %
            %USAGE:
            %   parStruct = sObj.getParameters(mObj)
            %
            %INPUT PARAMETERS:
            %   sObj : Signal object instance
            %   mObj : Manager instance containing the processor responsible
            %   for computing the signal of interest
            %
            %OUTPUT PARAMETERS:
            %   parStruct : Parameter structure
            
            % Find the processor that computed sObj
            proc = sObj.findProcessor(mObj);
            
            % Get the parameters under which this processor was running
            if ~isempty(proc)
                parStruct = proc.getCurrentParameters;
            else
                % Couldn't find the processor in charge
                parStruct = struct;
                % Return a warning, unless sObj is the original ear signal
                if ~strcmp(sObj.Name,'input')
                    warning('Could not find the processor that computed the signal ''%s.''',sObj.Name)
                end
            end
            
        end
        
    end
    
    methods (Access = protected)

        function cpObj = copyElement(obj)
            %copyElement	Override of Copyable method. Needed because Buf is handle.
            %               If Subclasses of Signal have private properties, override
            %               copyElement in those classes!

            cpObj = copyElement@matlab.mixin.Copyable(obj);
            if ~Signal.doShallowCopy
                if isa( cpObj.Buf, 'circVBuf' ) && cpObj.Buf.isvalid()
                    cpObj.setBufferSize( ceil( size( obj.Buf.dat, 1 ) / obj.FsHz ) );
                    cpObj.setData( obj.Data(:) );
                end
            end
        end
      
    end
    
    methods (Static)
        
        function b = doShallowCopy( bSet, newValue )
            persistent dsc;
            if isempty( dsc )
                dsc = false;
            end
            if nargin > 0  &&  bSet
                dsc = newValue;
            end
            b = dsc;
        end

        
        function sList = signalList()
            
            % Signal directory
            signalDir = mfilename('fullpath');
            
            % Get file information
            fileList = listFiles(signalDir(1:end-7),'*.m',-1);
            
            % Extract name only
            sList = cell(size(fileList));
            for ii = 1:size(fileList)
                % Get file name
                [~,fName] = fileparts(fileList(ii).name);
                
                % Check if it is a valid signal
                try
                    s = feval(str2func(fName));
                    if isa(s,'Signal')
                        sList{ii} = fName;
                    else
                        sList{ii} = [];
                    end
                catch
                    sList{ii} = [];
                end
                
            end
                
            % Remove empty elements
            sList = sList(~cellfun('isempty',sList));
            
        end
       
        function signalName = findSignalFromParameter(parameterName,no_warning)
            %Signal.findSignalFromParameter   Finds the signal that uses a given plotting
            %parameter
            %
            %USAGE:
            %   signalName = Signal.findSignalFromParameter(parName)
            %   signalName = Signal.findSignalFromParameter(parName, no_warning)
            %
            %INPUT ARGUMENT:
            %      parName : Name of the parameter
            %   no_warning : Set to 1 (default: 0) to suppress warning message
            %
            %OUTPUT ARGUMENT:
            %   signalName : Name of the signal using that parameter

            if nargin<2||isempty(no_warning); no_warning = 0; end
            
            % Get a list of processor
            signalList = Signal.signalList;

            % Add the general plot properties to the list
            signalList = ['Signal'; signalList];
            
            % Loop over each processor
            for ii = 1:size(signalList,1)
                try
                    sigParNames = feval([signalList{ii} '.getPlottingParameterInfo']);
                    
                    if ismember(parameterName,sigParNames)
                        signalName = signalList{ii};
                        return
                    end
                    
                catch
                    % Do not return a warning here, as this is called in a loop
                end

            end

            % If still running, then we haven't found it
            if ~no_warning
                warning('Could not find a signal which uses plotting parameter ''%s''',...
                        parameterName)
            end
            signalName = [];
            
            
        end
        
        function [names, defaultValues, descriptions] = getPlottingParameterInfo()
            %GETPLOTTINGPARAMETERINFO   Stores plot parameters that are common to all
            %signals.
            
            
            names = {'ftype',...
                    'fsize_label',...
                    'fsize_title',...
                    'fsize_axes',...
                    'color',...
                    'colors',...
                    'linewidth_s',...
                    'linewidth_m',...
                    'linewidth_l'};
                 
            descriptions = {'Plots font name', ...
                    'Labels font size', ...
                    'Titles font size', ...
                    'Axes font size', ...
                    'Main plot color', ...
                    'Multiple plot colors', ...
                    'Small linewidth', ...
                    'Medium linewidth', ...
                    'Large linewidth'};
                
            defaultValues = {'Helvetica', ...
                    12, ...
                    14, ...
                    10, ...
                    'b', ...
                    {'b', 'r', 'g', 'c'}, ...
                    1, ...
                    2, ...
                    3};
            
        end
        
    end
    
end