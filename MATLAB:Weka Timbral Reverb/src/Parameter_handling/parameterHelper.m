function  parameterHelper(procName)
%parameterHelper     Extensive and user friendly listing of parameters involved
%                    in the Two!Ears Auditory Front-End.
%
%USAGE:
%    parameterHelper



% Load the parameter info file
% path = fileparts(mfilename('fullpath'));
% load([path filesep 'parameterInfo.mat'])
% 
% % List the categories
% cats = fieldnames(pInfo);
% cats = sort(cats);

if nargin == 0 || isempty(procName)
    % Display the list of processors with an introductory header
    
    % Get a list of presumably valid processors
    procList = Processor.processorList;
    
    % Access their description
    procDescription = cell(size(procList));
    for ii = 1:size(procList,1)
        try
            pInfo = feval([procList{ii} '.getProcessorInfo']);
			procDescription{ii} = pInfo.label;
        catch
            procDescription{ii} = [];
        end
    end
    
    % Remove invalid elements
    procList = procList(~cellfun('isempty',procDescription));
    procDescription = procDescription(~cellfun('isempty',procDescription));
    
    % Get an index for alphabetical ordering
    [~,idx] = sort(procDescription);
    
    
    % Display header
    fprintf('\nParameter handling in the Two!Ears Auditory Front-End')
    fprintf('\n-------------------------------------------------\n')
    fprintf(['The extraction of various auditory representations '...
        'performed by the Two!Ears auditory front-end (AFE) software involves many parameters.\n'])
    fprintf('Each parameter is given a unique name and a default value. ')
    fprintf(['When placing a request for Two!Ears auditory front-end processing that\n'...
        'uses one or more non-default parameters, a specific structure of non-default parameters needs to be provided as input.\n'])
    fprintf('Such structure can be generated from <a href="matlab: help genParStruct">genParStruct</a>, using pairs of parameter name and chosen value as inputs.\n' )
    fprintf('\nParameters names for each processors are listed below:\n')

    % Display the categories
%     for ii = 1:size(cats,1)
%         category = cats{ii};
%         link = ['<a href="matlab:parameterHelper(''' category ''')">'];
%         fprintf(['\t' link pInfo.(cats{ii}).label '</a>\n'])
%     end
    
    % Display the categories
    for ii = 1:size(idx,1)
        hLink = ['<a href="matlab:parameterHelper(''' procList{idx(ii)} ''')">'];
        fprintf(['\t' hLink procDescription{idx(ii)} '</a>\n'])
    end

    % Add a plotting parameter category
    hLink = '<a href="matlab:parameterHelper(''plotting'')">';
    fprintf(['\t' hLink 'Plotting parameters' '</a>\n'])
    
    fprintf('\n')
else
    
    % Load the parameter names, default values and description 
    if ~strcmp(procName,'plotting')
        try 
            [names,values,description] = feval([procName '.getParameterInfo']);
        catch
            warning(['There is no %s processor, or its getParameterInfo static method '...
                     'is not implemented.'])
        end
    else
        % General plotting properties
        [names,values,description] = feval('Signal.getPlottingParameterInfo');
        
        % Specific properties
        sigList = Signal.signalList;
        
        for ii = 1:size(sigList,1)
            [names_sp,values_sp,description_sp] = ...
                feval([sigList{ii} '.getPlottingParameterInfo']);
            
            % Add to the list if non-empty
            if ~isempty(names_sp)
                names = [names, names_sp];                      %#ok<AGROW>
                values = [values, values_sp];                   %#ok<AGROW>
                description = [description, description_sp];    %#ok<AGROW>
            end
            
        end
        
    end
        
    % TODO: Do something specific for processors without parameters?
    
    % Find appropriate columns widths
    name_size = 4;  % Size of string 'Name'
    for ii = 1:size(names,2)
        name_size = max(name_size,size(names{ii},2));
    end
    
    valueStr = cell(1,size(names,2));
    val_size = 7;   % Size of string 'Default'
    
    % Text formatting for the parameter default value
    for ii = 1:size(names,2)
        if iscell(values{ii})
            % Then it's multiple strings concatenated in a cell array, add braces
            val = '{';
            for jj = 1:size(values{ii},2)-1  % TODO: Use strjoin instead?
                val = [val '''' values{ii}{jj} ''',']; %#ok<AGROW>
            end
            valueStr{ii} = [val '''' values{ii}{jj+1} '''}'];
        elseif ischar(values{ii})
            % Then it's a single string
            valueStr{ii} = ['''' values{ii} ''''];
        elseif size(values{ii},2)>1
            % Then it's an numerical array
            val = ['['];
            for jj = 1:size(values{ii},2)-1
                val = [val num2str(values{ii}(jj)) ' '];
            end
            valueStr{ii} = [val num2str(values{ii}(jj+1)) ']'];
        else
            % It's a single numeral
            valueStr{ii} = num2str(values{ii});
        end

        % Keep track of larger string for table formatting
        val_size = max(val_size,size(valueStr{ii},2));
    end
    
    % Display the processor name
    pInfo = feval([procName '.getProcessorInfo']);
    if ~isempty(names)
        fprintf([pInfo.label ' parameters:\n\n'])
    else
        fprintf([pInfo.label ' has no parameters.\n\n'])
    end
    
    % Display a header
    if ~isempty(names)
        fprintf(['  %-' int2str(name_size+2) 's  %-' int2str(val_size+1) 's  %-s\n'],'Name','Default','Description')
        fprintf(['  %-' int2str(name_size+2) 's  %-' int2str(val_size+1) 's  %-s\n'],'----','-------','-----------')
    end

    for ii = 1:size(names,2)
        % Display on command window
        fprintf(['  %-' int2str(name_size+2) 's  %-' int2str(val_size+1) 's  %-s\n'],...
                names{ii},valueStr{ii},description{ii})
    end
    fprintf('\n')
    
    
end
    
