function requestList
%requestList  Returns a list of currently implemented valid requests for 
%             Two!Ears Auditory Front-End processing.
%
%USAGE:
%   requestList()

% Get a list of processor
procList = Processor.processorList;

list = cell(size(procList,1),3);

l_req = 12;
l_lab = 5;
l_prc = 0;

for ii = 1:size(list,1)
    pInfo = feval([procList{ii} '.getProcessorInfo']);
    list{ii,1} = pInfo.requestName;
    list{ii,2} = pInfo.requestLabel;
    list{ii,3} = procList{ii};
    % Getting column width for display
    l_req = max(size(pInfo.requestName,2),l_req);
    l_lab = max(size(pInfo.requestLabel,2),l_lab);
    l_prc = max(size(procList{ii},2),l_prc);
end

% Sort alphabetically
[~,I] = sort(list(:,1));
list = list(I,:);

% Displaying it in three columns
fprintf('\n')
fprintf(['  %-' int2str(l_req+2) 's  %-' int2str(l_lab+1) 's  %-s\n'],...
        'Request name','Label','Corresponding processor')
fprintf(['  %-' int2str(l_req+2) 's  %-' int2str(l_lab+1) 's  %-s\n'],...
        '------------','-----','-----------------------')    

for ii = 1:size(list,1)
    % Display on command window
    fprintf(['  %-' int2str(l_req+2) 's  %-' int2str(l_lab+1) 's  %-s\n'],...
            list{ii,1},list{ii,2},list{ii,3})
end
fprintf('\n')

