function bSelect = selectCells(fullList,exp)

if ~iscell(fullList)
   fullList = {fullList}; 
end
if ~iscell(exp)
   exp = {exp}; 
end

% Number of expressions
nExp = numel(exp);

% Allocate memory
bSelect = false(size(fullList));

% Loop over number of entries
for ii = 1 : nExp
    bExp = strcmp(fullList,exp{ii});
    
    bSelect = xor(bSelect,bExp);
end