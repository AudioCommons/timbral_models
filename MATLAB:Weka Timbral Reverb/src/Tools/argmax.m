function maxidx = argmax(input, dim)

if nargin < 2 || isempty(dim)
    [temp, maxidx] = max(input);
else
    [temp, maxidx] = max(input,[],dim);
end
