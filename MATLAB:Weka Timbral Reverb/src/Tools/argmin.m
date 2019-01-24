function minidx = argmin(input, dim)

if nargin < 2 || isempty(dim)
    [temp, minidx] = min(input);
else
    [temp, minidx] = min(input,[],dim);
end
