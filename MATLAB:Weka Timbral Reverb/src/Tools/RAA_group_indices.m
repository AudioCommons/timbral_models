function [ consecIdxGroups ] = RAA_group_indices( idx_passing_threshold )
%UNTITLED Summary of this function goes here
%   Group consecutive indices found to be above or below a threshold

% Differentiate
diff_idx = diff(idx_passing_threshold);
% Find differences larger than 1 - this means a gap between indices
idxGap = find(diff_idx>1);

% Attach 0 and the last given index value for grouping
idxGap = [0; idxGap; length(idx_passing_threshold)];

consecIdxGroups = cell(length(idxGap)-1, 1);

for i=1:length(idxGap)-1
    % Return groups of consecutive indices
    consecIdxGroups{i} = idx_passing_threshold(idxGap(i)+1:idxGap(i+1));
end



end

