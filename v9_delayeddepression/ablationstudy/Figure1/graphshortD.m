function [md,md_nanclust] = graphshortD(A)
%GRAPHSHORTD This function returns the mean path distance between any (all) 
%two node combinations in directed graph A

% A is an adjacenecy matrix
% md is the mean distance between any two nodes
% md_nanclust means I dont look at path lengths where d is infinity.
% md can return nan if there is more than one cluster (lone nodes)


[~,catcher] = bfs(A);

if catcher ~= 1
    %error('There are two clusters.. might have to account for this. How do we calculate mean distance if one of the paths is infinity?')
end

N = length(A);

% ensure format
A = double(A);

% initialize funcction output
md_pre = zeros(N,1);
md_pre_nc = zeros(N,1);

% loop through each node...

for i = 1:length(A)
    
    % initialize distances array (distance from node i to all nodes)
    distances = zeros(N,1);
    % initialize array of uncounted nodes
    uncounted_nodes = 1:N;
    % remove node i from uncounted nodes, as that is itself.
    uncounted_nodes(i) = [];
    % increment all uncounted nodes (all but i)
    distances(uncounted_nodes) = distances(uncounted_nodes)+1;
    % initialize next_nodes parameter for loop
    next_nodes = i;
    
    iter = 0;
    % while there are nodes still uncounted...
    while ~isempty(uncounted_nodes) && iter < 1e3
        % determine which neurons the ones we have go to...
        [~,next_nodes_pre] = find(A(next_nodes,:));
        next_nodes = unique(next_nodes_pre);
        % removes those nodes from uncounted nodes
        for k = 1:length(next_nodes)
            uncounted_nodes(uncounted_nodes==next_nodes(k)) = [];
        end
        % incriment nodes again
        distances(uncounted_nodes) = distances(uncounted_nodes)+1; 
        
        iter = iter+1;
        if iter == 1e3
            distances(uncounted_nodes) = NaN;
            %warning('This network appears to have more than one cluster')
        end
    end
    
    md_pre(i) = mean(distances);
    md_pre_nc(i) = mean(distances,'omitnan');
    
end

md = mean(md_pre);
md_nanclust = mean(md_pre_nc,'omitnan');



end

