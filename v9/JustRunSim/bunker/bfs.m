function [Astruct,nclusts,N_G] = bfs(A,opts)
%BFS1 Breadth-First-Search for number of clusters in a graph

% A is the connectivity matrix
% formated so a node row-wise recieves input from those columns
% a node column-wise sends inputs to those rows.

% get number of nodes
N = length(A);

% set up structures
Astruct.A = A;
% set up structures for components inside matrix
Astruct.component = [];

% initialize B, for keeping track of total nodes counted
B = zeros(N,1); 
% initialize c, for number of actual components
c=0;

% while any node of B hasn't been counted
while any(B==0)
    % update component counter
    c = c+1;
    % pick an uncounted node
    l = find(B==0);
    li = l(randperm(length(l),1));
    
    % update B
    B(li) = 1;
    
    % reset n_nodes
    n_nodes = li;
    
    % initialize component structure bfs array
    Astruct.component(c).bfs = zeros(N,1);
    % set first node to 1
    Astruct.component(c).bfs(li) = 1;
    % begin bfs counter
    q = 1;
    
    % repeat until no new nodes are in component
    while q<N^2
        % update bfs counter
        q = q+1;      
        
        % get rows (outputs) of this iters nodes
        [a1,~] = find(A(:,n_nodes));
        % remove redundencies
        a2 = unique(a1);
        
        % update B, set changed char to no
        B(a2) = 1; ischanged = 0;
        
        % loop through nodes
        for i = 1:length(a2)
            % if the node is not yet recorded
           if Astruct.component(c).bfs(a2(i)) == 0
               % apply the bfs counter to it
               Astruct.component(c).bfs(a2(i)) = q;
               % dont break the loop yet
               ischanged = 1;
           end
        end
        
        % reset the nodes for what iter we are on
        n_nodes = a2;

        % if no new node was added this iter
        if ischanged == 0
            break
        end
        
        if q == N^2
            disp('woah')
        end
        
    end
    
end


% number of clusters
nclusts = length({Astruct.component.bfs});

% nodes in giant component
if nargout >= 3
    N_C = zeros(c,1);
    
    for i = 1:c
        
        N_C(i) = sum(Astruct.component(i).bfs~=0);
      
    end
    
    N_G = max(N_C);
end


end

