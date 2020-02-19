function [A] = MakeNetwork(Parameters,Options)
%MAKENETWORK
disp(['Making network. simtime (s) = ' mat2str(round(toc,2)) ])

N = Parameters.n1 + Parameters.n2;

% assign variables
switch Options.conmat
    case 1
        prob = Parameters.er_prob;
    case 2
        sw_p = Parameters.sw_p;
        sw_M = Parameters.sw_M;
    case 3
        sf_mo = Parameters.sf_mo;
        sf_m = Parameters.sf_m;
        sf_d = Parameters.sf_d;
    case 5
        ke_mo = Parameters.ke_mo;
        ke_mu = Parameters.ke_mu;
        ke_d = Parameters.ke_d;
end
% check if network already exists,
if isfield(Parameters,'A')
    disp('I see Parameters.A exists, therefore we will load that loaded network')
    disp('Not remaking a network!!! Note to user! Parameters.conmat is being ignored!')
    A = Parameters.A;
else

switch Options.conmat
    case 1 %E-R network
        
        A=(rand(N,N)<prob);
        for i=1:N, A(i,i)=0; end %makes sure no neuron is connected to itself
        
    case 2 %small world
        
        A = makesmallworld(N,sw_p, sw_M);
        
    case 3 %scale free
        
        A = makeBAscalefree(N,sf_mo, sf_m,sf_d,Options);
        
    case 4 %load scale-free
        
        load(Parameters.A2load)
        
        %     case 5 %directed clique
        %         A = tril(ones(N), -1);
        %         A = A';
        
    case 5 %Klemm and Eguilez
        
        A = makeKEscalefree(N,ke_mo,ke_mu,ke_d,Options);
        
    case 6 % GH lattice
        
        A = GHLattice(Parameters.n1+Parameters.n2,Parameters.s);
        
end

end

%% connmatrix visualization %
% check if option exists,
if ~isfield(Options,'doAplot')
    Options.doAplot = 0;
end

if Options.doAplot
    figure
    try
        plot(graph(A))
    catch
        plot(digraph(A))
    end
    title('network A', 'FontSize', 15)
    figure
    spy(A)
    title('Adjacency matrix', 'FontSize', 15)
end

disp(['Done network making. simtime (s) = ' mat2str(round(toc,2)) ])

%% Connectivity Matrix functions %%

%% small world %

    function A = makesmallworld(N,sw_p,sw_M)
        A = gallery('circul', N);
        %sw_X = sw_p*N; %number of "short cuts" in small world network
        
        val = zeros(1,N);
        val(2:sw_M+1) = 1;
        val(N-sw_M+1:N) = 1;
        A = val(A);
        % somehow that made the circular matrix w/ correct number of
        % neighbors
        
        Apicture = A;
        [row,col] = find(A);
        mutations = 0;
        
        for s = 1:length(row)
            if row(s) ~= col(s)
                row(s);
                col(s);
                if rand<=sw_p
                    mutations = mutations+1;
                    a = find(A(row(s),:)==0);
                    a(a==row(s))=[];
                    
                    XX = a(randperm(length(a),1));
                    
                    A(row(s),col(s)) = 0;
                    A(row(s),XX) = 1;
                end
                
            end
        end
        
        
        %         for sw_i=1:sw_X
        %             row = randi(N);
        %             col = find(A(row,:)==0);
        %             %col = randi(N);
        %             %a = find(A(row,:));
        %             A(row,a(randi(length(a)))) = 0;
        %             if row ~= col
        %                 A(row, col) = 1;
        %                 A(col, row) = 1;
        %             end
        %         end
        
    end

%% scale free Barabasi and Albert

    function A = makeBAscalefree(N,sf_mo,sf_m,d,opts)
        % Makes a BA scale free network
        % N is number of nodes in network
        % sf_mo is number of nodes to start with in network
        % sf_m is number of edges to add each time step
        
        % opts.startstruct deteremines starting structure of nodes
        % opts.startstruct = 0 - fully conected clique (default)
        % opts.startstruct = 1 - randomly connected (ER)
        % opts.startstruct.p = p of ER network (only if startstruct = 1)
        
        
        
        opts.forcespawn = [];
        if ~isfield(opts,'startstruct')
            opts.startstruct.s = 0;
            disp('going with default starting network (clique)')
        end
        
        % initialize matrix
        A = zeros(N);
        % begin network with fully connected network of sf_mo
        
        
        % design initial network architecture
        if opts.startstruct.s == 0
            % clique
            
            A(1:sf_mo,1:sf_mo) = 1;
            for i = 1:sf_mo; A(i,i) = 0; end
            
            
        elseif opts.startstruct.s == 1
            % ER
            
            A(1:sf_mo,1:sf_mo)=(rand(sf_mo)<opts.startstruct.p);
            A(sf_mo+1:end,sf_mo+1:end) = 0;
            for i=1:N; A(i,i)=0; end
            
        else
            error('network starting structure not specified')
        end
        
        % edge counter
        E = length(find(A==1));
        
        
        % waitbar
        wb = waitbar(0,'building scale free network');
        
        % for each remaining node to add
        for i = sf_mo + 1: N
            % keep track of the degree
            Current_Degree = 0;
            % while the degree is less than m
            while Current_Degree<sf_m
                % find node j, uniformly randomly chosen from the set of
                % all nodes, excluding i and nodes adjacent to i (from i to
                % j)
                jpre = find(A(i,:)==0);
                jpre(jpre==i) = [];
                j = jpre(randperm(length(jpre),1));
                
                % set b1 = number of nodes adjacent to j / E
                % adjacent to j == j --> another node
                b1 = sum(A(j,:)) / E;
                
                % set Chance = a uniform random number between 0 and 1
                Chance1 = rand;
                
                if b1>Chance1
                    
                    % flip to see if connection is reciprical or not
                    Chance2 = rand;
                    if d >= Chance2
                        % reciprical
                        A(i,j) = 1;
                        A(j,i) = 1;
                        E = E+2;
                    else
                        % not recripical
                        A(i,j) = 1;
                        E = E+1;
                        
                        No_Connection = true;
                        while No_Connection
                            % set node h = uniformly randomly chosen node
                            % from the set of all nodes, excluding i and
                            % nodes adjacent to i, exluding j
                            hpre = find(A(:,i)==0);
                            hpre(hpre==i) = [];
                            %hpre(hpre==j) = [];
                            h = hpre(randperm(length(hpre),1));
                            
                            b2 = sum(A(h,:)) / E;
                            
                            Chance3 = rand;
                            if b2 > Chance3
                                A(h,i) = 1;
                                E = E + 1;
                                No_Connection = false;
                            end
                        end
                    end
                end
                Current_Degree = sum(A(i,:));
            end
            waitbar(i/N)
        end
        
        % figures?
        %         figure;
        %         chunks = sum(A,2);
        %         %histogram(chunks)
        %         [a,b] = histcounts(chunks);
        %         temp = movmean(b,2);
        %         loglog(temp(2:end),a,'x','markersize',16)
        %         xlabel('k degree')
        %         ylabel('count of nodes with degree k');
        
        
        
        close(wb)
        
        
        
        
        
        
    end

%% Klemm and Eguilez scale-free & small-world

    function A = makeKEscalefree(N,ke_mo,mu,d,opts)
        
        % opts.startstruct deteremines starting structure of nodes
        % opts.startstruct = 0 - fully conected clique (default)
        % opts.startstruct = 1 - randomly connected (ER)
        % opts.startstruct.p = p of ER network (only if startstruct = 1)
        
        
        opts.forcespawn = [];
        if ~isfield(opts,'startstruct')
            opts.startstruct.s = 0;
            disp('going with default starting network (clique)')
        end
        % initialize matrix
        A = zeros(N);
        % begin network with fully connected network of sf_mo
        % design initial network architecture
        if opts.startstruct.s == 0
            % clique
            A(1:ke_mo,1:ke_mo) = 1;
            for i = 1:ke_mo; A(i,i) = 0; end
        elseif opts.startstruct.s == 1
            % ER
            A(1:ke_mo,1:ke_mo)=(rand(ke_mo)<opts.startstruct.p);
            A(ke_mo+1:end,ke_mo+1:end) = 0;
            for i=1:N; A(i,i)=0; end
        else
            error('network starting structure not specified')
        end
        
        Active_Nodes = false(N,1);
        
        
        % first, create a fully connected initial network, and add its
        % nodes to active nodes.
        
        Active_Nodes(1:ke_mo) = true;
        Deactivated_Nodes = (~Active_Nodes);
        
        % waitbar
        wb = waitbar(0,'building KE network');
        
        % second, iteratively connect remaining nodes to active nodes or a
        % rangom node with probability mu, then remove one active node.
        
        for i = ke_mo + 1 : N
            waitbar(i/N)
            for j = find(Active_Nodes)'
                Chance = rand;
                % mu is the chance edges are added to non-active nodes
                if Chance > mu || i == ke_mo+1
                    A(i,j) = 1;
                    A(j,i) = 1;
                else
                    Connected = false;
                    while Connected == false
                        nodej_pre = find(Deactivated_Nodes);
                        nodej_pre(nodej_pre==i) = [];
                        nodej = datasample(nodej_pre,1);
                        
                        Chance = rand;
                        % E = sum((out) Degree of all deactivated nodes)
                        E = sum(A(nodej_pre,:),'all');
                        if E == 0
                            E = eps;
                        end
                            
                        
                        % kj = the degree of node k
                        kj = sum(A(nodej,:));
                        
                        if kj/E > Chance
                            A(i,nodej) = 1;
                            A(nodej,i) = 1;
                            Connected = true;
                        end
                    end
                end
                
            end
                
                % Replace an active node with node i. Active nodes with
                % lower degree are more likely to be replaced.
                
                Active_Nodes(i) = true;
                Deactivated_Nodes(i) = false;
                
                nodejchosen = false;
                while nodejchosen == false
                    nodej_pre = find(Active_Nodes);
                    nodej = datasample(nodej_pre,1);
                    
                    % kj = the degree of node k
                    kj = sum(A(nodej,:));
                    
                    % E = sum((out) Degree of all active nodes)
                    E = sum(A(Active_Nodes,:),'all');
                    
                    % pd = probability a node will be deactived
                    % % % moving away from paper psuedo-code here... %
                    % in paper: pd = (1/kj) / sum(1/kj)....
                    % how the fuck would that work? kj is not an array!
                    % it is the degree of a single node!
                    % i think what ive done works too.
                    % nodes will high degree have less chance of being
                    % picked...
                    pd = (1/kj) / E;
                    
                    
                    Chance = rand;
                    
                    if pd>Chance
                        nodejchosen = true;
                        % remove j from active nodes
                        Active_Nodes(nodej) = false;
                        Deactivated_Nodes(nodej) = true;
                    end
                    
                end
        end
        
        % Third, if directed, use d to set undirected equivalency ( i think
        % text is missing an endfor before this...?
        
        for i = 1:N
            for j = find(A(i,:))
                Chance = rand;
                if d > Chance
                    nodeh_pre = find(A(i,:)==0);
                    nodeh_pre(nodeh_pre==i) = [];
                    nodeh_pre(nodeh_pre==j) = [];
                    
                    if ~isempty(nodeh_pre)
                        nodeh = datasample(nodeh_pre,1);
                        
                        A(i,j) = 0;
                        A(i,nodeh) = 1;
                    end
                    
                    
                end
            end
        end
    close(wb)    
    end

%% GH Lattice

    function A = GHLattice(N,s)
        % a is the baby matrix
        a = zeros(N);
        % first we define the total grid size (sqrt(N))
        gridMAX = sqrt(N);  if floor(gridMAX)~=gridMAX; 
                            error('N must have an integer square root'); 
                            end
        
        % getting the array which houses the positions of the nodes ready
        xpos = zeros(gridMAX,gridMAX,2);

        xpos = repmat([1:gridMAX],[gridMAX,1]);
        ypos = xpos';
        
        for myn = 1:N
            
            myx = xpos(myn); myy = ypos(myn);
            
            for othern = 1:N
                if myn ~= othern
                    
                    otherx = xpos(othern); othery = ypos(othern);
                    
                    dx = abs(myx-otherx); dy = abs(myy-othery);
                    
                    d = sqrt(dx^2 + dy^2);
                    prob = exp(-d^2/(2*s^2));
                    
                    if prob>rand
                        % connection
                        a(myn,othern) = 1;
                    end
                    
                end
            end
        end
        
       
        A = a;
        
        givestats = 1;
        if givestats
            meankin = mean(sum(a))
            meankout = mean(sum(a,2))
            
        end
        
    end


end

