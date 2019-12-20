%% Algorithms 1 Exercise 4: Compute min cut for an undirected graph
%                           represented using an adjacency list.

% Initialization
clear ; close all; clc

% Read the input file 
%[V, L, n] = ReadAdjList('SimplerTest.txt');
%[V, L, n] = ReadAdjList('SimpleTest5.txt');
[V, L, n] = ReadAdjList('kargerMinCut.txt');
fprintf('Computing min cut of graph with %i vertices\n', n);

% Save a copy of the original adjacency list
origV = V;
origL = L;
origN = n;

%DisplayGraph(origV, origL, origN);

rng('shuffle');  % reset the random number generator with a new seed
%rng(314159258);
    
fprintf('Number of crossing edges:\n');
bestCrossingNumber = n*2;  % anything larger than n
for count = 1:5000
    V = origV;
    L = origL;
    n = origN;

    % Karger's random contraction algorithm for min cuts:
    % While there are more than 2 vertices:
    %   - pick a remaining edge (u,v) uniformly at random
    %   - merge (or ?contract? ) u and v into a single vertex
    %   - remove self-loops
    % return cut represented by final 2 vertices.
    
    while n > 2
        % 1) Pick a remaining edge (u,v) uniformly at random
        u = randi(n);
        v = randi(n);
        while v == u
            v = randi(n);
        end;
        
        % 2) Merge/contract u and v into a single vertex.
        %    Create new V, L, and n values to reflect the merged vertices.
        %    (MATLAB is call-by-value for arrays :-( )
        [mergedV, mergedL, mergedN] = MergeVertices(V, L, n, u, v);
        
        % 3) Remove self-loops in V.
        trimmedV = RemoveSelfLoops(mergedV, mergedN);
        
        V = trimmedV;
        L = mergedL;
        n = mergedN;
        
        %DisplayGraph(V, L, n);
    end;
    
    % Display the number of crossing edges for the min cut we found
    crossingNum = length(V{1});
    if mod(count, 30) == 0
        fprintf(' %i\n', crossingNum);
    else
        fprintf(' %i', crossingNum);
    end;
    bestCrossingNumber = min(crossingNum, bestCrossingNumber);
    %DisplayGraph(V, L, n);
    if length(V{1}) ~= length(V{2})
        fprintf('  --ERROR: length(V{1})=%i  ~=  length(V{2})=%i\n', length(V{1}), length(V{2}));
        error('Inconsistent vertex counts');
    end;
end;
fprintf('\nFewest number of crossing edges = %i\n', bestCrossingNumber);




