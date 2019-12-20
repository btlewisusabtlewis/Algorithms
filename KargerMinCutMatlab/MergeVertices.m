function [mergedV, mergedL, mergedN] = MergeVertices(V, L, n, u, v)
    % Merge/contract vertices u and v into a single vertex in the 
    % adjacency list with vertex set V, label set L, and vertex count n.
    
    assert(n == length(V));
    assert(n == length(L));
    assert(u ~= v);
    
    % First, make sure that u < v (swap if necessary).
    % This means we always merge v into u, and then discard v afterwards.
    if u > v
        temp = u;  u = v;  v = temp;
    end;
    %fprintf('  Merging vertex %i into %i\n', v, u);
    
    % Merge vertex v's edges and labels into vertex u in the original array
    V{u} = [V{u}; V{v}];
    L{u} = [L{u}; L{v}];
    
    % Clear out v to aid debugging
    V{v} = zeros(1, 1);
    L{v} = zeros(1, 1);
    
    % Update any reference to vertex v to now refer to u. 
    for i = 1:n
        numAdjVertices = length(V{i});
        for j = 1:numAdjVertices
            if V{i}(j) == v
                %fprintf('    Correcting vertex %i ref from %i to %i\n', i, v, u);
                V{i}(j) = u;
            end;
        end;
    end;
    
    if v < n
        % Move V's and L's elements after v (those in [(v+1):n]) up by one
        for i = v:(n-1)  % 5:6
            %fprintf('    Moving vertex %i up by one\n', (i+1));
            numAdjVertices = length(V{i+1});
            V{i} = zeros(numAdjVertices, 1);
            for j = 1:numAdjVertices
                V{i}(j) = V{i+1}(j);
            end;
            
            numLabels = length(L{i+1});
            L{i} = zeros(numLabels, 1);
            for j = 1:numLabels
                L{i}(j) = L{i+1}(j);
            end;
        end;
        
        % Clear out V{n}/L{n} to aid debugging
        V{n} = zeros(1, 1);
        L{n} = zeros(1, 1);
    
        % Correct any references to vertex (i+1) to i.
        for i = 1:(n-1)
            numAdjVertices = length(V{i});
            for j = 1:numAdjVertices
                if V{i}(j) > v
                    %fprintf('    Correcting vertex %i ref from %i to %i\n', i, V{i}(j), (V{i}(j) - 1));
                    V{i}(j) = (V{i}(j) - 1);
                end;
            end;
        end;
    end;
    
    % Copy over the non-v elements to the new merged edge and label arrays.
    mergedN = (n - 1);
    mergedV = V(1:mergedN);
    mergedL = L(1:mergedN);
end

