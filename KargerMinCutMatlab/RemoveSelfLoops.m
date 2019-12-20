function [trimmedV] = RemoveSelfLoops(V, n)
    % Remove any self-loops in the vertex array V.
    % Probably don't need to consider all elements of V, just u and v in
    % MinCut.m.
    trimmedV = cell(length(V), 1);
    for i = 1:n
        numAdjVertices = length(V{i});
        trimmedV{i} = zeros(numAdjVertices, 1);
        k = 0;
        for j = 1:numAdjVertices
            if V{i}(j) == i   % a self loop to be removed
                %fprintf('      Removing self-ref in vertex %i\n', i);
            else                 
                k = k+1;
                trimmedV{i}(k) = V{i}(j);
            end;
        end;
        trimmedV{i} = trimmedV{i}(1:k);  %sort(trimmedV{i}(1:k));
    end; 
    trimmedV = trimmedV;
end

