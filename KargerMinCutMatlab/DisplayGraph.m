function DisplayGraph(V, L, n)
    % Display the graph's vertices and labels.
    fprintf('      -----------------\n');
    fprintf('      Graph: %i nodes\n', n);
    fprintf('        Labels:\n', n);
    for i = 1:n
        fprintf('        %i: ', i);
        numLabels = length(L{i});
        for j = 1:numLabels
            fprintf('%i ', L{i}(j));
        end;
        fprintf('\n');
    end;
    fprintf('        Adjacent Vertices:\n', n);
    for i = 1:n
        fprintf('        %i: ', i);
        numAdjVertices = length(V{i});
        for j = 1:numAdjVertices
            fprintf('%i ', V{i}(j));            
        end;
        fprintf('\n');
    end;
    fprintf('      -----------------\n');
end

