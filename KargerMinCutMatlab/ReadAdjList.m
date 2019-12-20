function [V, L, n] = ReadAdjList(fileName)
    % Read the input file containing an adjacency list for a graph.
    % For each input vertex v, return 
    %    1) a set V indexed by vertex number that contains an array of v's 
    %       adjacent vertices and
    %    2) a set L indexed by vertex number that contains an array of v's 
    %       labels (initially holding just v).
    
    % First, read each line into an array and get the number of vertices n
    n = 0;
    fileID = fopen(fileName, 'r');
    tline = fgetl(fileID);
    while ischar(tline)
        n = n+1;
        tline = fgetl(fileID);
    end;
    fclose(fileID);
    
    lines = strings(n, 1);
    
    i = 0;
    fileID = fopen(fileName, 'r');
    tline = fgetl(fileID);
    while ischar(tline)
        i = i+1;
        lines(i) = tline;
        tline = fgetl(fileID);
    end;    
    fclose(fileID);
    
    % Now scan each line for the vertex number 
    V = cell(n, 1);
    L = cell(n, 1);
    for i = 1:n
        % Read a vertex number then the edges incident upon it
        lineContents = sscanf(lines(i), '%d');
        vertexNum = lineContents(1);
        assert(vertexNum == i);
        V{vertexNum} = lineContents(2:end);
        L{vertexNum} = vertexNum;  % only one label for now
    end;
end

