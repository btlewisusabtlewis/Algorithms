function [D, Z] = MergeAndCountSplitInv(B, C, lengthB, lengthC)
    %% MergeAndCountSplitInv(B, C, lengthB, lengthC) Merge the two sorted
    % arrays B (of length lengthB) and C (of length lengthC) creating a 
    % new sorted array D. Return D and the count of its split inversions Z.
    
    n = (lengthB + lengthC);
    D = zeros(n, 1);
    Z = 0;  % number of split inversions in sorted result D
    
    i = 1;  % index into B
    j = 1;  % index into C
    for k = 1:n
        if i > lengthB
            % We exhausted array B, so copy C(j) to D if it exists.
            if j <= lengthC
                % Copy C(j) to D.
                D(k) = C(j);
                j = j+1;
                % The number of split inversions involving C(j) are the
                % number of elements left in the first array B. That is,
                % C(j) is smaller (an inversion) than remaining elements of B.
                Z = Z + (lengthB - i + 1);
            end;
        elseif j > lengthC
            % We exhausted array C, so copy B(i) to D if it exists.
            if i <= lengthB
                D(k) = B(i);
                i = i+1;
            end;
        elseif B(i) < C(j)
            % Copy B(i) to D.
            D(k) = B(i);
            i = i+1;
        elseif C(j) < B(i)
            % Copy C(j) to D.
            D(k) = C(j);
            j = j+1;
            % The number of split inversions involving C(j) are the 
            % number of elements left in the first array B. That is,
            % C(j) is smaller (an inversion) than remaining elements of B.
            Z = Z + (lengthB - i + 1);
        else
            assert(0);
        end;
    end;
