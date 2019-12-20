function [resultArr, numComparisons] = QuickSortAndCount(A, l, r, pivotKind)
    % QuickSort the array A of length n and count the number of comparisons
    % using different pivot choices. Return the sorted subarray and the
    % comparison count.
    % NB: Matlab passes arrays by value, so the array "A" here is a copy!
    
    numElements = (r-l+1);
    if numElements <= 1
        resultArr = A(l:r);
        numComparisons = 0;
    else
        % Choose pivot value p (at element pIndex)
        [p, pIndex] = ChoosePivot(A, l, r, pivotKind);

        % Partition A around p
        %   1) swap A(l) and A(pIndex)
        CheckIndex(l, r, pIndex);
        t = A(l);  A(l) = A(pIndex);  A(pIndex) = t;
        %   2) scan through A doing the partition
        i = (l+1);       % i points just after <p part = at first of >p part
        for j = (l+1):r  % j points just after >p part = after all partioning
            if A(j) < p
                % swap A(i) and A(j)
                CheckIndex(l, r, i);  CheckIndex(l, r, j);
                t = A(i);  A(i) = A(j);  A(j) = t;
                i = i+1;
            end;
        end;
        %   3) swap pivot into its right place (swap A(l) and A(i-1)).
        CheckIndex(l, r, (i-1));
        t = A(l);  A(l) = A(i-1);  A(i-1) = t;
        
        % Recursively sort first part (elements < p, note: p is at (i-1))
        firstLen = ((i-2) - l + 1);
        if firstLen > 0
            [sortedFirst, numFirstComp] = QuickSortAndCount(A, l, (i-2), pivotKind);
        else
            numFirstComp = 0;
        end;
        
        % Recursively sort second part (elements > p)
        secondLen = (r - i + 1);
        if secondLen > 0
            [sortedSecond, numSecondComp] = QuickSortAndCount(A, i, r, pivotKind);
        else
            numSecondComp = 0;
        end;
        
        % Result subarray
        if firstLen > 0
            resultArr(l:(i-2)) = sortedFirst;
        end;
        resultArr(i-1) = p;
        if secondLen > 0
            resultArr(i:r) = sortedSecond;
        end;
        % But, we only want the subarray slice that we sorted
        resultArr = resultArr(l:r);
        
        % Total number of comparisons
        %   1) the pivot p was compared to each of A[l..r]'s other elements
        numCompHere = max((numElements-1), 0);  
        %   2) and the recursive calls each did their own comparisons
        numComparisons = (numCompHere + numFirstComp + numSecondComp);
        
        % Debugging
        assert(IsSorted(resultArr, 1, length(resultArr)) == 1);
    end;    
end

