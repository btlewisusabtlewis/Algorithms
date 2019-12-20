function [p, pIndex] = ChoosePivot(A, l, r, pivotKind)
    % Choose pivot value p (at element pIndex) for array A[l..r].
    
    numElements = (r-l+1);
    assert(numElements >= 2);
    % Choose pivot value p (at element pIndex)
    switch pivotKind
        case 1  % Use the first element of the array
            pIndex = l;
        case 2  % Use the last element of the array
            pIndex = r;
        case 3  % Use the median-of-three element of the array
            if mod(numElements, 2) == 0  % is even
                midIndex = l + (numElements/2) - 1;
            else
                midIndex = l + ((numElements+1)/2) - 1;
            end;
            theThree = [A(l), A(midIndex), A(r)];
            med = median(theThree);
            pIndex = find(A==med);
            assert(length(pIndex) == 1);
            %if numElements <= 10
            %    fprintf('  For array [');
            %    fprintf('%i ', A(l:r));
            %    fprintf('], l=%i, r=%i, med=%i, pIndex=%i, A(pIndex)=%i\n', ...
            %            l, r, med, pIndex, A(pIndex));
            %end;
        otherwise
            fprintf('ERROR: bad pivotKind value %i\n', pivotKind);
            assert(0);
    end;
    p = A(pIndex);
end

