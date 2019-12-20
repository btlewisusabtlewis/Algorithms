function inOrder = IsSorted(A, l, r)
    % Return 1 if A[l..r] is sorted in increasing order, else 0.
    inOrder = 1;
    prev = A(l);
    for i = (l+1):r
        if A(i) < prev
            inOrder = 0;
            return;
        else
            prev = A(i);
        end;
    end;
end

