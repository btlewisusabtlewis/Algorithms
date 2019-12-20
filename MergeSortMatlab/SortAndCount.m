function [SortedA, inversionCount] = SortAndCount(A, n)
    %% SortAndCount(A, n) Sort the input array A and count the number of its 
    %  inversions, then return both the sorted array and the inversion count.

    if n <= 1
        SortedA = A;
        inversionCount = 0;
        return;
    else
        firstHalf  = floor(n/2);
        secondHalf = (n - firstHalf);
        [B, X] = SortAndCount(A(1:firstHalf),     firstHalf);
        [C, Y] = SortAndCount(A((firstHalf+1):n), secondHalf);
        [SortedA, Z] = MergeAndCountSplitInv(B, C, firstHalf, secondHalf);
        inversionCount = (X + Y + Z);
    end;
    