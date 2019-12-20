function [res, isNeg] = subNumbers(X, Y)
    % subNumbers(X, Y) Compute difference of two n-digit numbers X and Y 
    % stored as arrays of n digits, returning a new array "res" with the
    % absolute value of the result |X-Y| and "isNeg" set 1 if the result is
    % negative. We assume X and Y are positive, and depend on the caller to
    % do appropriate add or subNumbers calls.
    
    % Get the largest length (number of digits) of the input arrays
    n = max(length(X), length(Y));
    res = [zeros(n, 1)];
    isNeg = 0;
    
    % Flip the input arrays L<->R to simplify indexing from the lowest digit
    % to the highest one. So, [7;8;9] -> [9; 8; 7].
    % If needed, append zeros so both arrays are length n. So [4;5] -> [5;4;0]
    a = flipud(X);
    while length(a) < n
        a = [a; 0];
    end;
    b = flipud(Y);
    while length(b) < n
        b = [b; 0];
    end;
    
    % Subtract a-b from lowest-order (index 1) digit to highest-order (index n)
    for i = 1:n
        if b(i) <= a(i)
            % Easy case: no borrowing needed
            res(i) = a(i) - b(i);
        else 
            % Ugh: Borrow from left digit a(i+1) if possible. Otherwise, 
            % borrow from some further left digit of a. If that's 
            % not possible (we run out of digits of a to borrow from), 
            % the result a-b is negative (special case).
            j = i+1;
            while (j <= n) && (a(j) == 0)
                j = j+1;
            end;
            if j > n
                % Never found a digit to borrow from, result is negative
                isNeg = 1;
                break;
            end;
            % We borrow from the digit at j, and adjust b's digits between
            % j and i, all of which must be zero.
            a(j) = a(j)-1;
            k = j-1;
            while k > i 
                assert(a(k) == 0);
                a(k) = 9;
                k = k-1;
            end;
            % Now increment digit a(i) with the borrowed value 
            a(i) = a(i) + 10;
            % Subtract b(i) from the adjusted, larger a(i)
            res(i) = a(i) - b(i);
        end;
    end;
    
    if isNeg == 1
        % The result is negative. Set res to b-a (absolute value of a-b).
        [res, wontBeIsNeg] = subNumbers(flipud(b), flipud(a));
        assert(wontBeIsNeg == 0);
    else
        res = flipud(res); 
    end;    
end

