function sum = addNumbers(X, Y)
    %% addNumbers(X, Y) Compute sum of two n-digit numbers X and Y stored as
    % arrays of n digits, returning a new array sum.

    % Get the length (number of digits) of the input arrays
    n = max(length(X), length(Y));
    % The result might have one more digit
    res   = [zeros(n+1, 1)];
    % We need a carry input value for all digits but the lowest-order one
    carry = [zeros(n+1, 1)];
    
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
    
    % Add from lowest-order (index 1) digit to highest-order (index n)
    for i = 1:n
        res(i) = carry(i) + a(i) + b(i);
        if res(i) >= 10
            carry(i+1) = 1;
            res(i) = (res(i) - 10);
        end;
    end;
    % Check for carry out to the high-order (n+1)th digit
    if carry(n+1) ~= 0
        res(n+1) = 1;
        sum = flipud(res);
    else
        sum = flipud(res(1:n));
    end
end
