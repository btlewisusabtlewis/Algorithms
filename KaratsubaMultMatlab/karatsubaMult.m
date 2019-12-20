function xy = karatsubaMult(x, y)
    %% karatsubaMult(x, y) Compute product of two n-digit numbers 
    % represented as arrays of n digits, returning a new digit array xy.
    
    n = max(length(x), length(y));
    assert(n > 0);
    if n == 1
        % Base case: just multiply the digits and return
        xy = [zeros(2, 1)];
        temp = (x(1) * y(1));
        if temp >= 10
            xy(1) = floor(temp/10);
            xy(2) = mod(temp, 10);
        else
            xy(2) = temp;
        end;
        return;
    end;
    
    n = 2^nextpow2(n);
    assert(mod(n,2) == 0);
    
    % Enlarge x or y to n digits if necessary
    if length(x) < n
        temp = [zeros(n, 1)];
        temp((n-length(x)+1) : end) = x(1:end);
        x = temp;
    end;
    if length(y) < n
        temp = [zeros(n, 1)];
        temp((n-length(y)+1) : end) = y(1:end);
        y = temp;
    end
    
    % Divide the input number arrays x and y into [a . b] and [c . d] subarrays.
    % I.e., x = 10^(n/2)*a + b and y = 10^(n/2)*c + d.
    half = n/2;         % integer since n is even
    a = x(1:half);
    b = x((half+1):n);
    c = y(1:half);
    d = y((half+1):n);
    
    % Recursively compute a*c
    ac = karatsubaMult(a, c);                 % (1)
    
    % Recursively compute b*d
    bd = karatsubaMult(b, d);                 % (2)
    
    % Recursively compute (a+b)*(c+d)
    aPlusB = addNumbers(a, b);
    cPlusD = addNumbers(c, d);
    multSums = karatsubaMult(aPlusB, cPlusD);  % (3)
    
    % Gauss' trick: form (a*d + b*c) as gauss = ((3) - (1)) - (2).
    [threeSubOne, is31Neg] = subNumbers(multSums, ac);
    assert(is31Neg == 0);
    isGaussNeg = 0;
    if is31Neg
        gauss = addNumbers(threeSubOne, bd);
        isGaussNeg = 1;
    else
        [gauss, isGaussNeg] = subNumbers(threeSubOne, bd);
        assert(isGaussNeg == 0);
    end;
    
    % Now form result xy using shifts and addition
    res = [zeros(n*2, 1)];
    res(1:n)         = ac;    % high-order n digits of (2n) digit result
    res((n+1):(2*n)) = bd;    % low-order n digits 
    % Add the partial res and gauss shifted half digits left to form x*y
    res2 = [zeros(n*2, 1)];
    lg = length(gauss);
    lenDiff = (2*n - lg);
    temp = [zeros(n*2, 1)];
    temp((lenDiff+1) : (lenDiff+lg)) = gauss(1:end);
    res2(1:(end-half)) = temp((half+1):end); 
    xy = addNumbers(res, res2);
end
    
