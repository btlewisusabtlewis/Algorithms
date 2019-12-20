%% Algorithms 1 Exercise 1: Fast arbitrary digit multiplication using 
%%                          Karatsuba Multiplication

% Initialization
clear ; close all; clc

% Create arrays x and y holding the digits of the numbers to be multiplied
S1 = '3141592653589793238462643383279502884197169399375105820974944592';
S2 = '2718281828459045235360287471352662497757247093699959574966967627';
l1 = length(S1);
l2 = length(S2);
n = max(l1, l2);
n = ceil(n/2) * 2;  % round n up to the closest even number
assert(mod(n,2) == 0);
x = [zeros(n, 1)];
y = [zeros(n, 1)];
for i = 1:l1
    x(i) = str2num(S1(i));
end
for i = 1:l2
    y(i) = str2num(S2(i));
end

% Multiply x and y to create the product in xy
xy = karatsubaMult(x, y);
fprintf('Computed the final product:\n');
resStr = num2str(xy);
for i = 1:128
    fprintf('%c', resStr(i));
end;
fprintf('\n');
pause; %----------------------


