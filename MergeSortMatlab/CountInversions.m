%% Algorithms 1 Exercise 2: Count inversions using divide and conquer
%%                          by piggy-backing on merge sort.

% Initialization
clear ; close all; clc

% Read the input file
fileID = fopen('IntegerArray.txt', 'r');
A = fscanf(fileID, '%i');
fclose(fileID);
n = length(A);

[SortedA, numInversions] = SortAndCount(A, n);

fprintf('Computed the number of inversions = %i\n', numInversions);
pause; %----------------------
