%% Algorithms 1 Exercise 3: Count QuickSort comparisons using different 
%%                          pivot choices.

% Initialization
clear ; close all; clc

% Read the input file
%A = [3; 8; 2; 5];
%A = [3; 8; 2; 5; 1; 4; 7; 6];
%fileID = fopen('10.txt', 'r');
%fileID = fopen('100.txt', 'r');
%fileID = fopen('1000.txt', 'r');
fileID = fopen('QuickSort.txt', 'r');
A = fscanf(fileID, '%i');
fclose(fileID);
n = length(A);

% Always use the first element of the array as the pivot element
[sortedA1, numComparisons] = QuickSortAndCount(A, 1, n, 1);
fprintf('With pivot=first element, number of comparisons = %i\n', ...
        numComparisons);

% Always use the final element of the array as the pivot element
[sortedA2, numComparisons] = QuickSortAndCount(A, 1, n, 2);
fprintf('With pivot=last element, number of comparisons = %i\n', ...
        numComparisons);

% Always use the median-of-three (median of first, middle, and last 
% elements of the array) value as the pivot element.
[sortedA3, numComparisons] = QuickSortAndCount(A, 1, n, 3);
fprintf('With pivot=median-of-three, number of comparisons = %i\n', ...
        numComparisons);
