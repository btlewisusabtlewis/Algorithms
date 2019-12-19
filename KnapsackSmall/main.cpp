//
//  main.cpp
//  KnapsackSmall
//
//  Created by Brian Lewis on 1/28/17.
//  Copyright © 2017 Brian Lewis. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>    // for std::min
#include "assert.h"

typedef std::vector< std::vector<int> >   Matrix;


// The number and list of edges read from the input file.
int n = 0;
std::vector<int> value;
std::vector<int> weight;

// Variables for the dynamic programming Knapsack algorithm
int W = 0;      // the knapsack's capacity

bool dbg = false;


//----------------------------------------------------------------------------------------------------------------

// Read items from the input file
void readInputFile(const char *fileName) {
    int numItemsRead = 0;
    int val, wt, i;
    
    // Read the edge descriptions from the input file.
    std::cout << "Input file: " << fileName << "\n";
    std::ifstream inputFile;
    inputFile.open(fileName);
    if (inputFile.is_open()) {
        inputFile >> W >> n;
        if (inputFile.eof()) {
            return;
        }
        if (dbg) std::cout << "Knapsack size " << W << ", number of items " << n << "\n";
        
        // Make room for (n+1) entries in the value and weight vectors, and initialize them to all "unset" entries.
        // The extra entry is because the entries are labelled 1..n (not starting at 0).
        value.assign((n+1), 0);
        weight.assign((n+1), 0);
        
        for (i = 1;  i <= n;  i++) {
            inputFile >> val >> wt;
            if (inputFile.eof()) {
                break;
            }
            numItemsRead++;
            if (dbg && (numItemsRead < 20)) std::cout << "  Item " << i << ": value " << val << ", weight " << wt << "\n";
            value[i]  = val;
            weight[i] = wt;
        }
    }
    inputFile.close();
    std::cout << "Number of items actually read = " << numItemsRead << "\n";
    if (numItemsRead != n) {
        std::cout << "ERROR: number of items read " << numItemsRead << " doesn't match expected item count " << n << "!\n";
    }
    std::cout << "-------------------------\n";
}


//----------------------------------------------------------------------------------------------------------------

int solveKnapsackProblem() {
    int x, i;
    
    // Allocate the 2D array A used to record solutions to (values of) smaller subproblems.
    // A[n+1][W+1] is a row vector of (n+1) elements, where each element is a column vector of (W+1) ints.
    Matrix A(n+1, std::vector<int>(W+1));  // can be very large for many items n and especially, a large capacity W
    
    // Initialize A
    for (x = 0;  x <= W;  x++) {
        A[0][x] = 0;
    }
    
    // Now use the knapsack recurrence to systematically solve all the increasingly large subproblems.
    for (i = 1;  i <= n;  i++) {
        if (dbg) std::cout << "  Item " << i << ": ";
        for (x = 0;  x <= W;  x++) {
            if (x >= weight[i]) {
                A[i][x] = std::max(A[i-1][x], (A[i-1][x-weight[i]] + value[i]));
            } else {
                A[i][x] = A[i-1][x];
            }
            if (dbg) std::cout << "(" << x << ")" << A[i][x] << ", ";
        }
        if (dbg) std::cout << "\n";
    }
    return A[n][W];
}


//----------------------------------------------------------------------------------------------------------------

int main(int argc, const char * argv[]) {
    if (argc < 2) {
        std::cout << "Usage: knapsacksmall <input file>\n";
        return 1;
    }
    std::cout << "The dynamic programming Knapsack algorithm\n";
    
    // Read the edge descriptions from the input file.
    readInputFile(argv[1]);
    
    // Run the dynamic programming algorithm from lecture to solve the Knapsack problem
    int optimalValue = solveKnapsackProblem();
    
    // Report the overall cost of T
    std::cout << "\nOptimal value of the knapsack = " << optimalValue << "\n";
    std::cout << "\nProgram finished\n";
    return 0;
}
