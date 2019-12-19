//
//  main.cpp
//  KnapsackBig
//
//  Created by Brian Lewis on 1/28/17.
//  Copyright © 2017 Brian Lewis. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>    // for std::min
#include "assert.h"
#include <iostream>

// The number and list of edges read from the input file.
int n = 0;
std::vector<int> value;
std::vector<int> weight;

// Variables for the dynamic programming Knapsack algorithm
int W = 0;      // the knapsack's capacity
std::vector<int> one;
std::vector<int> two;
std::vector<int> *Curr = NULL;   // refers to either one or two, alternating
std::vector<int> *Prev = NULL;   // refers to the other vector

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
        std::cout << "Knapsack size " << W << ", number of items " << n << "\n";
        
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
    
    // We can't allocate a 2D array A to record solutions to values of smaller subproblems since it'd be too big.
    // So, we allocate two vectors to each hold a column of A. These are referenced through pointers, and contain
    // the current column of A being computed ("Curr"), and the previous column ("Prev") used for that computation.
    // Current[x] holds the values of column A[*,x].
    one.resize(W+1, 0);
    two.resize(W+1, 0);
    
    // Initialize the vector pointers
    Curr = &one;
    Prev = &two;
    
    // Initialize Prev to the 0-th column of A, which is all zero
    for (x = 0;  x <= W;  x++) {
        (*Prev)[x] = 0;
    }
    
    // Now use the knapsack recurrence to systematically solve all the increasingly large subproblems.
    for (i = 1;  i <= n;  i++) {
        int weightI = weight[i];
        if (i > 1) {
            std::vector<int> *tmp = Curr;
            Curr = Prev;
            Prev = tmp;
        }
        //if (dbg) std::cout << "  Item " << i << ": ";
        for (x = 0;  x <= W;  x++) {
            if (x >= weightI) {
                (*Curr)[x] = std::max((*Prev)[x], ((*Prev)[x-weightI] + value[i]));
            } else {
                (*Curr)[x] = (*Prev)[x];
            }
            //if (dbg) std::cout << "(" << x << ")" << (*Curr)[x] << ", ";
        }
        //if (dbg) std::cout << "\n";
    }
    return (*Curr)[W];
}


//----------------------------------------------------------------------------------------------------------------

int main(int argc, const char * argv[]) {
    if (argc < 2) {
        std::cout << "Usage: knapsackbig <input file>\n";
        return 1;
    }
    std::cout << "The big dynamic programming Knapsack algorithm\n";
    
    // Read the edge descriptions from the input file.
    readInputFile(argv[1]);
    
    // Run the dynamic programming algorithm from lecture to solve the Knapsack problem
    int optimalValue = solveKnapsackProblem();
    
    // Report the overall cost of T
    std::cout << "\nOptimal value of the knapsack = " << optimalValue << "\n";
    std::cout << "\nProgram finished\n";
    return 0;
}
