//
//  main.cpp
//  AllPairShortestPath - solution using the Floyd-Warshall algorithm
//
//  Created by Brian Lewis on 2/6/17.
//  Copyright © 2017 Brian Lewis. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>    // for std::min
#include <limits>       // for INT_MAX
#include "assert.h"
#include <iostream>

typedef std::vector< std::vector<int> >   Matrix;


// The number of vertices and edges.
int n = 0;
int m = 0;

// Variables for the Floyd-Warshall all-pair shortest-path algorithm
Matrix  one;
Matrix  two;
Matrix *Curr = NULL;          // refers to either one or two, alternating
Matrix *Prev = NULL;          // refers to the other matrix
bool    hasNegCycle = false;  // set true if the input graph has a negative cycle

bool dbg = false;


//----------------------------------------------------------------------------------------------------------------

void resizeMatrix(Matrix &vec, const unsigned int rows, const unsigned int columns) {
    vec.resize(rows);
    for (Matrix::iterator it = vec.begin();  it != vec.end();  ++it) {
        it->resize(columns);
    }
}


// Read items from the input file
void readInputAndInitArrayA(const char *fileName) {
    int numEdgesRead = 0;
    int i, j, e, length;
    
    // We can't allocate a 3D array A to record solutions to values of smaller subproblems since it'd be too big.
    // So, we allocate two 2D arrays one and two to each hold a "slice" of A. These are referenced through pointers,
    // and contain the current and previous 2D array slices of A being read or written ("Curr" and "Prev").
    std::cout << "Input file: " << fileName << "\n";
    std::ifstream inputFile;
    inputFile.open(fileName);
    if (inputFile.is_open()) {
        inputFile >> n >> m;
        if (inputFile.eof()) {
            return;
        }
        if (dbg) std::cout << "Number of vertices " << n << ", number of edges " << m << "\n";
        
        resizeMatrix(one, (n+1), (n+1));    // the array dimensions are indexed from 1 to n
        resizeMatrix(two, (n+1), (n+1));
        
        // Initialize the matrix pointers
        Curr = &one;
        Prev = &two;
        
        // Base cases: Initialize array A's k=0 slice, and set Prev to refer to it.
        for (i = 1;  i <= n;  i++) {
            // Set A[i,j,0] = +infinity for all A elements for now
            for (j = 1;  j <= n;  j++) {
                (*Prev)[i][j] = INT_MAX;
            }
            // Set A[i,j,0] = 0 for i=j
            (*Prev)[i][i] = 0;
        }
        for (e = 1;  e <= m;  e++) {
            // Now set A[i,j,0] = length of edge (i,j) if it exists in the input
            // NB: Not the same i and j used in the loop above
            inputFile >> i >> j >> length;
            if (inputFile.eof()) {
                break;
            }
            numEdgesRead++;
            if (dbg && (numEdgesRead < 20)) std::cout << "  Edge " << e << ": (" << i << "," << j << ") length " << length << "\n";
            (*Prev)[i][j] = length;
        }
    }
    inputFile.close();
    std::cout << "Number of edges actually read = " << numEdgesRead << "\n";
    if (numEdgesRead != m) {
        std::cout << "ERROR: number of edges read " << numEdgesRead << " doesn't match expected edge count " << m << "!\n";
    }
    
    if (dbg) {
        std::cout << "\nInitialized A[*,*,0]:\n";
        for (i = 1;  i <= n;  i++) {
            std::cout << "  Row " << i << ":\n";
            for (j = 1;  j <= n;  j++) {
                std::cout << "    (" << i << "," << j << "): " << (*Prev)[i][j] << "\n";
            }
        }
        std::cout << "\n";
    }
    std::cout << "-------------------------\n";
}


//----------------------------------------------------------------------------------------------------------------

int solveAllPairShortestPath() {
    int i, j, k;
    
    // Use the Floyd-Warshall recurrence to systematically solve all the increasingly large APSP subproblems.
    for (k = 1;  k <= n;  k++) {
        if (k > 1) {
            Matrix *tmp = Curr;
            Curr = Prev;
            Prev = tmp;
        }
        
        if (dbg) std::cout << "Slice k = " << k << ":\n";
        for (i = 1;  i <= n;  i++) {
            if (dbg) std::cout << "  Row " << i << ":\n";
            for (j = 1;  j <= n;  j++) {
                int case1Len = (*Prev)[i][j];  // Case 1: vertex k is not internal to shortest i-j path
                int l2 = (*Prev)[i][k];
                int l3 = (*Prev)[k][j];
                int case2Len =                 // Case 2: k is internal to P
                    ((l2 == INT_MAX) || (l3 == INT_MAX))? INT_MAX : (l2 + l3);
                (*Curr)[i][j] = std::min(case1Len, case2Len);
                if (dbg) std::cout << "    (" << i << "," << j << "): " << (*Curr)[i][j] << "\n";
            }
        }
    }
    
    // Check for a negative cycle in the input graph
    hasNegCycle = false;
    for (i = 1;  i <= n;  i++) {
        if ((*Curr)[i][i] < 0) {
            hasNegCycle = true;
        }
    }
    
    // If no negative cycle, compute the shortest shortest path dist(u,v) where u,v are in V and u != v.
    int minShortestPath = INT_MAX;
    if (!hasNegCycle) {
        for (i = 1;  i <= n;  i++) {
            for (j = 1;  j <= n;  j++) {
                minShortestPath = std::min((*Curr)[i][j], minShortestPath);
            }
        }
    }
    
    return minShortestPath;
}


//----------------------------------------------------------------------------------------------------------------

int main(int argc, const char * argv[]) {
    if (argc < 2) {
        std::cout << "Usage: allpairsshortestpath <input file>\n";
        return 1;
    }
    std::cout << "The Floyd-Warshall all-pair shortest-path algorithm\n";
    
    // Read the edge descriptions from the input file and initialize the two 2-D arrays that represent slices of the 3-D array A.
    readInputAndInitArrayA(argv[1]);
    
    // Run the Floyd-Warshall dynamic programming algorithm to solve the APSP problem, and return the shortest shortest path
    int minShortestPath = solveAllPairShortestPath();
    
    // Report the shortest shortest path, if one exists
    if (hasNegCycle) {
        std::cout << "\nInput graph has a negative cycle!\n";
    } else {
        std::cout << "\nShortest shortest path in the input graph = " << minShortestPath << "\n";
    }
    std::cout << "\nProgram finished\n";
    return 0;
}
