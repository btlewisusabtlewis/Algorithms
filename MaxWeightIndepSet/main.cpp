//
//  main.cpp
//  MaxWeightIndepSet
//
//  Created by Brian Lewis on 1/21/17.
//  Copyright © 2017 Brian Lewis. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <climits>
#include "assert.h"

//--------------------------------------------------------------------------------------------------------

#define NUM_SPECIFIED_VERTICES 8

// The vertex weights and their count.
unsigned n = 0;
std::vector<unsigned> w;

// Values of the max-weight independent sets for the initial subsequences (subgraphs) of the path's vertices
std::vector<long long unsigned> A;

// Set containing the max-weight independent set itself (vertex indices), not the MWIS value
// for the whole path of n vertices
std::set<unsigned> S;

// Vertex indices to look up in the computed MWIS
unsigned specifiedVertices[NUM_SPECIFIED_VERTICES] = {1, 2, 3, 4, 17, 117, 517, 997};

bool dbg = false;

//--------------------------------------------------------------------------------------------------------

// Read vertices (i.e., their weights) from the input file and add to the "w" vector
int readInputFile(const char *fileName) {
    // Read the point descriptions from the input file.
    std::cout << "Input file: " << fileName << "\n";
    std::ifstream inputFile;
    inputFile.open(fileName);
    if (inputFile.is_open()) {
        inputFile >> n;
        if (inputFile.eof()) {
            return -1;
        }
        std::cout << "Number of vertices = " << n << "\n";
        
        w.resize(n+1);      // the 1 is for the unused zero-th element
        w[0] = 0;           // unused: the vertices are numbered [1..n]
        
        // Read the weights of each vertex
        unsigned weight;
        if (dbg) std::cout << "Input vertex weights:\n";
        for (int i = 1;  i <= n;  i++) {
            inputFile >> weight;
            if (inputFile.eof()) {
                return -1;  // unexpected EOF
            }
            if (dbg && (i < 20)) std::cout << "  (" << i << ") " << weight << "\n";
            w[i] = weight;
        }
        if (dbg) std::cout << "\n";
        inputFile.close();
    }
    
    std::cout << "-------------------------\n";
    return 0;
}

//--------------------------------------------------------------------------------------------------------

void printBitset(unsigned bitset, unsigned numBits) {
    unsigned i;
    std::vector<unsigned> bit;
    bit.resize(numBits);
    for (i = 0;  i < numBits;  i++) {
        bit[i] = (bitset & 0x1);
        bitset = (bitset >> 1);
    }
    for (i = numBits;  i >= 1;  i--) {
        std::cout << (bit[i-1]? "1" : "0");
    }
}

void printVertexSetMembers() {
    bool prev = false;
    std::cout << "\nMax-weight IS vertex indices: [";
    for (unsigned i=1;  i <= n;  i++) {
        if (S.count(i) > 0) {
            if (prev) std::cout << ", ";
            std::cout << i;
            prev = true;
        }
    }
    std::cout << "]\n";
}

//--------------------------------------------------------------------------------------------------------

void computeMWISValuesOfSubpaths() {
    if (dbg) std::cout << "\nComputing max-weight IS values of the subpaths\n";
    A.resize(n+1);     // the A elements go from [0..n]
    // Initialization
    A[0] = 0;
    A[1] = w[1];
    // Main loop of the bottom-up iterative MWIS dynamic programming algorithm
    for (unsigned i=2;  i <= n;  i++) {
#if 0
        A[i] = std::max(A[i-1], (A[i-2] + w[i]));
#else
        if (A[i-1] >= (A[i-2] + w[i])) {  // case 1 wins
            A[i] = A[i-1];
            if (dbg) std::cout << "  Case 1: A[" << i << "] = " << A[i] << "\n";
        } else {                          // case 2 wins
            A[i] = (A[i-2] + w[i]);
            if (dbg) std::cout << "  Case 2: A[" << i << "] = " << A[i] << "\n";
        }
#endif
    }
}

void reconstructOptimalMWISSoln() {
    if (dbg) std::cout << "\nReconstructing the vertex indices of the max-weight IS\n";
    // A vertex vi belongs to a max-weight IS of Gi iff
    //    ((max-weight IS of Gi-2) + wi) >= (max-weight IS of Gi-1)
    S.clear();
    unsigned i = n;
    while (i >= 2) {
        if (A[i-1] >= (A[i-2] + w[i])) {  // case 1 wins
            if (dbg) std::cout << "  Case 1: skipping vertex index " << i << ", decr i by 1\n";
            i--;
        } else {                          // case 2 wins
            S.insert(i);
            if (dbg) std::cout << "  Case 2: adding vertex index " << i << " to S, decr i by 2\n";
            i = i-2;
        }
    }
    if (i == 1) {
        // case 2 wins
        S.insert(i);
        if (dbg) std::cout << "  Case 2: adding vertex index " << i << " to S, special test\n";
    }
}

// Look up the vertex indices from "specifiedVertices" in the MWIS set S and return an unspecified
// with a bit set to 1 if the corresponding specified vertex is found in S.
unsigned lookupSpecifiedVertices() {
    unsigned bitset = 0;
    for (unsigned i = 0;  i < NUM_SPECIFIED_VERTICES;  i++) {
        bitset = (bitset << 1);
        if (S.count(specifiedVertices[i]) > 0) {
            bitset = (bitset | 0x1);
            if (dbg) {
                std::cout << "  specified vertex " << i << " = " << specifiedVertices[i]
                          << " is in S, bitset now " << bitset << "\n";
                std::cout << "  bitset = ";  printBitset(bitset, NUM_SPECIFIED_VERTICES);  std::cout << "\n";
            }
        }
    }
    return bitset;
}

unsigned mwis() {
    computeMWISValuesOfSubpaths();
    std::cout << "Max sum: " << A[n] << "\n";
    reconstructOptimalMWISSoln();
    printVertexSetMembers();
    return lookupSpecifiedVertices();
}

//--------------------------------------------------------------------------------------------------------

int main(int argc, const char * argv[]) {
    std::cout << "Max-weight independent set algorithm\n";
    if (argc < 2) {
        std::cout << "Usage: mwis <input file>\n";
        return -1;
    }
    
    // Read the weights for the path vertices from the input file
    readInputFile(argv[1]);

    // Compute a max-weight independent set for the given vertices and look up the specified
    // vertices for membership in that set.
    unsigned bitset = mwis();
    std::cout << "\nFor the given vertices, the membership bitset for the specified vertices = " << bitset << "\n";
    std::cout << "The bitset = ";  printBitset(bitset, NUM_SPECIFIED_VERTICES);  std::cout << "\n";
    std::cout << "\nProgram finished\n";
    return 0;
}
