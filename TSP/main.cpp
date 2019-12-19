//
//  main.cpp
//  TSP - The traveling salesman problem, implemented using the O(n^2 2^n) dynamic programming algorithm from the class
//
//  Created by Brian Lewis on 2/15/17.
//  Copyright © 2017 Brian Lewis. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>    // for std::min
#include <cfloat>       // for FLT_MAX
#include <math.h>       // for sqrt, floor, and pow
#include <cstdint>      // for int32_t and uint32_t
#include <chrono>
#include "assert.h"
#include <iostream>


// Max number of cities for our tightly-encoded space-efficient TSP solution
// We encode A[S,j] in 32 bits: set S is size n bits while j is (32-n) bits.
// With 27 cities, we encode S in 27 bits, and 2<=j<=27 in 5 bits.
const int MAX_NUM_CITIES = 27;

typedef std::vector<float>                 Vector32;
typedef std::vector< std::vector<float> >  Matrix;


// The number of cities, and their x- and y-coordinates.
unsigned int n = 0;          // n must be <= 27 for our space-efficient solution
Vector32 xcoord;
Vector32 ycoord;

// Variables for the dynamic programming traveling salesman algorithm
Matrix   dist;              // dist[i][j] = distance between cities i and j, a float

Vector32  one;
Vector32  two;
Vector32 *Curr = NULL;       // refers to either one or two, alternating
Vector32 *Prev = NULL;       // refers to the other matrix

uint32_t maskNBits = 0;      // Used to index into the array A with its packed fields for the vertex set and destination vertex

const bool dbg = false;


//----------------------------------------------------------------------------------------------------------------

float readA(Vector32* &A, uint32_t S, uint32_t j) {
    uint32_t idx = ((j << n) | (S & maskNBits));
    return (*A)[idx];
}

void writeA(Vector32* &A, uint32_t S, uint32_t j, float value) {
    uint32_t idx = ((j << n) | (S & maskNBits));
    (*A)[idx] = value;
}

//----------------------------------------------------------------------------------------------------------------

void resizeMatrix(Matrix &vec, const unsigned int rows, const unsigned int columns) {
    vec.resize(rows);
    for (Matrix::iterator it = vec.begin();  it != vec.end();  ++it) {
        it->resize(columns);
    }
}

// Read cities from the input file
void readInput(const char *fileName) {
    int numCitiesRead = 0;
    int i, j;
    float x, y;
    
    std::cout << "Input file: " << fileName << "\n";
    std::ifstream inputFile;
    inputFile.open(fileName);
    if (inputFile.is_open()) {
        inputFile >> n;
        if (inputFile.eof()) {
            return;
        }
        std::cout << "Number of cities " << n << "\n";
        if (n > MAX_NUM_CITIES) {
            std::cout << "ERROR: Number of cities " << n << " is too large (more than " << MAX_NUM_CITIES << ")\n";
            return;
        }
        
        // Use size n+1 since our indices go from [1..n]
        xcoord.resize((n+1), 0);
        ycoord.resize((n+1), 0);

        for (i = 1;  i <= n;  i++) {
            inputFile >> x >> y;
            if (inputFile.eof()) {
                break;
            }
            numCitiesRead++;
            //if (dbg && (numCitiesRead < 20)) std::cout << "  City " << i << ": (" << x << "," << y << ")\n";
            xcoord[i] = x;
            ycoord[i] = y;
        }
    }
    inputFile.close();
    //std::cout << "Number of cities actually read = " << numCitiesRead << "\n";
    if (numCitiesRead != n) {
        std::cout << "ERROR: number of cities read " << numCitiesRead << " doesn't match expected city count " << n << "!\n";
    }
    
    // Compute inter-city distances
    resizeMatrix(dist, (n+1), (n+1));  // use size n+1 since the indices go from [1..n]
    for (i = 1;  i <= n;  i++) {
        //if (dbg) std::cout << "  Dist from " << i << " to:\n";
        for (j = 1;  j <= n;  j++) {
            dist[i][j] = sqrt( pow((xcoord[i] - xcoord[j]), 2) +
                               pow((ycoord[i] - ycoord[j]), 2) );
            //if (dbg) std::cout << "    " << j << ": " << dist[i][j] << "\n";
        }
    }
    
    // Compute needed constants and masks
    maskNBits = ((1 << n) - 1);  // only least n bits are set 1
    std::cout << "-------------------------\n";
}

void allocateAndInitArrayA() {
    int sliceSizes = ((1024 * 1024 * 1024) + 1);  // 32M + 1 since we index from 1..end
    
    auto t1 = std::chrono::high_resolution_clock::now();
    
    // We can't allocate a 3D array A to record solutions to values of all smaller subproblems since it'd be too big.
    // So, we allocate two 2D arrays one and two to each hold a "slice" of A. These are referenced through pointers,
    // and contain the current and previous 2D array slices of A being read or written ("Curr" and "Prev").
    one.resize(sliceSizes, FLT_MAX);
    two.resize(sliceSizes, FLT_MAX);
    
    // Initialize the matrix pointers
    Curr = &one;
    Prev = &two;
    
    // Base case: subproblem size m=1. Initialize A[S,1]=0 if S = {1}, otherwise +inf. Set Prev to refer to A with m=1.
    writeA(Prev, /*S*/ 1, /*j*/ 1, 0);
    
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Allocating A took " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << " ms\n";
}


//----------------------------------------------------------------------------------------------------------------

// Return true iff the bit of bitArray at bitNumber (in [1..32]) is 1.
bool isBitSet(uint32_t bitArray, uint32_t bitNumber) {
    return (bitArray & (1 << (bitNumber-1)));
}

// Return the argument bitArray with the bit at bitNumber (in [1..32]) cleared.
uint32_t clearBit(uint32_t bitArray, uint32_t bitNumber) {
    uint32_t bitSetOn = (1 << (bitNumber-1));
    return (bitArray & ~bitSetOn);
}

float solveTSPTour() {
    float minCostTour = FLT_MIN;
    const int twoToN = pow(2, n);   // 2^n
    uint32_t S;
    uint32_t j;
    auto totalt1 = std::chrono::high_resolution_clock::now();
    
    std::cout << "Starting the TSP solution\n";
    for (int m = 2;  m <= n;  m++) {
        std::cout << "\n  Processing subproblems of size m=" << m << "\n";
        auto t1 = std::chrono::high_resolution_clock::now();
        
        // Flip Curr and Prev if appropriate
        if (m > 2) {
            Vector32 *tmp = Curr;
            Curr = Prev;
            Prev = tmp;
        }
        
        // For each m-sized subset S of {1,2,...,n} containing vertex 1
        // Treat S as an array of n bits, so S = {1,4,5} = 0x19 = 0b11001 (with C++14 binary literals)
        for (S = 0x1;  S <= twoToN;  S += 2) {
            std::bitset<32> sBits(S);
            if (sBits.count() == m) {
                //if (dbg) std::cout << "    Subset S=" << S << " or " << sBits << "\n";
                // For each j != 1 in S...   (so j is in {2,...,n})
                for (j = 2;  j <= n;  j++) {
                    if (isBitSet(S, j)) {
                        // j is in S
                        uint32_t k;
                        //if (dbg) std::cout << "      With j=" << j << "\n";
                        // Find k in S, k != j, to minimize path length from 1 to j via next-to-last k: A[S-{j}, k] + dist[k][j]
                        uint32_t sMinusJ = clearBit(S, j);  // NB: the set sMinusJ has size m-1
                        //std::bitset<32> sMinusJBits(sMinusJ);
                        //if (dbg) std::cout << "        S-{j}=" << sMinusJ << " or " << sMinusJBits << "\n";
                        float minDistToJ = FLT_MAX;
                        uint32_t bestK = 0;
                        for (k = 1;  k <= n;  k++) {
                            if ((k != j) && isBitSet(S, k)) {
                                float pathViaSMinusJToK = readA(Prev, sMinusJ, k);
                                //if (dbg) std::cout << "          k=" << k << ": read A[S-{j},k] = A[" << sMinusJ << "," << k << "] = " << pathViaSMinusJToK << "\n";
                                //if (dbg) std::cout << "               dist[k][j] = dist[" << k << "," << j << "] = " << dist[k][j] << "\n";
                                float distViaSMinusJToJ = ((pathViaSMinusJToK == FLT_MAX)? FLT_MAX : (pathViaSMinusJToK + dist[k][j]));
                                //if (dbg) std::cout << "               A[S-{j},k] + dist[k,j] = " << distViaSMinusJToJ << "\n";
                                if (distViaSMinusJToJ < minDistToJ) {
                                    minDistToJ = distViaSMinusJToJ;
                                    bestK = k;
                                }
                            }
                        }
                        writeA(Curr, S, j, minDistToJ);
                        //if (dbg) std::cout << "        Min path from 1..j=" << j << " via best k=" << bestK << " is " << minDistToJ << "\n";
                        //if (dbg) std::cout << "        Wrote A[S,j] = A[" << S << "," << j << "] = " << minDistToJ << "\n";
                    }
                }
            }
        }
        
        auto t2 = std::chrono::high_resolution_clock::now();
        std::cout << "    Size m=" << m << " subproblems took " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << " ms\n";
    }
    
    // Return min{A[{1,2,...,n}, j] + dist[j][1]} for j=2,...,n
    // This is min cost from 1 to j of (cost of visiting all vertices once to dest j) + (cost of final hop back to start of tour)
    std::cout << "\n  Finding min path visiting all vertices\n";
    minCostTour = FLT_MAX;
    S = (twoToN - 1);   // all n bits set
    uint32_t bestDest = 0;
    for (j = 2;  j <= n;  j++) {
        if (isBitSet(S, j)) {
            float pathViaAllSToJ = readA(Curr, S, j);
            //if (dbg) std::cout << "    j=" << j << ": read A[allS,j] = A[" << S << "," << j << "] = " << pathViaAllSToJ << "\n";
            float fullTourDistViaJ = ((pathViaAllSToJ == FLT_MAX)? FLT_MAX : (pathViaAllSToJ + dist[j][1]));
            //if (dbg) std::cout << "      A[allS,j] + dist[j,1] = " << pathViaAllSToJ << " + " << dist[j][1] << " = " << fullTourDistViaJ << "\n";
            if (fullTourDistViaJ < minCostTour) {
                minCostTour = fullTourDistViaJ;
                bestDest = j;
            }
        }
    }
    //if (dbg) std::cout << "  Min TSP tour distance is " << minCostTour << " via " << bestDest << "\n";
    
    auto totalt2 = std::chrono::high_resolution_clock::now();
    std::cout << "Total TSP computation took " << std::chrono::duration_cast<std::chrono::milliseconds>(totalt2-totalt1).count() << " ms\n";
    return minCostTour;
}


//----------------------------------------------------------------------------------------------------------------

int main(int argc, const char * argv[]) {
    if (argc < 2) {
        std::cout << "Usage: tsp <input file>\n";
        return 1;
    }
    std::cout << "The traveling salesman problem, implemented using a dynamic programming algorithm\n";
    
    // Read the city coordinates from the input file.
    readInput(argv[1]);
    
    // Initialize the two 2-D arrays that represent slices of the 3-D array A.
    allocateAndInitArrayA();
    
    // Run the dynamic programming TSP algorithm and return the minimum cost of a traveling salesman tour.
    float minTourCost = solveTSPTour();
    
    // Report the minimum cost of a traveling salesman tour, rounded down to closest integer
    std::cout << "\nMinimum cost of a traveling salesman tour = " << minTourCost << " = " << floor(minTourCost) << "\n";
    std::cout << "\nProgram finished\n";
    return 0;
}
