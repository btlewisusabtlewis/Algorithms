//
//  main.cpp
//  HeuristicTSP
//
//  Created by Brian Lewis on 2/19/17.
//  Copyright © 2017 Brian Lewis. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>    // for std::min
#include <cfloat>       // for FLT_MAX
#include <math.h>       // for sqrt, floor, and pow
#include <cstdint>      // for int32_t and uint32_t
#include <chrono>
#include "assert.h"
#include <iostream>

// Cities have (x,y) coordinates, an index, and a alreadyVisited boolean
struct City
{
    uint32_t index;       // the integer identifying the edge
    bool     alreadyVisited;
    double   x;
    double   y;
};

// The number of cities
unsigned int n = 0;

// Variables for the heuristic TSP algorithm
double tourCost = 0.0;              // sum of Euclidean distances for cities on the tour (tourCities)
std::list<City*> unvisitedCities;  // initially all the cities, which are moved if selected to tourCities
std::list<City*> tourCities;       // cities in the tour, in order of from the first to the last

const bool dbg = false;


//----------------------------------------------------------------------------------------------------------------

City *createCity(int i, double cityX, double cityY, bool visited) {
    City *c = new City;
    c->index = i;
    c->alreadyVisited = visited;
    c->x     = cityX;
    c->y     = cityY;
    return c;
}

// Read city descriptions from the input file
void readInputFile(const char *fileName) {
    int numCitiesRead = 0;
    int i, idx;
    double x, y;
    
    // Read the edge (point pair) descriptions from the input file.
    std::cout << "Input file: " << fileName << "\n";
    std::ifstream inputFile;
    inputFile.open(fileName);
    if (inputFile.is_open()) {
        inputFile >> n;
        std::cout << "Number of cities = " << n << "\n";
        if (inputFile.eof()) {
            return;
        }
        
        // Read and create the cities
        for (i = 1;  i <= n;  i++) {
            inputFile >> idx >> x >> y;
            if (inputFile.eof()) {
                break;
            }
            ++numCitiesRead;
            
            City *c = createCity(idx, x, y, false);
            unvisitedCities.push_back(c);
        }
    }
    inputFile.close();
    std::cout << "Number of cities actually read = " << numCitiesRead << "\n";
    if (numCitiesRead != n) {
        std::cout << "ERROR: number of cities read doesn't match expected city count!\n";
    }
    std::cout << "-------------------------\n";
}


//----------------------------------------------------------------------------------------------------------------

double heuristicTSPTour() {
    double minCostTour = 0.0;
    unsigned int numUnvisited = n;
    City *firstCity = NULL;
    City *lastCity  = NULL;
    
    auto totalt1 = std::chrono::high_resolution_clock::now();
    std::cout << "Starting the heuristic TSP solution\n";
    
    // The first city is the start of the tour
    firstCity = lastCity = unvisitedCities.front();
    unvisitedCities.pop_front();
    tourCities.push_back(firstCity);
    numUnvisited--;
    
    while (numUnvisited > 0) {
        if (dbg || ((numUnvisited % 1000) == 0)) std::cout << "  Cities left = " << numUnvisited << ", minCostTour = " << minCostTour << "\n";
        assert(numUnvisited == unvisitedCities.size());
        auto t1 = std::chrono::high_resolution_clock::now();
        
        // Find the closest city the tour hasn't visited yet. On a tie, use the closest city with the lowest index.
        // Ugly O(n) linear search, but I can't think of a better way...
        City    *bestCity = NULL;
        double   bestDist = DBL_MAX;
        uint32_t bestIdx  = UINT32_MAX;
        std::list<City*>::iterator bestIt;
        for (std::list<City*>::iterator it = unvisitedCities.begin();  it != unvisitedCities.end();  ++it) {
            City *candidate = *it;
            assert(!candidate->alreadyVisited);
            //double distToCandidateSq = (pow((candidate->x - lastCity->x), 2.0) + pow((candidate->y - lastCity->y), 2.0));
            double diffx = (candidate->x - lastCity->x);
            double diffy = (candidate->y - lastCity->y);
            double distToCandidateSq = (diffx*diffx + diffy*diffy);
            if (distToCandidateSq < bestDist) {
                if (dbg) std::cout << "      Better city " << candidate->index << " with dist " << distToCandidateSq << "\n";
                bestCity = candidate;
                bestDist = distToCandidateSq;
                bestIdx  = candidate->index;
                bestIt   = it;
            }
        }
        
        // Add the closest city to the tour cities
        assert(bestCity != NULL);
        lastCity = createCity(bestCity->index, bestCity->x, bestCity->y, /*visited*/ true);
        tourCities.push_back(lastCity);
        minCostTour += sqrt(bestDist);
        
        // Remove (and delete) the closest city from the unvisitedCities list
        bestIt = unvisitedCities.erase(bestIt);
        numUnvisited--;
        
        if (dbg) {
            auto t2 = std::chrono::high_resolution_clock::now();
            std::cout << "    Finding best city " << bestIdx << " with closest distance^2 " << bestDist
                      << " took " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << " ms\n";
            std::cout << "      Min tour cost so far = " << minCostTour << "\n";
        }
    }
    assert(unvisitedCities.size() == 0);
    
    // Now that every city has been visited, return to the first city to complete the tour.
    //float distToFirstCitySq = (pow((lastCity->x - firstCity->x), 2.0) + pow((lastCity->y - firstCity->y), 2.0));
    double lastDiffx = (lastCity->x - firstCity->x);
    double lastDiffy = (lastCity->y - firstCity->y);
    double distToFirstCitySq = (lastDiffx*lastDiffx + lastDiffy*lastDiffy);
    minCostTour += sqrt(distToFirstCitySq);
    if (dbg) std::cout << "\nAdded dist " << sqrt(distToFirstCitySq) << " from city " << lastCity->index
                       << " back to first city " << firstCity->index << "\n";
    
    auto totalt2 = std::chrono::high_resolution_clock::now();
    std::cout << "Min TSP tour distance is " << minCostTour << "\n";
    std::cout << "Total TSP computation took " << std::chrono::duration_cast<std::chrono::milliseconds>(totalt2-totalt1).count() << " ms\n";
    return minCostTour;
}


//----------------------------------------------------------------------------------------------------------------

int main(int argc, const char * argv[]) {
    if (argc < 2) {
        std::cout << "Usage: htsp <input file>\n";
        return 1;
    }
    std::cout << "The traveling salesman problem, implemented using a heuristic algorithm\n";
    
    // Read the city coordinates from the input file.
    readInputFile(argv[1]);
    
    // Run the heuristic TSP algorithm and return the minimum cost of a traveling salesman tour.
    double minTourCost = heuristicTSPTour();
    
    // Report the minimum cost of a traveling salesman tour, rounded down to closest integer
    std::cout << "\nMinimum cost of a traveling salesman tour = " << floor(minTourCost) << " = " << uint64_t(floor(minTourCost)) << "\n";
    std::cout << "\nProgram finished\n";
    return 0;
}
