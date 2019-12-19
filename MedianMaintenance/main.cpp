//
//  main.cpp
//  MedianMaintenance
//
//  Created by Brian Lewis on 12/19/16.
//  Copyright © 2016 Brian Lewis. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <functional>     // std::greater
#include <climits>
#include "assert.h"

bool dbg = false;

// Read integers from the input file, determine their running median, and return the median sum at eof.
uint64_t medianMaintenance(const char *fileName) {
    // The input number and count of numbers read so far
    int inputNumber;
    int numberCount = 0;
    
    // Two near-equal-sized heaps holding the numbers read so far.
    std::priority_queue<int, std::vector<int>, std::less<int>>    Hlow;  // top() returns heap's max
    std::priority_queue<int, std::vector<int>, std::greater<int>> Hhigh; // top() returns heap's min
    
    // Sum of the medians read so far
    uint64_t medianSum = 0;
    
    // Read the numbers from the input file, one at a time.
    std::cout << "Input file: " << fileName << "\n";
    std::ifstream inputFile;
    inputFile.open(fileName);
    if (inputFile.is_open()) {
        for (;;) {
            inputFile >> inputNumber;
            if (inputFile.eof()) {
                break;
            }

            numberCount++;
            if (dbg) std::cout << "Input number " << numberCount << ": " << inputNumber << "\n";
            
            // Add the input number to the appropriate heap
            unsigned HlowMax  = ((Hlow.size() > 0)?  Hlow.top()  : 0);
            unsigned HhighMin = ((Hhigh.size() > 0)? Hhigh.top() : INT_MAX);
            if (inputNumber <= HlowMax) {
                // Insert into Hlow
                Hlow.push(inputNumber);
                if (dbg) std::cout << "  Added to Hlow\n";
            } else if (inputNumber >= HhighMin) {
                // Insert into Hhigh
                Hhigh.push(inputNumber);
                if (dbg) std::cout << "               Added to Hhigh\n";
            } else {
                // Insert into Hlow (could equally be Hhigh)
                Hlow.push(inputNumber);
                if (dbg) std::cout << "  Added to Hlow (Hhigh would be okay too)\n";
            }
            
            // If the heaps became unbalanced, move the appropriate min/max from one to the other
            size_t HlowSize  = Hlow.size();
            size_t HhighSize = Hhigh.size();
            if (HlowSize > (HhighSize + 1)) {  // take Hlow's max and insert it into Hhigh
                HlowMax  = Hlow.top();
                Hlow.pop();     // remove that top element
                Hhigh.push(HlowMax);
                if (dbg) std::cout << "  Removed " << HlowMax << " from Hlow, inserted into Hhigh\n";
            } else if (HhighSize > (HlowSize + 1)) {
                HhighMin  = Hhigh.top();
                Hhigh.pop();    // remove that top element
                Hlow.push(HhighMin);
                if (dbg) std::cout << "               Removed " << HhighMin << " from Hhigh, inserted into Hlow\n";
            }
            
            HlowSize  = Hlow.size();
            HhighSize = Hhigh.size();
            HlowMax   = ((HlowSize  > 0)? Hlow.top()  : 0);
            HhighMin  = ((HhighSize > 0)? Hhigh.top() : INT_MAX);
            if (dbg) std::cout << "    Hlow size = " << HlowSize << ", max = " << HlowMax << "\n";
            if (dbg) std::cout << "               Hhigh size = " << HhighSize << ", min = " << HhighMin << "\n";
            assert((HlowSize == HhighSize) || (HlowSize == HhighSize+1) || (HhighSize == HlowSize+1));
            
            // Determine the median of the numbers read so far
            unsigned median;
            if ((numberCount % 2) == 0) {   // even
                median = HlowMax;
            } else {
                median = ((HlowSize > HhighSize)? HlowMax : HhighMin);
            }
            medianSum += median;
            if (dbg) std::cout << "  Median = " << median << ", medianSum = " << medianSum << "\n";
        }
    }
    inputFile.close();
    std::cout << "Total numbers read = " << numberCount << ", final medianSum = " << medianSum << "\n";
    return medianSum;
}

int main(int argc, const char * argv[]) {
    if (argc < 2) {
        std::cout << "Usage: median <input file>\n";
        return -1;
    }
    std::cout << "Median Maintenance Using Heaps\n";
    
    // Read and process numbers from the input file.
    uint64_t finalMedianSum = medianMaintenance(argv[1]);
    std::cout << "\nFinal median sum mod 10000 = " << (finalMedianSum % 10000) << "\n";
    
    std::cout << "\nProgram finished\n";
    return 0;
}
