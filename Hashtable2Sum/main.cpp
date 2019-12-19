//
//  main.cpp
//  Hashtable2Sum
//
//  Created by Brian Lewis on 12/24/16.
//  Copyright © 2016 Brian Lewis. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <unordered_set>   // for hashtable with no mapping: we want to insert and test for existence of keys

#define NUM_VALUES             1000000
#define PRIME_ABOVE_NUM_VALUES 1000171

bool dbg = true;

void twoSum(const char *fileName) {
    // The count of input values read so far and the last input value
    int valueCount = 0;
    // Array holding the input numbers
    long long int A[NUM_VALUES];
    // Hashtable holding the input numbers
    std::unordered_set<long long int> H;
    // The number of input values t in [-10,000..10,000] such that there are distinct x,y where x+y = t.
    int numFoundTValues = 0;
    int i;
    
    // Read the numbers from the input file, one at a time, into the vector A.
    std::cout << "Input file: " << fileName << "\n";
    std::ifstream inputFile;
    inputFile.open(fileName);
    if (inputFile.is_open()) {
        for (;;) {
            inputFile >> A[valueCount];
            if (inputFile.eof()) {
                break;
            }
            //if (dbg) std::cout << "Input number " << valueCount << ": " << A[valueCount] << "\n";
            valueCount++;
        }
    }
    inputFile.close();
    std::cout << "Total numbers read = " << valueCount << "\n";
    
    // Insert vector A's elements (the input values) into the hashtable (unordered set) H
    H.rehash(PRIME_ABOVE_NUM_VALUES);    // set the number of H's buckets to PRIME_ABOVE_NUM_VALUES or larger
    for (i = 0;  i < NUM_VALUES;  i++) {
        H.insert(A[i]);
    }
    
    // Count the target values t in [-10,000..10,000] where there are distinct input values x,y where t = (x+y).
    long long int x, y, t;
    int tCount;
    for (t = -10000;  t <= 10000;  t++) {
        tCount = 0;
        for (i = 0;  i < NUM_VALUES;  i++) {
            x = A[i];
            y = (t-x);
            if (x != y) {
                std::unordered_set<long long int>::const_iterator fiter = H.find(y);
                if (fiter != H.end()) {   // a match
                    numFoundTValues++;
                    tCount++;
                    //if (dbg) std::cout << "  For t = " << t << " with x = " << x << " and y = " << y << ", x+y = t\n";
                }
            }
        }
        if (tCount > 0) {
            std::cout << "  For t = " << t << ", found " << tCount << " pairs\n";
        }
    }
    if (dbg) std::cout << "Found " << numFoundTValues << " t in [-10,000..10,000] st there are distinct x,y where x+y = t\n";
}

int main(int argc, const char * argv[]) {
    if (argc < 2) {
        std::cout << "Usage: twosum <input file>\n";
        return -1;
    }
    std::cout << "A modified 2-SUM algorithm\n";
    
    twoSum(argv[1]);
    return 0;
}
