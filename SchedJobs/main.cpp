//
//  main.cpp
//  SchedJobs
//
//  Created by Brian Lewis on 1/7/17.
//  Copyright © 2017 Brian Lewis. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>

// Comment this line out to use the nonoptimal (weight-length) instead of the optimal (weight/length) job scores
//#define USE_RATIO_SCORE  1

#define DBG_PRINT       10


struct Job
{
    int weight;
    int length;
#ifdef USE_RATIO_SCORE
    float score;            // = (weight/length)
#else
    int score;            // = (weight - length)
#endif
    int completionTime;
};

// The number and vector of jobs read from the input file.
int numJobs = 0;
std::vector<Job> jobs;

bool dbg = true;


// Compare jobs by (weight - length) for sorting in decreasing order.
bool sortByReverseScore(const Job &lhs, const Job &rhs) {
    bool comp;
    if (lhs.score == rhs.score) {
        // Break tie in favor of the job with higher weight
        comp = (lhs.weight > rhs.weight);
    } else {
        comp = (lhs.score > rhs.score);
    }
    return comp;
}

void schedule(const char *fileName) {
    int weight, length;
    int jobsRead = 0;
    int i;
    
    // Read the jobs from the input file, one at a time, into the vector jobs.
    std::cout << "Input file: " << fileName << "\n";
    std::ifstream inputFile;
    inputFile.open(fileName);
    if (inputFile.is_open()) {
        inputFile >> numJobs;
        std::cout << "Number of input jobs = " << numJobs << "\n";
        
        // Make room for numJobs jobs in the jobs vector
        jobs.reserve(numJobs+1);
        for (i = 0;  i <= numJobs;  i++) {
            inputFile >> weight >> length;
            if (inputFile.eof()) {
                break;
            }
            jobsRead++;
            jobs.push_back(Job());      // use the default constructor
            jobs[i].weight = weight;
            jobs[i].length = length;
#ifdef USE_RATIO_SCORE
            jobs[i].score  = (float(weight) / float(length));
#else
            jobs[i].score  = (weight - length);
#endif
            jobs[i].completionTime = 0;
            
            if (dbg && (i < DBG_PRINT)) {
                std::cout << "  Job " << i << ": weight " << jobs[i].weight << ", length "
                          << jobs[i].length << ", score " << jobs[i].score << "\n";
            }
        }
    }
    inputFile.close();
    if (dbg) std::cout << "Total jobs actually read = " << jobsRead << "\n\n";
    
    // Sort the jobs by decreasing order of the difference (weight - length)
    std::sort(jobs.begin(), jobs.end(), sortByReverseScore);
    
    // Now assign completion times to each job and total the sum of those times
    int prevFinishTime = 0;
    long long int sumOfWeightedCompTimes = 0;
    for (i = 0;  i <= numJobs;  i++) {
        if (dbg && (i < DBG_PRINT)) {
            std::cout << "  Job " << i << ": weight " << jobs[i].weight << ", length "
                      << jobs[i].length << ", score " << jobs[i].score << ", prevFinishT " << prevFinishTime << "\n";
        }
        
        jobs[i].completionTime = (prevFinishTime + jobs[i].length);
        prevFinishTime = jobs[i].completionTime;
        
        // Add j's weighted completion time to the running total
        sumOfWeightedCompTimes += (jobs[i].weight * jobs[i].completionTime);
        
        if (dbg && (i < DBG_PRINT)) {
            std::cout << "    j->compTime " << jobs[i].completionTime << ", sumWCompT " << sumOfWeightedCompTimes << "\n";
        }
    }
    
    std::cout << "\nResult sumOfWeightedCompTimes = " << sumOfWeightedCompTimes << "\n";
}



int main(int argc, const char * argv[]) {
    if (argc < 2) {
        std::cout << "Usage: schedjobs <input file>\n";
        return -1;
    }
    std::cout << "Schedule jobs using a greedy algorithm\n";
    
    schedule(argv[1]);
    return 0;
}
