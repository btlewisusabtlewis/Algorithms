//
//  main.cpp
//  BigKClustering
//
//  Created by Brian Lewis on 1/13/17.
//  Copyright © 2017 Brian Lewis. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <set>
#include <queue>          // std::queue
#include <algorithm>      // for std::binary_search
#include <climits>
#include "assert.h"

struct Cluster;

// Points refer to the edges incident upon them
struct Point
{
    unsigned label;      // the point's labelBits-bit integer label
    Cluster *cluster;    // the cluster currently containing the point
};

// Clusters contain one or more points
struct Cluster {
    unsigned label;      // the integer identifying the cluster
    std::list<Point*> points;
};


// Number of bits for the labels of the input nodes
unsigned labelBits = 0;

// The count and vector of point labels read from the input file. Note there may be duplicates.
unsigned numInputLabels = 0;
std::vector<unsigned> inputLabels;

// The input labels after sorting and duplicate removal
unsigned sortedLabelCount = 0;
std::vector<unsigned> sortedLabels;

// The vector of points created from the sorted/deduped input node labels.
std::vector<Point*> points;

// The vector of clusters. Initially one per sorted/deduped point label. These contain one or more points.
unsigned clusterCount = 0;
std::list<Cluster*> clusters;   // initially one per point, but these are merged by the algorithm

bool dbg = false;


Point *createPoint(unsigned pLabel, Cluster *containingCluster) {
    Point *p   = new Point;
    p->label   = pLabel;
    p->cluster = containingCluster;
    return p;
}

Cluster *createCluster(unsigned cLabel) {
    Cluster *c = new Cluster;
    c->label = cLabel;
    return c;
}

// Read point descriptions from the input file
int readInputFile(const char *fileName) {
    unsigned numLabelsRead = 0;
    
    // Read the point descriptions from the input file.
    std::cout << "Input file: " << fileName << "\n";
    std::ifstream inputFile;
    inputFile.open(fileName);
    if (inputFile.is_open()) {
        inputFile >> numInputLabels >> labelBits;
        if (inputFile.eof()) {
            return -1;
        }
        std::cout << "Number of points = initial clusters = " << numInputLabels << "\n";
        std::cout << "Number of bits per label = " << labelBits << "\n";
        
        // Read the labels of the input points, each a labelBits bit vector
        if (dbg) std::cout << "Input point labels:\n";
        for (int i = 0;  i < numInputLabels;  i++) {
            unsigned label = 0;
            for (int b = 0;  b < labelBits;  b++) {
                unsigned bit;
                inputFile >> bit;
                if (inputFile.eof()) {
                    return -1;  // unexpected EOF
                }
                label = ((label << 1) | bit);
            }
            if (dbg) std::cout << "  (" << i << ") " << label << "\n";
            inputLabels.push_back(label);       // add label to the inputLabels array
        }
        inputFile.close();
    }
    
    numLabelsRead = (unsigned)inputLabels.size();
    std::cout << "Number of point labels actually read = " << numLabelsRead << "\n";
    if (numLabelsRead != numInputLabels) {
        std::cout << "ERROR: number of node labels read don't match expected node count!\n";
    }
    std::cout << "-------------------------\n";
    return 0;
}
    
        
//----------------------------------------------------------------------------------------------------------------

void radixSort() {
    std::queue<unsigned> zeroQ;
    std::queue<unsigned> oneQ;
    // Working copy of the input labels: holds partially-sorted labela
    std::vector<unsigned> workingVec;
    std::vector<unsigned>::iterator it;
    
    // Initialize the workingVec as a copy of the input labels
    workingVec.insert(workingVec.end(), inputLabels.begin(), inputLabels.end());
    
    // Radix sort: Iterate over the bits of the labels from the LSB to the MSB.
    for (int b = 0;  b < labelBits;  b++) {
        // Enqueue the labels based on their bth bit into one of the two queues.
        for (int i = 0;  i < numInputLabels;  i++) {
            unsigned l = workingVec[i];
            unsigned bit = ((l >> b) & 0x1);
            if (bit == 0) {
                zeroQ.push(l);
            } else {
                oneQ.push(l);
            }
        }
        // Dequeue the queues into the workingVec in increasing order of the bits processed so far
        workingVec.clear();
        while (!zeroQ.empty()) {
            workingVec.push_back(zeroQ.front());
            zeroQ.pop();
        }
        while (!oneQ.empty()) {
            workingVec.push_back(oneQ.front());
            oneQ.pop();
        }
        if (dbg) {
            std::cout << "  Pass " << b << ": partially-sorted labels are:\n";
            std::cout << "    ";
            for (it = workingVec.begin();  it != workingVec.end();  ++it) {
                unsigned l = *it;
                std::cout << l << " ";
            }
            std::cout << "\n";
        }
    }
    
    // Copy the now-sorted labels to sortedLabels and eliminate duplicates.
    bool firstTime = true;
    unsigned prevLabel = UINT_MAX;
    for (it = workingVec.begin();  it != workingVec.end();  ++it) {
        unsigned l = *it;
        if (firstTime) {
            sortedLabels.push_back(l);
            firstTime = false;
        } else {
            if (l != prevLabel) {
                sortedLabels.push_back(l);
            }
        }
        prevLabel = l;
    }
    sortedLabelCount = (unsigned)sortedLabels.size();
    workingVec.clear();  // no longer needed
    if (dbg) {
        std::cout << "The " << sortedLabelCount << " sorted and deduped labels:\n";
        std::cout << "    ";
        for (it = sortedLabels.begin();  it != sortedLabels.end();  ++it) {
            unsigned l = *it;
            std::cout << l << " ";
        }
        std::cout << "\n";
    }
}


//----------------------------------------------------------------------------------------------------------------

// Create points and clusters from the sorted and deduped input node labels
void createClusters() {
    // Create the points and initial clusters.
    for (int i = 0;  i < sortedLabelCount;  i++) {
        unsigned l = sortedLabels[i];
        Cluster *c = createCluster(i);
        clusters.push_back(c);            // add c to the clusters list
        clusterCount++;

        Point *p = createPoint(l, c);
        points.push_back(p);              // add p to the points vector
        c->points.push_back(p);           // the new cluster c contains the new point p
    }
    assert(clusters.size() == clusterCount);
    assert(sortedLabelCount == clusterCount);
    std::cout << "Created " << clusterCount << " initial clusters from sorted/deduped labels\n";
    std::cout << "-------------------------\n";
}

void displayPointsAndClusters() {
    int i;
    std::cout << "\nInput points and clusters\n";
    std::cout << sortedLabelCount << " Sorted/deduped Points:\n";
    for (i = 0;  i < sortedLabelCount;  i++) {
        Point *p = points.at(i);  //[i];
        std::cout << "  (" << i << "): label = " << p->label << ", cluster = " << p->cluster->label << "\n";
    }
    std::cout << clusters.size() << " Clusters:\n";
    i = 0;
    for (std::list<Cluster*>::iterator it = clusters.begin();  it != clusters.end();  ++it) {
        Cluster *c = *it;
        i++;
        std::cout << "  (" << i << "): label = " << c->label << ", points:\n";
        for (std::list<Point*>::iterator pit = c->points.begin();  pit != c->points.end();  ++pit) {
            Point *p = *pit;
            std::cout << "    label = " << p->label << ", cluster = " << p->cluster->label << "\n";
        }
    }
    std::cout << "-------------------------\n";
}

void sanityCheckPointPointsAndClusters() {
    bool problem = false;
    int i;
    
    // Make sure we've created every Point.
    for (i = 0;  i < sortedLabelCount;  i++) {
        Point *p = points.at(i);  //[i];
        if (p->label == -1) {
            std::cout << "*** Point[" << i << "] has an unset label (not an edge endpoint).\n";
            problem = true;
        }
        if (p->cluster == NULL) {
            std::cout << "*** Point[" << i << "] has an unset cluster\n";
            problem = true;
        }
        std::list<Point*> cPoints = p->cluster->points;
        std::list<Point*>::iterator res = std::find(cPoints.begin(), cPoints.end(), p);
        if (res == cPoints.end()) {
            std::cout << "*** Point[" << i << "] has the cluster " << p->cluster->label
                      << " whose point list doesn't include the point\n";
            problem = true;
        }
    }
    if (!problem) {
        std::cout << "Points and clusters check out okay\n";
    }
}


//----------------------------------------------------------------------------------------------------------------

bool lookupLabel(unsigned l) {
    bool found = std::binary_search(sortedLabels.begin(), sortedLabels.end(), l);
    //if (dbg) std::cout << "      Looking up label " << l << ": " << (found? "found\n" : "not found\n");
    //if (dbg && found) std::cout << "      Found label " << l << "\n";
    return found;
}

// Merge the clusters pc & qc into a single cluster, the one with a lower label.
void mergeClusters(Cluster *pc, Cluster *qc) {
    Cluster *lowc, *highc;
    if (pc->label < qc->label) {
        lowc  = pc;
        highc = qc;
    } else {
        lowc  = qc;
        highc = pc;
    }
    if (dbg) std::cout << "    Merging clusters " << lowc << "(" << lowc->points.size() << ") and " << highc << "(" << highc->points.size() << ")\n";
    // Set the cluster for highc's points to lowc, then merge them into lowc's point list.
    for (std::list<Point*>::iterator pit = highc->points.begin();  pit != highc->points.end();  ++pit) {
        Point *pInHighC = *pit;
        pInHighC->cluster = lowc;
    }
    lowc->points.splice(lowc->points.begin(), highc->points);
    // Now delete the cluster with higher label from the clusters list
    std::list<Cluster*>::iterator it = std::find(clusters.begin(), clusters.end(), highc);
    if (it == clusters.end()) {   // not found
        std::cout << "\n*** Didn't find cluster " << highc->label << " in clusters list!\n";
        return;
    }
    clusters.erase(it);
    clusterCount--;
}

int bigKClustering() {
    // At this point, we've sorted and deduped the input labels. There are sortedLabels of these.
    // And we've created the resulting sortedLabelCount points and clusters.
    // The variable clusterCount is the initial number of clusters, each with one point having a unique label.
    // Note that clusterCount may be smaller than the number of original labels if there were duplicates.
    
    // Scan through the unique (sorted) labels, looking for ones that are within Hamming distance 2 of
    // other labels. This can be checked by a binary search using lookupLabel().
    // These clusters of these nearby labels, if separated, should be merged.
    if (dbg) std::cout << "\nDoing clustering on " << sortedLabelCount << " nodes\n";
    if (dbg) std::cout << "Initial clusterCount = " << clusterCount << "\n";
    for (int i = 0;  i < sortedLabelCount;  i++) {
        Point *p = points.at(i);  //points[i];
        unsigned l = p->label;
        assert(l == sortedLabels.at(i));
        if (dbg) std::cout << "\n  Looking for cousin labels close to label " << l << "\n";
        
        // Are there any labels that are Hamming distance 1 away from l? This is a single-bit difference.
        // There are labelBits of these. Call these "cousins".
        for (int b = 0;  b < labelBits;  b++) {
            unsigned cousin = (l ^ (0x1 << b));
            if (lookupLabel(cousin)) {  // this close "cousin" label is among the input labels
                std::vector<unsigned>::iterator low = std::lower_bound(sortedLabels.begin(), sortedLabels.end(), cousin);
                assert(low != sortedLabels.end());
                unsigned idx = (unsigned)(low - sortedLabels.begin());
                //if (dbg) std::cout << "    Found cousin label " << cousin << " at index " << idx << "\n";
                Point *q = points.at(idx);  //[cousin];
                // Are p and q separated?
                Cluster *pc = p->cluster;
                Cluster *qc = q->cluster;
                if (pc != qc) {
                    // Merge the clusters con  taining p & q into a single cluster, the one with a lower label.
                    if (dbg) std::cout << "  Found dist-1 cousin label " << cousin << "\n";
                    mergeClusters(pc, qc);
                }
            }
        }
        
        // Are there any labels that are Hamming distance 2 away from l? This is a double-bit difference.
        // There are (labelBits choose 2) of these. Call these "cousins" also.
        for (int b1 = 0;  b1 < labelBits;  b1++) {
            unsigned cousin1 = (l ^ (0x1 << b1));
            for (int b2 = (b1+1);  b2 < labelBits;  b2++) {
                unsigned cousin = (cousin1 ^ (0x1 << b2));
                if (lookupLabel(cousin)) {  // this close "cousin" label is among the input labels
                    std::vector<unsigned>::iterator low = std::lower_bound(sortedLabels.begin(), sortedLabels.end(), cousin);
                    assert(low != sortedLabels.end());
                    unsigned idx = (unsigned)(low - sortedLabels.begin());
                    //if (dbg) std::cout << "    Found cousin label " << cousin << " at index " << idx << "\n";
                    Point *q = points.at(idx);  //[cousin];
                    // Are p and q separated?
                    Cluster *pc = p->cluster;
                    Cluster *qc = q->cluster;
                    if (pc != qc) {
                        // Merge the clusters containing p & q into a single cluster, the one with a lower label.
                        if (dbg) std::cout << "    Found dist-2 cousin label " << cousin << "\n";
                        mergeClusters(pc, qc);
                    }
                }
            }
        }
    }
    return clusterCount;
}


//----------------------------------------------------------------------------------------------------------------

int main(int argc, const char * argv[]) {
    int maxClusters = 0;
    std::cout << "Big and Implicit Max-spacing K-Clustering algorithm\n";
    
    if (argc < 2) {
        std::cout << "Usage: bigkcluster <input file>\n";
        return -1;
    }
    
    // Read the edge descriptions from the input file (first argument).
    readInputFile(argv[1]);
    
    // Sort and dedup the input labels, and create the points and clusters
    radixSort();
    createClusters();
    
    // Display the newly-created points and clusters.
    if (dbg) displayPointsAndClusters();
    
    // Doublecheck the newly-created graph data structures.
    sanityCheckPointPointsAndClusters();
    
    // Execute the implicit max-spacing K-clustering algorithm and get the max clusters k found with Hamming distance <= 2.
    maxClusters = bigKClustering();
    std::cout << "\nFor Hamming distance <= 2, max clusters = " << maxClusters << "\n";
    std::cout << "\nProgram finished\n";
    return 0;
}
