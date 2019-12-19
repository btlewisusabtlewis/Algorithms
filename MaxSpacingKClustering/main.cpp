//
//  main.cpp
//  MaxSpacingKClustering
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
#include <algorithm>    // for std::min
#include <climits>
#include "assert.h"

struct Point;
struct Cluster;

// Edges (point pairs) refer to the points they connect.
struct Edge
{
    int    label;       // the integer identifying the edge
    Point *lhs;
    Point *rhs;
    int    cost;
};

// Points refer to the edges incident upon them
struct Point
{
    int  label;          // the integer identifying the point
    std::vector<Edge*> edges;
    Cluster *cluster;    // the cluster currently containing the point
};

// Clusters contain one or more points
struct Cluster {
    int label;          // the integer identifying the cluster
    std::list<Point*> points;
};


// The number of clusters we want to end up with
int k = 0;

// The number and list of edges (point pairs) read from the input file.
int edgeCount = 0;
std::vector<Edge> edges;

// The number and vector of points read from the input file.
int pointCount = 0;
std::vector<Point*> points;

// Variables for the max-spacing k-clustering algorithm
int clusterCount = 0;
std::list<Cluster*> clusters;   // initially one per point, but these are merged by the algorithm

bool dbg = true;


Point *createPoint(int pLabel, Cluster *containingCluster) {
    Point *p   = new Point;
    p->label   = pLabel;
    p->cluster = containingCluster;
    return p;
}

Cluster *createCluster(int cLabel) {
    Cluster *c = new Cluster;
    c->label = cLabel;
    return c;
}

// Read edge (point pair) descriptions from the input file
void readInputFile(const char *fileName) {
    int numEdgesRead = 0;
    int endPt1, endPt2, cost, i;
    
    // Read the edge (point pair) descriptions from the input file.
    std::cout << "Input file: " << fileName << "\n";
    std::ifstream inputFile;
    inputFile.open(fileName);
    if (inputFile.is_open()) {
        inputFile >> pointCount;
        if (inputFile.eof()) {
            return;
        }
        // The edgeCount is (n choose 2) = n(n-1)/2 since there is one edge (i,j) for each choice of 1≤i<j≤n for n=pointCount.
        edgeCount = (pointCount * (pointCount-1))/2;
        if (dbg) std::cout << "Number of points = clusters = " << pointCount << "\n";
        if (dbg) std::cout << "Number of edges should be " << edgeCount << "\n";
        
        // Create the points and initial clusters. The extra zero element for the points vector is because the points are
        // labelled 1..pointCount (not starting at 0).
        points.push_back(NULL);        // points[0] = NULL (it's unused)
        for (i = 1;  i <= pointCount;  i++) {
            Cluster *c = createCluster(i);
            Point   *p = createPoint(-1, c);  // -1 => unset label for now
            c->points.push_back(p);    // the new cluster c contains the new point p
            points.push_back(p);       // add p to the points array and...
            clusters.push_back(c);     // add c to the clusters list
            clusterCount++;
        }
        
        // Make room for edgeCount jobs in the edges vector; the +1 is because they are numbered from 1 on.
        edges.reserve(edgeCount+1);
        edges.push_back(Edge());      // use the default constructor for the unused zeroth element
        for (i = 1;  i <= edgeCount;  i++) {
            inputFile >> endPt1 >> endPt2 >> cost;
            if (inputFile.eof()) {
                break;
            }
            numEdgesRead++;
            if (dbg && (numEdgesRead < 20)) std::cout << "  Read point pair (" << endPt1 << ", " << endPt2 << "), cost " << cost << "\n";
            
            Point *lhs = points[endPt1];
            Point *rhs = points[endPt2];
            
            // Update the edge's endpoint Point objects.
            lhs->label = endPt1;
            rhs->label = endPt2;
            // Now create and record a new Edge object.
            edges.push_back(Edge());      // use the default constructor
            edges[i].label = i;
            edges[i].lhs   = lhs;
            edges[i].rhs   = rhs;
            edges[i].cost  = cost;
            // Add the new Edge to the edge lists of the endpoint points.
            lhs->edges.push_back(&edges[i]);
            rhs->edges.push_back(&edges[i]);
        }
    }
    inputFile.close();
    std::cout << "Number of edges actually read = " << numEdgesRead << "\n";
    if (numEdgesRead != edgeCount) {
        std::cout << "ERROR: number of edges read don't match expected edge count!\n";
    }
    if (clusterCount != pointCount) {
        std::cout << "ERROR: number of clusters " << clusterCount << " doesn't match point count " << pointCount << "!\n";
    }
    std::cout << "Created input graph\n";
    std::cout << "-------------------------\n";
}

void displayGraph() {
    int i;
    std::cout << "\nInput graph\n";
    std::cout << pointCount << " Nodes:\n";
    for (i = 1;  i <= pointCount;  i++) {
        Point *v = points[i];
        std::cout << "  (" << i << "): #edges = " << v->edges.size() << "\n";
    }
    std::cout << edges.size() << " Edges:\n";
    for (i = 1;  i <= edgeCount;  i++) {
        Edge &e = edges[i];
        Point *endPt1 = e.lhs;
        Point *endPt2 = e.rhs;
        std::cout << "  (" << e.label << "): (" << endPt1->label << ", " << endPt2->label << "), cost = " << e.cost << "\n";
    }
    std::cout << "-------------------------\n";
}

void sanityCheckPointPairsAndPoints() {
    bool problem = false;
    int i;
    
    // Make sure we've created every Point.
    for (i = 1;  i <= pointCount;  i++) {
        Point *p = points[i];
        if (p->label == -1) {
            std::cout << "*** Point[" << i << "] has an unset label (not an edge endpoint).\n";
            problem = true;
        }
        if (p->label != i) {
            std::cout << "*** Point[" << i << "] has an unexpected label " << points[i]->label << "\n";
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
    
    // For every Point v, if v has an incident edge e, e's lhs or rhs is v.
    std::vector<Edge*>::const_iterator iter;
    for (i = 1;  i <= pointCount;  i++) {
        Point *v = points[i];
        for (iter = v->edges.begin();  iter != v->edges.end();  ++iter) {
            Edge *e = *iter;
            if ((e->lhs != v) && (e->rhs != v)) {
                std::cout << "*** Point[" << v->label << "]->edges has an edge e where e.lhs and e.rhs aren't the Point!\n";
                problem = true;
            }
        }
    }
    if (!problem) {
        std::cout << "Edges and points check out okay\n";
    }
}


//----------------------------------------------------------------------------------------------------------------

// Compare edges by cost for sorting in increasing order
bool compareByCost(const Edge &first, const Edge &second) {
    return first.cost < second.cost;
}

int maxSpacingKClustering() {
    int maxSpacing = 0;
    Point *p, *q;
    
    // First, sort the edges (point pairs) in increasing order of cost (distance). Ignore first unused element.
    std::sort(edges.begin()+1, edges.end(), compareByCost);
    
    std::vector<Edge>::const_iterator iter = edges.begin()+1;  // Ignore first unused element.
    for (iter = edges.begin()+1;  iter != edges.end();  ++iter) {
        // Let (p, q) = closest pair of separated points (determines the current spacing).
        const Edge &e = *iter;
        p = e.lhs;
        q = e.rhs;
        // Are p and q separated?
        Cluster *pc = p->cluster;
        Cluster *qc = q->cluster;
        if (pc != qc) {
            if (e.cost < maxSpacing) {
                std::cout << "\n*** Found an edge with LESS cost " << e.cost << " than the largest so far "
                          << maxSpacing << "!\n";
                return 0;
            }
            maxSpacing = e.cost;
        
            if (clusterCount <= k) {
                break;
            }
            
            // Merge the clusters containing p & q into a single cluster, the one with a lower label.
            Cluster *lowc, *highc;
            if (pc->label < qc->label) {
                lowc  = pc;
                highc = qc;
            } else {
                lowc  = qc;
                highc = pc;
            }
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
                return 0;
            }
            clusters.erase(it);
            clusterCount--;
        }
    }
    
    return maxSpacing;
}


//----------------------------------------------------------------------------------------------------------------

int main(int argc, const char * argv[]) {
    int maxSpacing = 0;
    std::cout << "Max-spacing K-Clustering algorithm\n";
    
    if (argc < 3) {
        std::cout << "Usage: kcluster <num clusters> <input file>\n";
        return -1;
    }
    
    // Read the number of clusters k from the first argument
    std::istringstream ss(argv[1]);
    if (!(ss >> k)) {
        std::cout << "Invalid number " << argv[1] << '\n';
    }
    std::cout << "Number of clusters k = " << k << "\n";
    
    // Read the edge descriptions from the input file (second argument).
    readInputFile(argv[2]);
    
    // Display the newly-created graph.
    if (dbg) displayGraph();
    
    // Doublecheck the newly-created graph data structures.
    sanityCheckPointPairsAndPoints();
    
    // Execute the Max-spacing K-Clustering algorithm with target k clusters and get the max spacing found
    maxSpacing = maxSpacingKClustering();
    std::cout << "\nFor " << k << " clusters, the max spacing found was " << maxSpacing << "\n";
    
    std::cout << "\nProgram finished\n";
    return 0;
}
