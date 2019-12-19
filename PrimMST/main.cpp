//
//  main.cpp
//  PrimMST
//
//  Created by Brian Lewis on 1/7/17.
//  Copyright © 2017 Brian Lewis. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <set>
#include <algorithm>    // for std::min
#include <climits>
#include "assert.h"

struct Vertex;

// Edges refer to the Vertices they connect.
struct Edge
{
    int     label;       // the integer identifying the edge
    Vertex *lhs;
    Vertex *rhs;
    int     cost;
};

// Vertices refer to the edges incident on them
struct Vertex
{
    int  label;          // the integer identifying the vertex
    std::vector<Edge*> edges;
    bool inX;            // whether in set X yet or not
};


// The number and list of edges read from the input file.
int edgeCount = 0;
std::list<Edge*> edges;

// The number and vector of vertices read from the input file.
int vertexCount = 0;
std::vector<Vertex*> vertices;

// Variables for Prim's MST algorithm
std::set<int>   X;       // initially {s} for arb. vertex s. Loop adds endpoint v in VMinusX for cheapest edge from u in X.
std::set<int>   VMinusX; // initially set to all vertices minus {s}, shrinks as elements move to processed set X
std::set<Edge*> T;       // computed minimum cost spanning tree (edges) [invariant: X is vertices spanned by tree-so-far T]


bool dbg = false;


Edge *createEdge(int eLabel, Vertex *fromVertex, Vertex *toVertex, int cost) {
    Edge *e = new Edge;
    e->label = eLabel;
    e->lhs   = fromVertex;
    e->rhs   = toVertex;
    e->cost  = cost;
    return e;
}

Vertex *createVertex(int vLabel) {
    Vertex *v = new Vertex;
    v->label = vLabel;
    v->inX = false;
    return v;
}

// Read edge descriptions from the input file
void readInputFile(const char *fileName) {
    int endPt1, endPt2, cost;
    int numEdgesRead = 0;
    int i;
    
    // Read the edge descriptions from the input file.
    std::cout << "Input file: " << fileName << "\n";
    std::ifstream inputFile;
    inputFile.open(fileName);
    if (inputFile.is_open()) {
        inputFile >> vertexCount >> edgeCount;
        if (inputFile.eof()) {
            return;
        }
        if (dbg) std::cout << "Number of vertices " << vertexCount << ", number of edges " << edgeCount << "\n";
        
        // Make room for (vertexCount+1) entries in the vertices vector, and initialize vertices to all "unset" entries.
        // The extra entry is because the vertices are labelled 1..vertexCount (not starting at 0).
        vertices.assign((vertexCount+1), NULL);
        for (i = 1;  i <= vertexCount;  i++) {
            vertices[i] = createVertex(-1);   // unset label!
        }

        for (i = 1;  i <= edgeCount;  i++) {
            inputFile >> endPt1 >> endPt2 >> cost;
            if (inputFile.eof()) {
                break;
            }
            numEdgesRead++;
            if (dbg && (numEdgesRead < 20)) std::cout << "  Read edge (" << endPt1 << ", " << endPt2 << ", " << cost << ")\n";
            
            // Update the edge's endpoint Vertex objects.
            Vertex *lhsVertex = vertices[endPt1];
            lhsVertex->label = endPt1;
            Vertex *rhsVertex = vertices[endPt2];
            rhsVertex->label = endPt2;
            // Now create and record a new Edge object.
            Edge *newEdge = createEdge(i, lhsVertex, rhsVertex, cost);
            edges.push_back(newEdge);
            // Add the new Edge to the edge lists of the endpoint vertices.
            lhsVertex->edges.push_back(newEdge);
            rhsVertex->edges.push_back(newEdge);
        }
    }
    inputFile.close();
    std::cout << "Number of edges actually read = " << numEdgesRead << "\n";
    if (numEdgesRead != edgeCount) {
        std::cout << "ERROR: number of edges read don't match expected edge count!\n";
    }
    std::cout << "Created input graph\n";
    std::cout << "-------------------------\n";
}

void displayGraph() {
    std::cout << "\nInput graph\n";
    std::cout << vertexCount << " Nodes:\n";
    for (int i = 1;  i <= vertexCount;  i++) {
        Vertex *v = vertices[i];
        std::cout << "  (" << i << "): #edges = " << v->edges.size() << "\n";
    }
    std::cout << edges.size() << " Edges:\n";
    for (std::list<Edge*>::const_iterator iter = edges.begin();  iter != edges.end();  ++iter) {
        Edge *e = *iter;
        Vertex *endPt1 = e->lhs;
        Vertex *endPt2 = e->rhs;
        std::cout << "  (" << e->label << "): (" << endPt1->label << ", " << endPt2->label << "), cost = " << e->cost << "\n";
    }
    std::cout << "-------------------------\n";
}

void sanityCheckEdgesAndVertices() {
    bool problem = false;
    int i;
    
    // Make sure we've created every vertex.
    for (i = 1;  i <= vertexCount;  i++) {
        if (vertices[i]->label == -1) {
            std::cout << "*** Vertex[" << i << "] has an unset label (not an edge endpoint).\n";
            problem = true;
        }
        if (vertices[i]->label != i) {
            std::cout << "*** Vertex[" << i << "] has an unexpected label " << vertices[i]->label << "\n";
            problem = true;
        }
    }
    
    // For every vertex v, if v has an incident edge e, e's lhs or rhs is v.
    std::vector<Edge*>::const_iterator iter;
    for (i = 1;  i <= vertexCount;  i++) {
        Vertex *v = vertices[i];
        for (iter = v->edges.begin();  iter != v->edges.end();  ++iter) {
            Edge *e = *iter;
            if ((e->lhs != v) && (e->rhs != v)) {
                std::cout << "*** Vertex[" << v->label << "]->edges has an edge e where e.lhs and e.rhs aren't the vertex!\n";
                problem = true;
            }
        }
    }
    if (!problem) {
        std::cout << "Edges and vertices check out okay\n";
    }
}


//----------------------------------------------------------------------------------------------------------------

// Run brute-force Prim's MST algorithm to find a minimum cost spanning tree T and return its total cost
long long int computeMST() {
    int i;
    
    // Arbitrarily add the starting vertex 1 to the set X.
    X.insert(1);
    vertices[1]->inX = true;
    // Initialize V-X to all vertices except for 1.  V-X holds the set of vertices remaining to be processed.
    for (i = 2;  i <= vertexCount;  i++) {
        VMinusX.insert(i);
    }
    
    // Main loop:
    //   Let e = (u, v) be the cheapest edge of G with u in (already processed) X, and v in X-V.
    //   Then add e to T and v to X.
    unsigned round = 0;
    while (X.size() < vertexCount) {
        assert(VMinusX.size() > 0);
        round++;
        if (dbg) std::cout << "Round " << round << ": X.size() = " << X.size() << "\n";
        
        // Find the cheapest edge e = (u,v) where u in X and v in V-X
        int minCost = INT_MAX;
        Edge *cheapestEdge = NULL;
        for (std::list<Edge*>::const_iterator iter = edges.begin();  iter != edges.end();  ++iter) {
            Edge *e = *iter;
            Vertex *endPt1 = e->lhs;
            Vertex *endPt2 = e->rhs;
#if 0
            if (dbg) std::cout << "    edge " << e->label << ": cost " << e->cost
                               << ", (" << endPt1->label << ", " << endPt2->label
                               << "), lhs " << (endPt1->inX? "inX" : "!inX")
                               <<  ", rhs " << (endPt2->inX? "inX" : "!inX") << "\n";
#endif
            if ((endPt1->inX && !endPt2->inX) || (endPt2->inX && !endPt1->inX)) { // edge from endPt in X to an endPt in V-X
                if (e->cost < minCost) {
                    minCost = e->cost;
                    cheapestEdge = e;
                }
            }
        }
        if (dbg) std::cout << "  Cheapest edge across cut: " << cheapestEdge->label << ", cost " << cheapestEdge->cost << "\n";
        
        // Add e to T and add v to X.
        T.insert(cheapestEdge);
        Vertex *v;
        if (cheapestEdge->lhs->inX && !cheapestEdge->rhs->inX) {   // Here (u,v) = (lhs,rhs)
            v = cheapestEdge->rhs;
            if (dbg) std::cout << "    (u,v) = (" << cheapestEdge->lhs->label << ", " << cheapestEdge->rhs->label;
        } else if (cheapestEdge->rhs->inX && !cheapestEdge->lhs->inX) {   // Here (u,v) = (rhs,lhs)
            v = cheapestEdge->lhs;
            if (dbg) std::cout << "    (u,v) = (" << cheapestEdge->rhs->label << ", " << cheapestEdge->lhs->label;
        } else {
            assert(false);
        }
        VMinusX.erase(v->label);
        X.insert(v->label);
        vertices[v->label]->inX = true;
        if (dbg) std::cout << "), X.size " << X.size() << ", VMinusX.size " << VMinusX.size() << "\n";
    }
    
    // Now compute the overall cost of the minimum spanning tree T we just found
    long long int overallCostOfT = 0;
    if (dbg) std::cout << "\nFound an MST T with edges:\n";
    for (std::set<Edge*>::iterator it = T.begin();  it != T.end();  ++it) {
        Edge *e = *it;
        overallCostOfT += e->cost;
        if (dbg) std::cout << "  " << e->label << ": cost " << e->cost << ", ("
                           << e->lhs->label << ", " << e->rhs->label << ")\n";
    }
    
    return overallCostOfT;
}


//----------------------------------------------------------------------------------------------------------------

int main(int argc, const char * argv[]) {
    long long int overallCostOfT = 0;
    
    if (argc < 2) {
        std::cout << "Usage: primmst <input file>\n";
        return -1;
    }
    std::cout << “Brute-force Prim's MST algorithm\n";
    
    // Read the edge descriptions from the input file.
    readInputFile(argv[1]);
    
    // Display the newly-created graph.
    if (dbg) displayGraph();
    
    // Doublecheck the newly-created graph data structures.
    sanityCheckEdgesAndVertices();
    
    // Run Prim's MST algorithm to find a minimum cost spanning tree T
    overallCostOfT = computeMST();
    
    // Report the overall cost of T
    std::cout << "\nOverall cost of the MST T = " << overallCostOfT << "\n";
    std::cout << "\nProgram finished\n";
    return 0;
}
