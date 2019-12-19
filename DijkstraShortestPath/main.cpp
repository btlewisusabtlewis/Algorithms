//
//  main.cpp
//  DijkstraShortestPath
//
//  Created by Brian Lewis on 12/14/16.
//  Copyright © 2016 Brian Lewis. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <set>
#include <sstream>
#include <string>
#include "assert.h"

#define DIST_IF_NO_PATH_FROM_1  1000000

struct Vertex;

// Edges point to the Vertices they connect.
struct Edge
{
    Vertex *from;
    Vertex *to;
    unsigned length;         // length or weight of the arc from "from" to "to"
};

// Vertices contain two edge vectors: one for the incoming edges and another for the outgoing edges.
struct Vertex
{
    unsigned label;          // the integer identifying the vertex
    std::vector<Edge*> inEdges;
    std::vector<Edge*> outEdges;
    unsigned distFrom1;
};


// The vector of vertices we read from the input file.
std::vector<Vertex*> vertices;
unsigned vertexCount = 0;

// The list of the edges read from the input file.
std::list<Edge*> edges;
unsigned edgeCount = 0;

// Variables for the Dijkstra algorithm
std::set<int>    X;         // initially empty, grows as vertices are processed
std::set<int>    VMinusX;   // initially set to all vertices, shrinks as elements move to processed set X
std::vector<int> A;         // computed shortest path distances for each vertex

// Vertices to report the shortest path from vertex 1
#define DIST_REPORT_ARRAY_SIZE 10
int verticesForDistReport[DIST_REPORT_ARRAY_SIZE] = {7,37,59,82,99,115,133,165,188,197};


Edge *createEdge(unsigned fromVertex, unsigned toVertex, int edgeLength) {
    Edge *e = new Edge;
    if ((fromVertex == 0) || (fromVertex > vertexCount)) {
        std::cout << "*** from vertex " << fromVertex << " is out of range [1.." << vertexCount << "]\n";
        e->from = NULL;
    } else {
        e->from = vertices[fromVertex];
    }
    if ((toVertex == 0) || (toVertex > vertexCount)) {
        std::cout << "*** to vertex " << toVertex << " is out of range [1.." << vertexCount << "]\n";
        e->to = NULL;
    } else {
        e->to = vertices[toVertex];
    }
    e->length = edgeLength;
    return e;
}

Vertex *createVertex(int vLabel) {
    Vertex *v = new Vertex;
    v->label = vLabel;
    v->distFrom1 = DIST_IF_NO_PATH_FROM_1;   // if there is a path from 1, this will be reset to actual shortest path from 1
    return v;
}


//----------------------------------------------------------------------------------------------------------------

// Read a [<outVertex> "," <edgeLength>] pair of an outgoing edge from inVertex and its edge weight.
// Returns 0 on success and -1 on failure.
int readOutEdge(std::istringstream &iss, int fromVertex) {
    int toVertex, edgeLength;
    char delimeter;    // should be a comma
    if (!(iss >> toVertex >> delimeter >> edgeLength)) {    // error
        return -1;
    }
    //std::cout << "  [" << toVertex << delimeter << edgeLength << "]\n";
    if (toVertex <= vertexCount) {
        vertices[toVertex]->label = toVertex;   // set the "to" vertex's label
    }
    // Now create and record a new Edge object.
    Edge *newEdge = createEdge(fromVertex, toVertex, edgeLength);
    edges.push_back(newEdge);
    edgeCount++;
    // Add the new Edge to the appropriate in and out edge lists of the endpoint vertices.
    if (fromVertex <= vertexCount) {
        vertices[fromVertex]->outEdges.push_back(newEdge);
    }
    if (toVertex <= vertexCount) {
        vertices[toVertex]->inEdges.push_back(newEdge);
    }
    return 0;
}

int readInputLine(char *inputLine) {
    std::istringstream iss(inputLine);
    int fromVertex;
    int res;
    if (!(iss >> fromVertex)) {    // error
        return -1;
    }
    //std::cout << "  From vertex " << fromVertex << ": \n";
    vertices[fromVertex]->label = fromVertex;   // set the "from" vertex's label
    do {
        res = readOutEdge(iss, fromVertex);
    } while (res == 0);
    //std::cout << "\n";
    return 0;
}

// Read edge descriptions from the input file
void readInputFile(const char *fileName) {
    // Read the edge descriptions from the input file.
    //std::cout << "Input file: " << fileName << "\n";
    std::ifstream inputFile;

    char inLine[512];
    inputFile.open(fileName);
    if (inputFile.is_open()) {
        // Count the number of input lines == number of input vertices
        while (!inputFile.eof()) {
            inputFile.getline(inLine, 512);
            if (!inputFile.eof()) {
                vertexCount++;
            }
        }
        inputFile.close();
        
        // Create the vertex vector
        // Make room for (maxVertexLabel+1) entries in the vertices vector, and initialize vertices to all "unset" entries.
        // The extra entry is because the vertices are labelled 1..maxVertexLabel (not starting at 0).
        vertices.assign((vertexCount+1), NULL);
        for (int i = 1;  i <= vertexCount;  i++) {
            Vertex *v = createVertex(-1);   // unset label!
            vertices[i] = v;
        }
        
        // Now read and scan the input lines to get edge information
        inputFile.open(fileName);
        assert(inputFile.is_open());
        while (!inputFile.eof()) {
            inputFile.getline(inLine, 512);
            if (!inputFile.eof()) {
                readInputLine(inLine);
            }
        }
        inputFile.close();
    }
    std::cout << "Number of vertices = " << vertexCount << ", number of edges = " << edgeCount << "\n";
    std::cout << "Created input graph\n";
    std::cout << "-------------------------\n";
}

void displayGraph() {
    std::cout << "\nInput graph\n";
    std::cout << vertexCount << " Vertices:\n";
    for (int i = 1;  i <= vertexCount;  i++) {
        Vertex *v = vertices[i];
        std::cout << "  (" << i << "): distFrom1 = " << v->distFrom1 << ", #in edges = " << v->inEdges.size() << ", #out edges = " << v->outEdges.size() << "\n";
    }
    std::cout << edges.size() << " Edges:\n";
    for (std::list<Edge*>::const_iterator iter = edges.begin();  iter != edges.end();  ++iter) {
        Edge *e = *iter;
        if (e != NULL) {
            if (e->from == NULL) {
                std::cout << "  (NULL, ";
            } else {
                std::cout << "  (" << e->from->label << ", ";
            }
            if (e->to == NULL) {
                std::cout << "NULL)\n";
            } else {
                std::cout << e->to->label << ")\n";
            }
        }
    }
    std::cout << "-------------------------\n\n";
}

void sanityCheckEdgesAndVertices() {
    bool problem = false;
    // Make sure we've created every vertex.
    for (int i = 1;  i <= vertexCount;  i++) {
        if (vertices[i]->label == -1) {
            std::cout << "*** Vertex[" << i << "] has an unset label (not an edge endpoint).\n";
            problem = true;
        }
        if (vertices[i]->label != i) {
            std::cout << "*** Vertex[" << i << "] has an unexpected label " << vertices[i]->label << "\n";
            problem = true;
        }
    }
    
    // For every vertex v,
    //    a) if v has an in-edge e, e's to == v, and
    //    b) if v has an out-edge e, e's from == v.
    std::vector<Edge*>::const_iterator iter;
    for (int i = 1;  i <= vertexCount;  i++) {
        Vertex *v = vertices[i];
        for (iter = v->inEdges.begin();  iter != v->inEdges.end();  ++iter) {
            Edge *e = *iter;
            if (e->to != v) {
                std::cout << "*** Vertex[" << v->label << "]->inEdges has an edge e where e.to isn't the vertex!\n";
                problem = true;
            }
        }
        for (iter = v->outEdges.begin();  iter != v->outEdges.end();  ++iter) {
            Edge *e = *iter;
            if (e->from != v) {
                std::cout << "*** Vertex[" << v->label << "]->outEdges has an edge e where e.to isn't the vertex!\n";
                problem = true;
            }
        }
    }
    
    if (!problem) {
        //std::cout << "Edges and vertices check out okay\n";
    }
}


//----------------------------------------------------------------------------------------------------------------

void dijkstraShortestPath() {
    unsigned i;
    // Initialize A
    A.resize((vertexCount+1), DIST_IF_NO_PATH_FROM_1);   // +1 since our vertices go from [1..vertexCount]
    // The starting vertex is 1. Since it is already "processed" with shortest path length 0,
    // add 1 to the set of processed vertices X.
    A.at(1) = 0;
    X.insert(1);
    // Initialize V-X to all vertices except for 1. V-X holds the set of vertices remaining to be processed.
    for (i = 2;  i <= vertexCount;  i++) {
        VMinusX.insert(i);
    }
    
    // Main loop: choose the vertex wstar in V-X that minimizes A[vstar]+length(vstar,wstar) for some
    // vertex vstar in X, then add wstar to X.
    unsigned round = 0;
    while (X.size() < vertexCount) {
        assert(VMinusX.size() > 0);
        unsigned minShortestPath = DIST_IF_NO_PATH_FROM_1;
        int vstar = 0;   // vertex in X that through wstar, has shortest path length from 1
        int wstar = 0;
        round++;
        //std::cout << "\nRound " << round << ": X.size() = " << X.size() << "\n";
        // From each vertex v in X, consider total path lengths to adjacent vertices w in V-X.
        for (std::set<int>::iterator xiter = X.begin();  xiter != X.end();  ++xiter) {
            int v = *xiter;
            assert((v > 0) && (v <= vertexCount));
            Vertex *vVertex = vertices.at(v);
            assert(vVertex != NULL);
            //std::cout << "   Examining edges from X vertex v=" << v << "\n";
            for (std::vector<Edge*>::iterator edgeiter = vVertex->outEdges.begin();
                    edgeiter != vVertex->outEdges.end();  ++edgeiter) {
                Edge *e = *edgeiter;
                assert(e != NULL);
                Vertex *toVertex = e->to;
                assert(toVertex != NULL);
                int w = toVertex->label;
                assert((w > 0) && (w <= vertexCount));
                //std::cout << "      Considering V-X vertex w=" << w << "\n";
                // if w is in VMinusX, compute path length to w through v
                std::set<int>::iterator witer = VMinusX.find(w);
                if (witer != VMinusX.end()) {
                    unsigned pathLenToW = (A.at(v) + e->length);
                    //std::cout << "         w in V-X, (1,w) length=" << pathLenToW << ", (v,w) length=" << e->length << "\n";
                    if (pathLenToW < minShortestPath) {   // a shorter path through v to a vertex in V-X
                        minShortestPath = pathLenToW;
                        vstar = v;
                        wstar = w;
                    }
                }
            }
        }
        assert((vstar > 0) && (vstar <= vertexCount));
        assert((wstar > 0) && (wstar <= vertexCount));
        // The edge (vstar,wstar) is the one that minimizes the path langth to a vertex in V-X.
        // Add wstar to X and set its shortest path length.
        std::set<int>::iterator it;
#if 0
        std::cout << "   Before addition, set X size=" << X.size() << ":\n";
        for (it = X.begin();  it != X.end();  ++it) {
            int x = *it;
            assert((x > 0) && (x <= vertexCount));
            std::cout << "      " << x << ": length=" << A.at(x) << "\n";
        }
        std::cout << "   Before addition, set V-X size=" << VMinusX.size() << ":\n";
        for (it = VMinusX.begin();  it != VMinusX.end();  ++it) {
            int x = *it;
            assert((x > 0) && (x <= vertexCount));
            std::cout << "      " << x << "\n";
        }
#endif
        X.insert(wstar);
#if 0
        std::cout << "   Added vertex " << wstar << " with min path length=" << minShortestPath
            << " to set X (now size=" << X.size() << "):\n";
        for (it = X.begin();  it != X.end();  ++it) {
            int x = *it;
            assert((x > 0) && (x <= vertexCount));
            std::cout << "      " << x << "\n";
        }
#endif
        A.at(wstar) = minShortestPath;
        VMinusX.erase(wstar);
        std::cout << "   Round " << round << ": added vertex " << wstar << " with min path length=" << minShortestPath
                  << " to set X (now size=" << X.size() << ")\n";
#if 0
        for (it = X.begin();  it != X.end();  ++it) {
            int x = *it;
            assert((x > 0) && (x <= vertexCount));
            std::cout << "      " << x << ": " << A.at(x) << "\n";
        }
#endif
    }
    
    // Report the min path lengths to the problem set's required 10 vertices
    std::cout << "\nMin path lengths for the required vertices\n";
    for (i = 0;  i < DIST_REPORT_ARRAY_SIZE;  i++) {
        int v = verticesForDistReport[i];
        std::cout << A.at(v);
        if (i < (DIST_REPORT_ARRAY_SIZE-1)) {
            std::cout << ",";
        }
    }
    std::cout << "\n";
}

//----------------------------------------------------------------------------------------------------------------

int main(int argc, const char * argv[]) {
    if (argc < 2) {
        std::cout << "Usage: dijkstra <input file>\n";
        return -1;
    }
    std::cout << "Dijkstra shortest path algorithm\n";
    
    // Read the edge descriptions from the input file.
    readInputFile(argv[1]);
    
    // Display the newly-created graph.
    //displayGraph();
    
    // Doublecheck the newly-created graph data structures.
    sanityCheckEdgesAndVertices();
    
    dijkstraShortestPath();
    
    std::cout << "\nProgram finished\n";
    return 0;
}

