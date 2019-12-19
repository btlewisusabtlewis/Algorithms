//
//  main.cpp
//  KosarajuSCCS
//
//  Created by Brian Lewis on 12/10/16.
//  Copyright © 2016 Brian Lewis. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
// #include <map>


struct Vertex;

// Edges point to the Vertices they connect.
struct Edge
{
    Vertex *from;
    Vertex *to;
};

// Vertices contain two edge vectors: one for the incoming edges and another for the outgoing edges.
struct Vertex
{
    int     label;          // the integer identifying the vertex
    std::vector<Edge*> inEdges;
    std::vector<Edge*> outEdges;
    
    bool    explored;       // Has this vertex been explored by DFS yet?
    Vertex *leader;         // Kosaraju 2nd pass: what vertex is the leader for the SCCS containing this vertex?
    int     finishingTime;  // Finishing time for vertices found in 1st pass,
};

struct SCCS
{
    Vertex *leader;
    int     size;           // number of members (avoids repeated calls to members.size())
    std::list<Vertex*> members;
};


// The list of the edges we read from the input file.
std::list<Edge*> edges;

// The vector of vertices from the input file.
std::vector<Vertex*> vertices;
int maxVertexLabel = 0;

// Variables for the Kosaraju algorithm
int globalTime = 0;           // used to set finishing times for each vertex in 1st pass
Vertex *sourceVertex = NULL;  // used to set SCCS leaders in 2nd pass

// Other variables
std::vector<Vertex*> vertexWithFinishingTime;  // maps finishing time to vertex with that time
std::list<SCCS> sccssFound;

Edge *createEdge(Vertex *fromVertex, Vertex *toVertex) {
    Edge *e = new Edge;
    e->from = fromVertex;
    e->to   = toVertex;
    return e;
}

Vertex *createVertex(int vLabel) {
    Vertex *v = new Vertex;
    v->label = vLabel;
    v->explored = false;
    v->leader = NULL;
    v->finishingTime = 0;
    return v;
}


// Read edge descriptions from the input file
void readInputFile(const char *fileName) {
    // A list of the edge descriptions read from the input file.
    std::list<std::pair<int,int>> edgeDescList;
    
    // Number of input edges
    int edgeCount = 0;
    
    // Read the edge descriptions from the input file.
    std::cout << "Input file: " << fileName << "\n";
    std::ifstream inputFile;
    int endPt1, endPt2;
    inputFile.open(fileName);
    if (inputFile.is_open()) {
        while (!inputFile.eof()) {
            inputFile >> endPt1 >> endPt2;
            //std::cout << "  Read edge (" << endPt1 << ", " << endPt2 << ")\n";
            edgeCount++;
            
            // Record the labels of the new edge's endpoint vertices
            std::pair<int,int> newEdgeDesc = std::make_pair(endPt1, endPt2);
            edgeDescList.push_back(newEdgeDesc);
            // Remember the high watermark of the vertex labels we read
            maxVertexLabel = std::max(endPt1, maxVertexLabel);
            maxVertexLabel = std::max(endPt2, maxVertexLabel);
        }
    }
    inputFile.close();
    std::cout << "Number of edges read = " << edgeCount << ", largest vertex label = " << maxVertexLabel << "\n";
    std::cout << "-------------------------\n";

    // Make room for (maxVertexLabel+1) entries in the vertices vector, and initialize vertices to all "unset" entries.
    // The extra entry is because the vertices are labelled 1..maxVertexLabel (not starting at 0).
    vertices.assign((maxVertexLabel+1), NULL);
    for (int i = 1;  i <= maxVertexLabel;  i++) {
        Vertex *v = createVertex(-1);   // unset label!
        vertices[i] = v;
    }
    
    // Iterate through the edge descriptions creating the Vector and Edge objects
    for (std::list<std::pair<int,int>>::iterator listIter = edgeDescList.begin();
         listIter != edgeDescList.end();  listIter++) {
        std::pair<int,int> edgeDesc = *listIter;
        int from = edgeDesc.first;
        int to   = edgeDesc.second;
        // Look up or create the edge's endpoint Vertex objects.
        Vertex *fromVertex = vertices[from];
        fromVertex->label = from;
        Vertex *toVertex = vertices[to];
        toVertex->label = to;
        // Now create and record a new Edge object.
        Edge *newEdge = createEdge(fromVertex, toVertex);
        edges.push_back(newEdge);
        // Add the new Edge to the appropriate in and out edge lists of the endpoint vertices.
        fromVertex->outEdges.push_back(newEdge);
        toVertex->inEdges.push_back(newEdge);
    }
    std::cout << "Created input graph\n";
}

void displayGraph() {
    std::cout << "\nInput graph\n";
    std::cout << maxVertexLabel << " Nodes:\n";
    for (int i = 1;  i <= maxVertexLabel;  i++) {
        Vertex *v = vertices[i];
        std::cout << "  (" << i << "): #in edges = " << v->inEdges.size() << ", #out edges = " << v->outEdges.size() << "\n";
    }
    std::cout << edges.size() << " Edges:\n";
    for (std::list<Edge*>::const_iterator iter = edges.begin();  iter != edges.end();  ++iter) {
        Edge *e = *iter;
        Vertex *from = e->from;
        Vertex *to   = e->to;
        std::cout << "  (" << from->label << ", " << to->label << ")\n";
    }
    std::cout << "-------------------------\n";
}

void sanityCheckEdgesAndVertices() {
    bool problem = false;
    // Make sure we've created every vertex.
    for (int i = 1;  i <= maxVertexLabel;  i++) {
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
    for (int i = 1;  i <= maxVertexLabel;  i++) {
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
        std::cout << "Edges and vertices check out okay\n";
    }
}

void resetExploredFlags() {
    for (int i = 1;  i <= maxVertexLabel;  i++) {
        Vertex *v = vertices[i];
        v->explored = false;
    }
}

// Recursively explore child vertices starting from i.
// If reversed, process child edges in reverse order (i.e., use inEdges not outEdges).
void DFS(Vertex *i, bool reversed) {
    bool isPass1 = reversed;
    Edge   *e;
    Vertex *j;
#ifdef RECURSIVE_DFS
    i->explored = true;
    i->leader   = sourceVertex;
    if (!isPass1) {
        // Pass 2: Record the fact that i is a member of the SCCS with leader sourceVertex
        if (sccssFound.size() == 0) {
            std::cout << "\n*** Empty sccsFound list for leader " << sourceVertex->label << "\n";
        }
        sccssFound.back().size++;
        sccssFound.back().members.push_back(i);
    }
    if (reversed) {
        for (std::vector<Edge*>::const_iterator iter = i->inEdges.begin();  iter != i->inEdges.end();  ++iter) {
            e = *iter;
            j = e->from;
            if (!j->explored) {
                DFS(j, reversed);
            }
        }
    } else {
        for (std::vector<Edge*>::const_iterator iter = i->outEdges.begin();  iter != i->outEdges.end();  ++iter) {
            e = *iter;
            j = e->to;
            if (!j->explored) {
                DFS(j, reversed);
            }
        }
    }
    globalTime++;
    i->finishingTime = globalTime;
    if (isPass1) {
        vertexWithFinishingTime[globalTime] = i;
    }
#else  // !RECURSIVE_DFS
    std::list<Vertex*> S;   // work stack to replace DFS recursion with iteration
    S.push_back(i);
    i->explored = true;
    while (!S.empty()) {
        // Examine vertex v at the top of the stack
        Vertex *v = S.back();
        v->leader = sourceVertex;
        // If v has an unexplored child, push it and repeat the outer while loop;
        bool childFound = false;
        if (reversed) {
            for (std::vector<Edge*>::const_iterator iter = v->inEdges.begin();  ((iter != v->inEdges.end()) && !childFound);  ++iter) {
                e = *iter;
                j = e->from;
                if (!j->explored) {
                    S.push_back(j);
                    j->explored = true;
                    childFound = true;
                }
            }
        } else {
            for (std::vector<Edge*>::const_iterator iter = v->outEdges.begin();  ((iter != v->outEdges.end()) && !childFound);  ++iter) {
                e = *iter;
                j = e->to;
                if (!j->explored) {
                    S.push_back(j);
                    j->explored = true;
                    childFound = true;
                }
            }
        }
        // v has no children left, pop v from S and assign it a finishing time and if pass 2, add v to the SCCS with leader sourceVertex
        if (!childFound) {
            S.pop_back();
            globalTime++;
            v->finishingTime = globalTime;
            if (isPass1) {
                vertexWithFinishingTime[globalTime] = v;
            } else {
                // Pass 2: Record the fact that v is a member of the SCCS with leader sourceVertex
                sccssFound.back().size++;
                sccssFound.back().members.push_back(v);
            }
        }
    }
#endif // !RECURSIVE_DFS
}

void kosarajuFirstPass() {
    // Set global variables for Kosaraju algorithm
    globalTime = 0;       // used to set vertex finishing times in 1st pass
    sourceVertex = NULL;  // used to set SCCS leaders in 2nd pass
    
    // Initialize the vector of vertices indexed by their finishing time:
    // i.e., the ith entry is the vertex with finishing time i.
    vertexWithFinishingTime.assign((maxVertexLabel+1), NULL);
    
    // Run DFS on reversed graph in decreasing order of vertex number
    resetExploredFlags();
    for (int i = maxVertexLabel;  i >= 1;  i--) {
        Vertex *v = vertices[i];
        if (!v->explored) {
            sourceVertex = v;   // not really necessary in 1st pass
            DFS(v, /*reversed*/ true);
        }
    }
    
#if 0
    // For each finishing time N downto 1, display the vertex with that finishing time.
    std::cout << "\nFinishing times in decreasing order and the corresponding vertices\n";
    for (int i = maxVertexLabel;  i >= 1;  i--) {  // i is a finishing time
        Vertex *v = vertexWithFinishingTime[i];
        std::cout << "  (" << i << "): " << v->label << "\n";
    }
#endif
}

// Compare SCCS records by size for sorting.
bool compareBySize(const SCCS &first, const SCCS &second) {
    return first.size > second.size;
}

void kosarajuSecondPass() {
    // Set global variables for Kosaraju algorithm
    globalTime = 0;       // used to set vertex finishing times in 1st pass
    sourceVertex = NULL;  // used to set SCCS leaders in 2nd pass
    
    // Run DFS on forward graph in decreasing order of finishing time
    resetExploredFlags();
    for (int i = maxVertexLabel;  i >= 1;  i--) {  // i is a finishing time
        Vertex *v = vertexWithFinishingTime[i];
        if (!v->explored) {
            sourceVertex = v;   // needed in 2nd pass, v is going to be a leader
            SCCS newSCCS;
            newSCCS.leader = sourceVertex;
            newSCCS.size = 0;
            sccssFound.push_back(newSCCS);
            DFS(v, /*reversed*/ false);
        }
    }

#if 0
    // Display each SCCS found
    std::cout << "\nSCCSs found:\n";
    int i = 0;
    for (std::list<SCCS>::iterator iter = sccssFound.begin();  iter != sccssFound.end();  iter++) {
        SCCS s = *iter;
        i++;
        std::cout << "  (" << i << "): leader " << s.leader->label << ", size " << s.size << ", members:\n";
        for (std::list<Vertex*>::iterator memberIter = s.members.begin();  memberIter != s.members.end();  memberIter++) {
            Vertex *m = *memberIter;
            std::cout << "    " << m->label << "\n";
        }
    }
#endif
    
    // Display the SCCSs sorted by size
    sccssFound.sort(compareBySize);
    std::cout << "\nSCCSs sorted by size:\n";
    int i = 0;
    for (std::list<SCCS>::iterator iter = sccssFound.begin();  iter != sccssFound.end();  iter++) {
        SCCS s = *iter;
        i++;
        if (i > 10 ) {
            break;
        }
        std::cout << "  (" << i << "): leader " << s.leader->label << ", size " << s.size << "\n";
    }
}

void kosarajuSCCS() {
    // Pass 1 of Kosaraju's algorithm: find finishing times.
    kosarajuFirstPass();
    // Pass 2 of Kosaraju's algorithm: find SCCSs by traversing graph in descending order of finishing times.
    kosarajuSecondPass();
}

int main(int argc, const char * argv[]) {
    if (argc < 2) {
        std::cout << "Usage: kosaraju <input file>\n";
        return -1;
    }
    std::cout << "Kosaraju linear-time SCCS discovery\n";
    
    // Read the edge descriptions from the input file.
    readInputFile(argv[1]);
    
    // Display the newly-created graph.
    //displayGraph();
    
    // Doublecheck the newly-created graph data structures.
    sanityCheckEdgesAndVertices();
    
    kosarajuSCCS();
    
    std::cout << "\nProgram finished\n";
    return 0;
}

