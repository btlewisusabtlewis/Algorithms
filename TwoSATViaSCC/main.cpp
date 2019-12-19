//
//  main.cpp
//  TwoSATViaSCC
//
//  Created by Brian Lewis on 3/3/17.
//  Copyright © 2017 Brian Lewis. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include "assert.h"

// Uncomment to enable debugging output
//#define DEBUG_2SAT   1

// Uncomment to show more details about the found SCCSs
//#define SCCS_DETAILS 1

struct Vertex;

// Edges point to the Vertices they connect.
struct Edge
{
    int     label;          // the integer identifying the edge
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


// The number of clauses == number of variables
uint32_t numClauses   = 0;
uint32_t numVariables = 0;


// The number and list of the edges created for the implications added for each clause:
// That is, clause (lit1 or lit2) results in ~lit1 => lit2 and ~lit2 => lit1.
std::list<Edge*> edges;
uint32_t numEdges = 0;

// The number and vector of vertices created for each variable and its complement.
std::vector<Vertex*> vertices;
uint32_t numVertices = 0;


// Variables for the Kosaraju SCCS algorithm
int globalTime = 0;           // used to set finishing times for each vertex in 1st pass
Vertex *sourceVertex = NULL;  // used to set SCCS leaders in 2nd pass

// Other variables
std::vector<Vertex*> vertexWithFinishingTime;  // maps finishing time to vertex with that time
std::list<SCCS> sccssFound;


//----------------------------------------------------------------------------------------------------------------

Edge *createEdge(int eLabel, Vertex *fromVertex, Vertex *toVertex) {
    Edge *e = new Edge;
    e->label = eLabel;
    e->from  = fromVertex;
    e->to    = toVertex;
    return e;
}

Vertex *createVertex(int vLabel) {
    Vertex *v   = new Vertex;
    v->label    = vLabel;
    v->explored = false;
    v->leader   = NULL;
    v->finishingTime = 0;
    return v;
}

Edge *addEdge(Vertex *from, Vertex *to) {
    numEdges++;
    Edge *e = createEdge(numEdges, from, to);
    edges.push_back(e);
    // Add the new Edge e to the appropriate in and out edge lists of the endpoint vertices.
    from->outEdges.push_back(e);
    to->inEdges.push_back(e);
    return e;
}


//----------------------------------------------------------------------------------------------------------------

// The vertex for var itself is at vertices[2*var-1]: i.e., the 1st at 1, the second at 3, ...
int idxForVar(int var) {
    assert(var > 0);
    return (2*var - 1);
}

// The vertex for var's complement ~var is at vertices[2*var]: i.e., the 1st at 2, the second at 4, ...
int idxForComp(int var) {
    assert(var > 0);
    return (2*var    );
}

Vertex *litVertex(int lit) {
    int var, idx;
    if (lit > 0) {
        // the corresponding var is true, so varIdx == lit
        var = lit;
        idx = idxForVar(var);
    } else {
        // the corresponding var is false, so varIdx == -lit
        var = -lit;
        idx = idxForComp(var);
    }
    return vertices[idx];
}

Vertex *notLitVertex(int lit) {
    int var, idx;
    if (lit > 0) {
        // the corresponding var is true, so varIdx == lit, or (not lit) == var's complement
        var = lit;
        idx = idxForComp(var);
    } else {
        // the corresponding var is false, so varIdx == -lit and (not lit) == the var itself
        var = -lit;
        idx = idxForVar(var);
    }
    return vertices[idx];
}

void createImplicationEdges(int lit1, int lit2) {
    // Create two "implication" edges for each clause: i.e., clause (lit1 or lit2) results in
    // both ~lit1 => lit2 and ~lit2 => lit1 edges.
    Vertex *from, *to;
    Edge *e;
    assert((lit1 != 0) && (lit2 != 0));
    
    // Create the ~lit1 => lit2 edge.
    from = notLitVertex(lit1);
    to   = litVertex(lit2);
    e = addEdge(from, to);
#ifdef DEBUG_2SAT
    std::cout << "    Created edge " << e->label << " for " << from->label << " => " << to->label << "\n";
#endif
    
    // Create the ~lit2 => lit1 edge.
    from = notLitVertex(lit2);
    to   = litVertex(lit1);
    e = addEdge(from, to);
#ifdef DEBUG_2SAT
    std::cout << "    Created edge " << e->label << " for " << from->label << " => " << to->label << "\n";
#endif
}

// Read clauses from the input file
void readInputFile(const char *fileName) {
    std::cout << "Input file: " << fileName << "\n";
    std::ifstream inputFile;
    inputFile.open(fileName);
    if (inputFile.is_open()) {
        while (!inputFile.eof()) {
            inputFile >> numVariables;
            if (inputFile.eof()) {
                return;
            }
            numClauses = numVariables;
            std::cout << "Number of variables = number of clauses = " << numVariables << "\n";
            
            // Create the vertices vector
            // The "2*" is because we create two vertices for each variable v: v itself and the ~v complement
            numVertices = (2*numVariables);
            // The extra +1 entry is because the vertices are labelled 1..numVertices (not starting at 0).
            vertices.assign((numVertices + 1), NULL);
            for (int i = 1;  i <= numVariables;  i++) {
                Vertex *v = createVertex(i);
                vertices[idxForVar(i)]  = v;
                Vertex *vCompl = createVertex(-i);
                vertices[idxForComp(i)] = vCompl;
#ifdef DEBUG_2SAT
                std::cout << "  Created vertices " << v->label << " and " << vCompl->label
                          << " for variable " << i << " at idxs " << idxForVar(i) << " and " << idxForComp(i) << "\n";
#endif
            }
            
            for (int i = 1;  i <= numClauses;  i++) {
                int lit1, lit2;
                inputFile >> lit1 >> lit2;
                if (inputFile.eof()) {
                    return;
                }
#ifdef DEBUG_2SAT
#ifdef SCCS_DETAILS
                std::cout << "  Read clause (" << lit1 << " or " << lit2 << ")\n";
#else
                if (numClauses <= 10) std::cout << "  Read clause (" << lit1 << " or " << lit2 << ")\n";
#endif
#endif /* DEBUG_2SAT */
                numClauses++;
                // Create two "implication" edges for each clause
                createImplicationEdges(lit1, lit2);
            }
        }
    }
    inputFile.close();
    std::cout << "Number of clauses = " << numClauses << ", number of vertices = " << numVertices << "\n";
    std::cout << "Created input graph\n";
    std::cout << "-------------------------\n";
}

void displayGraph() {
#ifdef DEBUG_2SAT
    int i;
    std::cout << "\nInput graph\n";
    std::cout << numVertices << " Vertices:\n";
    for (i = 1;  i <= numVertices;  i++) {
        Vertex *v = vertices[i];
        std::cout << "  (" << i << "): variable " << v->label << ", in/out edges = " << v->inEdges.size() << "/" << v->outEdges.size() << "\n";
    }
    std::cout << edges.size() << " Edges:\n";
    i = 1;
    for (std::list<Edge*>::const_iterator iter = edges.begin();  iter != edges.end();  ++iter) {
        Edge *e = *iter;
        Vertex *from = e->from;
        Vertex *to   = e->to;
        std::cout << "  (" << i << "): (" << from->label << ", " << to->label << ")\n";
        i++;
    }
    std::cout << "-------------------------\n";
#endif
}

void sanityCheckEdgesAndVertices() {
    bool problem = false;
    // For every vertex v,
    //    a) if v has an in-edge e, e's to == v, and
    //    b) if v has an out-edge e, e's from == v.
    std::vector<Edge*>::const_iterator iter;
    for (int i = 1;  i <= numVertices;  i++) {
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
#ifdef DEBUG_2SAT
    if (!problem) {
        std::cout << "Graph checks out okay\n";
    }
#endif
    if (problem) {
        std::cout << "Graph check FAILED!\n";
    }
}


//----------------------------------------------------------------------------------------------------------------

void resetExploredFlags() {
    for (int i = 1;  i <= numVertices;  i++) {
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


//----------------------------------------------------------------------------------------------------------------

// Compare SCCS records by size for sorting.
bool compareBySize(const SCCS &first, const SCCS &second) {
    return first.size > second.size;
}

// For each finishing time N downto 1, display the vertex with that finishing time.
void displayVerticesByFinishingTime() {
#ifdef DEBUG_2SAT
    std::cout << "\nFinishing times in decreasing order and the corresponding vertices\n";
#ifdef SCCS_DETAILS
    for (int i = numVertices;  i >= 1;  i--) {  // i is a finishing time
        Vertex *v = vertexWithFinishingTime[i];
        std::cout << "  (" << i << "): " << v->label << "\n";
    }
#else
    int count = 0;
    for (int i = numVertices;  i >= 1;  i--) {  // i is a finishing time
        count++;
        if (count > 10) break;
        Vertex *v = vertexWithFinishingTime[i];
        std::cout << "  (" << i << "): " << v->label << "\n";
    }
#endif
#endif /* DEBUG_2SAT */
}

void displayFoundSCCSs() {
#ifdef DEBUG_2SAT
    int i = 0;
#ifdef SCCS_DETAILS
    // Display detailed information about ALL SCCSs found
    std::cout << "\nSCCSs found:\n";
    for (std::list<SCCS>::iterator iter = sccssFound.begin();  iter != sccssFound.end();  iter++) {
        SCCS s = *iter;
        i++;
        std::cout << "  (" << i << "): leader " << s.leader->label << ", size " << s.size << ", members:\n";
        for (std::list<Vertex*>::iterator memberIter = s.members.begin();  memberIter != s.members.end();  memberIter++) {
            Vertex *m = *memberIter;
            std::cout << "    " << m->label << "\n";
        }
    }
#else
    // Display the first 10 SCCSs sorted by size
    sccssFound.sort(compareBySize);
    std::cout << "\nSCCSs sorted by size:\n";
    for (std::list<SCCS>::iterator iter = sccssFound.begin();  iter != sccssFound.end();  iter++) {
        SCCS s = *iter;
        i++;
        if (i > 10 ) {
            break;
        }
        std::cout << "  (" << i << "): leader " << s.leader->label << ", size " << s.size << "\n";
    }
#endif
#endif /* DEBUG_2SAT */
}


//----------------------------------------------------------------------------------------------------------------

void kosarajuFirstPass() {
    // Set global variables for Kosaraju algorithm
    globalTime = 0;       // used to set vertex finishing times in 1st pass
    sourceVertex = NULL;  // used to set SCCS leaders in 2nd pass
    
    // Initialize the vector of vertices indexed by their finishing time:
    // i.e., the ith entry is the vertex with finishing time i.
    vertexWithFinishingTime.assign((numVertices+1), NULL);
    
    // Run DFS on reversed graph in decreasing order of vertex number
    resetExploredFlags();
    for (int i = numVertices;  i >= 1;  i--) {
        Vertex *v = vertices[i];
        if (!v->explored) {
            sourceVertex = v;   // not really necessary in 1st pass
            DFS(v, /*reversed*/ true);
        }
    }
    displayVerticesByFinishingTime();
}

void kosarajuSecondPass() {
    // Set global variables for Kosaraju algorithm
    globalTime = 0;       // used to set vertex finishing times in 1st pass
    sourceVertex = NULL;  // used to set SCCS leaders in 2nd pass
    
    // Run DFS on forward graph in decreasing order of finishing time
    resetExploredFlags();
    for (int i = numVertices;  i >= 1;  i--) {  // i is a finishing time
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
    displayFoundSCCSs();
}

void kosarajuSCCS() {
    // Pass 1 of Kosaraju's algorithm: find finishing times.
    kosarajuFirstPass();
    // Pass 2 of Kosaraju's algorithm: find SCCSs by traversing graph in descending order of finishing times.
    kosarajuSecondPass();
}


//----------------------------------------------------------------------------------------------------------------

bool isSatisfiable() {
    for (int i = 1;  i <= numVariables;  i++) {
        Vertex *var  = vertices[idxForVar(i)];
        Vertex *comp = vertices[idxForComp(i)];
        if (var->leader == comp->leader) {  // not satisfiable: var and its complement in same SCCS
            return false;
        }
    }
    return true;
}


//----------------------------------------------------------------------------------------------------------------

int main(int argc, const char * argv[]) {
    if (argc < 2) {
        std::cout << "Usage: 2SAT <input file>\n";
        return -1;
    }
    std::cout << "2SAT satisfiability check using SCCS reduction\n";
    
    // Read the edge descriptions from the input file.
    readInputFile(argv[1]);
    
    // Display the newly-created graph.
    displayGraph();
    
    // Doublecheck the newly-created graph data structures.
    sanityCheckEdgesAndVertices();
    
    // Find SCCSs of the graph
    kosarajuSCCS();
    
    if (isSatisfiable()) {
        std::cout << "2SAT instance is satisfiable\n";
    } else {
        std::cout << "2SAT instance is NOT satisfiable\n";
    }
    
    std::cout << "\nProgram finished\n";
    return 0;
}

