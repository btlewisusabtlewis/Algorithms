//
//  main.cpp
//  Huffman
//
//  Created by Brian Lewis on 1/19/17.
//  Copyright © 2017 Brian Lewis. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <queue>
#include <stack>
#include <functional>     // std::greater
#include <climits>
#include "assert.h"

// Symbols contain one or more points
struct Symbol {
    Symbol(unsigned symLabel, bool leaf, unsigned symWeight,
           Symbol *lhs = NULL, Symbol *rhs = NULL) :
        label(symLabel), weight(symWeight), isLeaf(leaf),
        left(lhs), right(rhs) { }
    unsigned label;      // the integer identifying the symbol
    unsigned weight;     // proportional to the symbol's probabiity
    bool     isLeaf;     // if a leaf, no children
    Symbol  *left;
    Symbol  *right;
};

struct LessThanByWeight {
    bool operator()(Symbol * const &s1, Symbol * const &s2) {
        // return "true" if "s1" is ordered before "s2"
        return s1->weight > s2->weight;
    }
};

Symbol *createLeafSym(unsigned symLabel, unsigned symWeight) {
    Symbol *s = new Symbol(symLabel, /*leaf*/ true, symWeight);
    return s;
}

Symbol *createMetaSym(unsigned symWeight, Symbol *lhs, Symbol *rhs) {
    Symbol *s = new Symbol(/*symLabel*/ UINT_MAX, /*leaf*/ false, symWeight, lhs, rhs);
    return s;
}

//--------------------------------------------------------------------------------------------------------

// The number of symbols. These symbols are merged over time to become meta-symbols of trees for child symbols
unsigned symbolCount = 0;

// A heap holding the symbols ordered by weight. top() returns the heap's min element.
std::priority_queue<Symbol*, std::vector<Symbol*>, LessThanByWeight> symbols;

// Current codeword length while traversing the Huffman binary tree
unsigned cwLength = 0;
unsigned minCodewordLen = 0;
unsigned maxCodewordLen = 0;

bool dbg = false;

//--------------------------------------------------------------------------------------------------------

// Read symbols (i.e., their weights) from the input file and add to the "symbols" heap
int readInputFile(const char *fileName) {
    // Read the point descriptions from the input file.
    std::cout << "Input file: " << fileName << "\n";
    std::ifstream inputFile;
    inputFile.open(fileName);
    if (inputFile.is_open()) {
        inputFile >> symbolCount;
        if (inputFile.eof()) {
            return -1;
        }
        std::cout << "Number of symbols = " << symbolCount << "\n";
        
        // Read the weights of each input symbol
        unsigned symWeight;
        if (dbg) std::cout << "Input symbol weights:\n";
        for (int i = 0;  i < symbolCount;  i++) {
            inputFile >> symWeight;
            if (inputFile.eof()) {
                return -1;  // unexpected EOF
            }
            if (dbg && (i < 20)) std::cout << "  (" << i << ") " << symWeight << "\n";
            Symbol *leaf = createLeafSym(i, symWeight);
            symbols.push(leaf);
        }
        if (dbg) std::cout << "\n";
        inputFile.close();
    }
    
    unsigned numSymbolsRead = (unsigned)symbols.size();
    std::cout << "Number of symbols actually read = " << numSymbolsRead << "\n";
    if (numSymbolsRead != symbolCount) {
        std::cout << "ERROR: number of symbols read don't match expected symbol count!\n";
    }
    std::cout << "-------------------------\n";
    return 0;
}

//--------------------------------------------------------------------------------------------------------

void assignCW(Symbol *s) {
    assert(s != NULL);
    if (s->isLeaf) {
        minCodewordLen = std::min(cwLength, minCodewordLen);
        maxCodewordLen = std::max(cwLength, maxCodewordLen);
        if (dbg && (s->label < 20)) {
            std::cout << "Symbol " << s->label << " has codeword of length " << cwLength << "\n";
            std::cout << "   minCWLen = " << minCodewordLen << ", maxCWLen = " << maxCodewordLen << "\n";
        }
    } else {
        if (s->left != NULL) {
            cwLength++;
            assignCW(s->left);
            cwLength--;
        }
        if (s->right != NULL) {
            cwLength++;
            assignCW(s->right);
            cwLength--;
        }
    }
}

void assignCodewords(Symbol *root) {
    cwLength = 0;
    maxCodewordLen = 0;
    minCodewordLen = UINT_MAX;
    assignCW(root);
}

void huffman() {
    // Do (symbolCount-1) merges of the symbols to create a binary tree with the symbols as leaves
    unsigned lastIter = (symbolCount - 1);
    for (unsigned iter = 0;  iter < lastIter;  iter++) {
        // Extract the two lowest-weight/frequency symbols from the heap
        Symbol *first = symbols.top();
        symbols.pop();     // remove that minimum symbol
        Symbol *second = symbols.top();
        symbols.pop();     // remove that next-to-minimum symbol
        // Create a new meta-symbol subtree containing the two symbols as leaves.
        // Its weight is the sum of their weights.
        unsigned symWeight = (first->weight + second->weight);
        Symbol *metaSym = createMetaSym(symWeight, first, second);
        symbols.push(metaSym);
    }
    
    // Do a depth-first traversal of the resulting Huffman symbol tree to assign codewords to symbols
    // (not actually implemented yet) and to find the minimum- and maximum-length codewords.
    Symbol *root = symbols.top();
    assignCodewords(root);
}

//--------------------------------------------------------------------------------------------------------

int main(int argc, const char * argv[]) {
    if (argc < 2) {
        std::cout << "Usage: huffman <input file>\n";
        return -1;
    }
    std::cout << "Linear-time Huffman coding using a heap\n";
    
    // Read the symbols (i.e., symbol weights) from the input file into the heap "symbols".
    readInputFile(argv[1]);
    
    // Run the linear-time Huffman encoding on the symbols
    huffman();
    std::cout << "\nMin codeword length = " << minCodewordLen << "\n";
    std::cout << "Max codeword length = "   << maxCodewordLen << "\n";
    
    std::cout << "\nProgram finished\n";
    return 0;
}
