#pragma once
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <queue>
#include <mutex>
#include "atomicQueue.h"

using namespace boost;
#define INTRODUCE 1;
#define FORGET 2;
#define JOIN 3;
#define LEAF 4;

struct Node {
    int type;
    size_t specialVertex;
    size_t inducedSubgraphSize;
    std::vector<size_t> bag;
    size_t bagSize;
    uint32_t* values;
    int height;
    size_t parent;
    std::vector<size_t> children;
    bool visited;
};


struct DecompositionUtility {
    bool isNice;
    bool tablesInitialized;
    int numberBags;
    size_t root;
    int numberVerticesOfDecomposedGraph;
    int maxBagSize;
    int** adjacencyMatrix;
};


typedef adjacency_list<vecS, vecS, undirectedS, Node, no_property, DecompositionUtility> TreeDecomposition;


struct CopiedNode {
    int type;
    //size_t specialVertex;
    size_t inducedSubgraphSize;
    //std::vector<size_t> bag;
    size_t bagSize;
    int height;
    size_t parent;
    std::vector<size_t> children;

    bool visited;
};


typedef adjacency_list<vecS, vecS, undirectedS, CopiedNode, no_property, no_property> CopiedTreeDecomposition;


struct vertex_copier {
    TreeDecomposition& from;
    CopiedTreeDecomposition& to;

    void operator()(size_t input, size_t output) const {
        to[output] = { from[input].type, from[input].inducedSubgraphSize, from[input].bag.size(), from[input].height, from[input].parent, from[input].children, from[input].visited };
    }
};