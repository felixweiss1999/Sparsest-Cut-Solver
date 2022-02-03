#pragma once
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>
#include <boost/dynamic_bitset.hpp>
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
    dynamic_bitset<uint64_t>* forget_bitset;
    uint16_t* join_j;
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


struct NodeIndexPair {
    NodeIndexPair(size_t n, uint64_t vI): node(n), valueIndex(vI) {}
    size_t node;
    uint64_t valueIndex;
};