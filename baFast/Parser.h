#pragma once
#include "TreeDecomposition.h"


class Parser {
public:
	TreeDecomposition* parse(std::istream&);
	void print(TreeDecomposition&, TreeDecomposition::vertex_descriptor);
	void exportDimax(TreeDecomposition&, std::ostream&);
	void debugAlgorithm(TreeDecomposition&);
	void fillAdjacencyMatrix(TreeDecomposition&, std::istream&, int width = 0);
	
private:
	void makeNice(TreeDecomposition&);
	std::queue<size_t> HOLLOWmakeNice(TreeDecomposition&, TreeDecomposition::vertex_descriptor);
	void insertChainBetweenUnequalNodes(TreeDecomposition&, TreeDecomposition::vertex_descriptor, TreeDecomposition::vertex_descriptor);
	void HOLLOWinsertChainBetweenUnequalNodes(TreeDecomposition&, TreeDecomposition::vertex_descriptor, TreeDecomposition::vertex_descriptor);
	uint64_t calculateNumberOfOperations(TreeDecomposition&, TreeDecomposition::vertex_descriptor, std::queue<size_t>&);
	size_t calculateOptimalRoot(TreeDecomposition&);
	std::vector<size_t>& removeElementFromBag(std::vector<size_t>&, size_t);
	std::queue<size_t> split(TreeDecomposition&, TreeDecomposition::vertex_descriptor, int);
	void traverseUpThread(TreeDecomposition&, TreeDecomposition::vertex_descriptor);
	std::string tableOfNode(TreeDecomposition&, TreeDecomposition::vertex_descriptor);
	uint64_t calculateCutWeight(TreeDecomposition& td, TreeDecomposition::vertex_descriptor root);
	void retraceCut(TreeDecomposition&, TreeDecomposition::vertex_descriptor, uint64_t);

	atomicQueue leafs;
	std::vector<std::atomic_flag> joinWasVisited;
};

