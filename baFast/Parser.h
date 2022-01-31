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
	std::queue<size_t> HOLLOWsplit(TreeDecomposition&, TreeDecomposition::vertex_descriptor, int);
	void traverseUpThread(TreeDecomposition&, TreeDecomposition::vertex_descriptor);
	uint64_t binomial(int, int);
	std::string tableOfNode(TreeDecomposition&, TreeDecomposition::vertex_descriptor);
	uint64_t getIndexOfSubset(int*, int);
	void perm2Indices(int*, uint64_t, int);
	uint32_t computeWeightIntroduceContained(TreeDecomposition&, size_t, const std::vector<size_t>&, int*, int);
	uint32_t computeWeightIntroduceNotContained(TreeDecomposition&, size_t, const std::vector<size_t>&, int*, int);
	uint32_t computeWeightJoin(TreeDecomposition&, const std::vector<size_t>, uint64_t);
	int calculateCutWeight(TreeDecomposition& td, TreeDecomposition::vertex_descriptor root);
	uint64_t mirror(uint64_t totalLength, uint64_t index);
	
	atomicQueue leafs;
	std::vector<std::atomic_flag> joinWasVisited;
};

