#include "Parser.h"
#include <thread>
#include <mutex>
#include <float.h>

#define THREADCOUNT 5

using namespace std;

TreeDecomposition* Parser::parse(istream& in) {

	string line;

	//read first line
	getline(in, line);
	istringstream iss(line);
	string isSolution, isTreeDecomp;
	int bagsNumber, maxBagSize, verticesNumber;
	iss >> isSolution >> isTreeDecomp >> bagsNumber >> maxBagSize >> verticesNumber;

	//create TreeDecomposition with bags starting at index 1
	TreeDecomposition* td;
	if (isSolution[0] != 's' || (isTreeDecomp.compare("td") != 0)) {
		throw std::exception("WARNING: error in input, returning unfinished td!");
	}
	//no error so far, proceed to create the actual object;
	td = new TreeDecomposition(bagsNumber + 1);
	(*td)[graph_bundle].tablesInitialized = false;
	(*td)[graph_bundle].maxBagSize = maxBagSize;
	(*td)[graph_bundle].numberBags = bagsNumber;
	(*td)[graph_bundle].numberVerticesOfDecomposedGraph = verticesNumber;

	//read bags
	for (int i = 1; i <= bagsNumber; i++) {
		getline(in, line);
		string b;
		int bagID, vertexID;
		istringstream lineStream(line);
		lineStream >> b >> bagID;
		if (b[0] != 'b' || i != bagID) {
			throw std::exception("WARNING: potential mismatch when parsing decomposition while trying to read bag contents");
			cout << "Debug Info: " << b.compare("b") << " " << i << " " << bagID << " line where error: " << line << endl;
		}
		while (lineStream >> vertexID)
			(*td)[bagID].bag.push_back(vertexID);
	}

	//read structure
	while (getline(in, line)) {
		if (line[0] == 'p') {

			fillAdjacencyMatrix(*td, in);
			break;
		}
		else {
			int parentNode, childNode;
			istringstream lineStream(line);
			if (!(lineStream >> parentNode >> childNode)) {
				throw std::exception("WARNING: some error occured while parsing edges between nodes! resulting treedecomp is incomplete!");
			}
			add_edge(vertex(parentNode, *td), vertex(childNode, *td), *td);
		}
	}

	//ensure all nodes are unvisited
	
	(*td)[graph_bundle].root = calculateOptimalRoot(*td);
	//(*td)[graph_bundle].root = 1;
	makeNice(*td);
	
	return td;
}
void Parser::fillAdjacencyMatrix(TreeDecomposition& td, istream& in, int width) {
	string line;
	if (width == 0) {
		getline(in, line);
		istringstream lineStream(line);
		string trash;
		if (!(lineStream >> trash >> trash >> width)) {
			cout << line << endl;
			throw std::exception("WARNING: some error occurred while parsing the width of the original graph!");
		}
	}
	//setup matrix to be filled, set it to 0
	int** matrix = new int* [width];
	for (int i = 0; i < width; i++) {
		matrix[i] = new int[width];
		for (int k = 0; k < width; k++)
			matrix[i][k] = 0;
	}
	cout << "setup adjacency matrix of width " << width << endl;
	while (getline(in, line)) {
		istringstream lineStream(line);
		int a, b;
		if (!(lineStream >> a >> b))
			throw std::exception("WARNING: some error occurred while parsing an edge of original graph!");
		matrix[a - 1][b - 1] = 1;
		matrix[b - 1][a - 1] = 1;
	}
	td[graph_bundle].adjacencyMatrix = matrix;
}

void Parser::print(TreeDecomposition& td, TreeDecomposition::vertex_descriptor node) {
	for (int i = 0; i < td[node].height; i++) {
		cout << "   ";
	}
	cout << "<" << node << "> [";
	for (int i = 0; i < td[node].bag.size(); i++) {
		if (i < td[node].bag.size() - 1) {
			cout << td[node].bag.at(i) << ", ";
		}
		else {
			cout << td[node].bag.at(i);
		}
	}
	string type;
	switch (td[node].type) {
	case 1: type = "INT " + std::to_string(td[node].specialVertex);
		break;
	case 2: type = "FRGT " + std::to_string(td[node].specialVertex);
		break;
	case 3: type = "JN";
		break;
	case 4: type = "LF";
		break;
	}
	cout << "] <" << "" << "> " << type << " values: " << tableOfNode(td, node) << endl;
	for (auto it = td[node].children.begin(); it != td[node].children.end(); it++) {
		print(td, *it);
	}
}
std::string Parser::tableOfNode(TreeDecomposition& td, TreeDecomposition::vertex_descriptor v) {
	if (!td[graph_bundle].tablesInitialized)
		return "[uninitialized]";
	string out;


	const int bagSize = td[v].bag.size();
	const int Y = td[v].inducedSubgraphSize;
	const int degOfFreedom = Y - bagSize + 1;
	uint64_t totalLength = ((uint64_t)1 << (bagSize - 1)) * degOfFreedom;

	out.append("[");
	for (int i = 0; i < totalLength; i++) {
		if (i == totalLength - 1) {
			out.append(to_string(td[v].values[i]));
		}
		else {
			out.append(to_string(td[v].values[i]) + ", ");
		}
	}
	out.append("]");

	if (td[v].type == 2) {
		out.append(" BITVECTOR: [");
		for (int i = 0; i < totalLength; i++) {
			if (i == totalLength - 1) {
				out.append(to_string((*td[v].forget_bitset)[i]));
			}
			else {
				out.append(to_string((*td[v].forget_bitset)[i]) + ", ");
			}
		}
		out.append("]");
	}
	else if (td[v].type == 3) {
		out.append(" JVECTOR: [");
		for (int i = 0; i < totalLength; i++) {
			if (i == totalLength - 1) {
				out.append(to_string(td[v].join_j[i]));
			}
			else {
				out.append(to_string(td[v].join_j[i]) + ", ");
			}
		}
		out.append("]");
	}

	return out;
}
void Parser::exportDimax(TreeDecomposition& td, std::ostream& out) {
	out << "s td " << td.vertex_set().size() - 1 << " " << td[graph_bundle].maxBagSize << " " << td[graph_bundle].numberVerticesOfDecomposedGraph << "\n";
	for (size_t i = 1; i < td.vertex_set().size(); i++) {
		out << "b " << i;
		for (int k = 0; k < td[i].bag.size(); k++) {
			out << " " << td[i].bag[k];
		}
		out << "\n";
	}
	for (size_t i = 1; i < td.vertex_set().size(); i++) {
		for (int k = 0; k < td[i].children.size(); k++) {
			if (i < td[i].children[k]) {
				out << i << " " << td[i].children[k] << "\n";
			}
			else {
				out << td[i].children[k] << " " << i << "\n";
			}
		}
	}
}


void Parser::traverseUpThread(TreeDecomposition& td, TreeDecomposition::vertex_descriptor root) {
	while (TreeDecomposition::vertex_descriptor node = leafs.pop()) {
		//cout << "thread with id " << std::this_thread::get_id() << " has obtained leaf " << node << endl;

		while (true) {
			if (td[node].type == 3) {//check if join
				if (!joinWasVisited[node].test_and_set()) {//if 0th child(always use 0th child) was not checked before, the other thread has not checked it yet, meaning it will then receive a true -> false such that it will not break and keep on going
					//cout << "thread with id " << std::this_thread::get_id() << " has now finished its chain at join node " << node << " because its sibling has not finished computing yet" << endl;
					break;
				}
				else {
					td[node].inducedSubgraphSize = td[td[node].children[0]].inducedSubgraphSize + td[td[node].children[1]].inducedSubgraphSize - td[node].bag.size();
				}
			} 
			else if (td[node].type == 1) {
				td[node].inducedSubgraphSize = td[td[node].children[0]].inducedSubgraphSize + 1;
			}
			else if (td[node].type == 2) {
				td[node].inducedSubgraphSize = td[td[node].children[0]].inducedSubgraphSize;
			}
			else {
				td[node].inducedSubgraphSize = 1;
			}

			const size_t bagSize = td[node].bag.size();
			const uint64_t Y = td[node].inducedSubgraphSize;
			const int degOfFreedom = Y - bagSize + 1;
			uint64_t totalLength = ((uint64_t)1 << (bagSize - 1)) * degOfFreedom;
			td[node].values = new uint32_t [totalLength];//want to be able to choose full bagSize as well
			uint32_t* vals = td[node].values;

			//cout << "thread with id " << std::this_thread::get_id() << " is calculating node " << node << " of type " << td[node].type << " at height " << td[node].height << " with " << ((uint64_t)1 << bagSize) << " * " << degOfFreedom << " values needing to be computed" << endl;


			if (td[node].type == 1) {//introduce
				uint64_t childTotalLength = totalLength >> 1;
				size_t introducedVertex = td[node].specialVertex;
				size_t child = td[node].children[0];
				int introducedVertexPos = 0;
				while (introducedVertexPos < bagSize) {
					if (td[node].bag[introducedVertexPos] == introducedVertex)
						break;
					introducedVertexPos++;
				}

				uint64_t counter = 0;
				uint64_t introducedVertexNotContainedCounter = 0;
				uint64_t introducedVertexContainedCounter = 0;
				int* indices = new int[bagSize];
				for (int k = 0; k < bagSize + 1; k++) {//for every k-class
					uint64_t perm = ((uint64_t)1 << k) - 1;
					uint64_t maxSubsets = binomial(bagSize, k);
					for (int s = 1; s <= maxSubsets; s++) {
						perm2Indices(indices, perm, k);
						if (perm & (((uint64_t)1) << introducedVertexPos)) {// x contained in S', never triggers when k = 0
							int cutWeight = computeWeightIntroduceContained(td, introducedVertex, td[node].bag, indices, k);
							for (int i = 0; i < degOfFreedom; i++) {
								vals[counter++] = td[child].values[introducedVertexContainedCounter++] + cutWeight;
								
								//cout << "node " << node << " executed CONTAINED vals[" << counter - 1 << "] = td[child].values[" << introducedVertexContainedCounter - 1 << "] + " << cutWeight << endl;
							}
						}
						else {//x not contained in S'
							int cutWeight = computeWeightIntroduceNotContained(td, introducedVertex, td[node].bag, indices, k);
							for (int i = 0; i < degOfFreedom; i++) {
								uint64_t uncontainedIndex = introducedVertexNotContainedCounter; //mirror
								if (introducedVertexNotContainedCounter >= childTotalLength) {
									uncontainedIndex = (totalLength) - introducedVertexNotContainedCounter++ - 1; //totalLength is childTotalLength * 2 already!
								}
								else {
									uncontainedIndex = introducedVertexNotContainedCounter++;
								}
								vals[counter++] = td[child].values[uncontainedIndex] + cutWeight;
								//cout << "node " << node << " executed UNCONTAINED vals[" << counter - 1 << "] = td[child].values[" << uncontainedIndex << "] + " << cutWeight << endl;
							}
						}

						//check that now out of bounds, necessary because often times not evenly split at some integer k, only need to check for every s, because those are always
						if (counter >= totalLength)
							goto finishedCalc;
						//calculate next perm
						uint64_t t = perm | (perm - 1);
						unsigned long temp;
						_BitScanForward64(&temp, perm); //compiler directive, temp is loaded with amount of trailing zeros of perm.
						perm = (t + 1) | (((~t & -~t) - 1) >> (temp + 1));
					}
				}
				delete[] indices;
			}
			else if (td[node].type == 2) {
				uint64_t childTotalLength = ((uint64_t)1 << (bagSize)) * (degOfFreedom - 1);
				size_t forgottenVertex = td[node].specialVertex;
				size_t child = td[node].children[0];

				td[node].forget_bitset = new boost::dynamic_bitset<uint64_t>(totalLength);

				int forgottenVertexPos = 0;
				while (forgottenVertexPos < bagSize) {
					if (forgottenVertex < td[node].bag[forgottenVertexPos])
						break;
					forgottenVertexPos++;
				}
				int counter = 0;
				int* indices = new int[bagSize + 1]; //need one more because we need to be able to store child's indices!
				uint64_t childHowManyBeforeKClassLate = 0;
				uint64_t childHowManyBeforeKClass = 0;
				for (int k = 0; k < bagSize + 1; k++) {//for every k-class

					uint64_t perm = ((uint64_t)1 << k) - 1;
					childHowManyBeforeKClass += binomial(bagSize + 1, k) * (degOfFreedom - 1); //important that this is up here, essentially making it one earlier

					uint64_t maxSubsets = binomial(bagSize, k);
					for (int s = 1; s <= maxSubsets; s++) {
						uint64_t leftPart = (UINT64_MAX << forgottenVertexPos) & perm;
						leftPart = leftPart << 1;
						uint64_t rightPart = (((uint64_t)1 << forgottenVertexPos) - 1) & perm;
						uint64_t permChild = leftPart | rightPart; //inserted 0 at forgottenVertexPos
						perm2Indices(indices, permChild, k);
						uint64_t baseIndexLeft = (getIndexOfSubset(indices, k) - 1) * (degOfFreedom - 1);
						permChild = permChild | ((uint64_t)1 << forgottenVertexPos); //set bit at forgottenVertexPos
						perm2Indices(indices, permChild, k + 1);
						uint64_t baseIndexRight = (getIndexOfSubset(indices, k + 1) - 1) * (degOfFreedom - 1);
						for (int i = 0; i < degOfFreedom; i++) {
							if (i == 0) {
								vals[counter++] = td[child].values[childHowManyBeforeKClassLate + baseIndexLeft];
								//cout << "node " << node << " calculated LEFT td[child].values[" << childHowManyBeforeKClassLate << " + " << baseIndexLeft << " + 0]" << endl;
							}
							else if (i == degOfFreedom - 1) {
								uint64_t rightIndex = childHowManyBeforeKClass + baseIndexRight + i - 1; //mirror
								if (rightIndex >= childTotalLength) {
									rightIndex = (childTotalLength << 1) - rightIndex - 1;
								}
								(*td[node].forget_bitset)[counter] = true;
								vals[counter++] = td[child].values[rightIndex];
								//cout << "node " << node << " calculated RIGHT------ td[child].values[" << rightIndex << "]" << endl;
							}
							else {
								uint32_t left = td[child].values[childHowManyBeforeKClassLate + baseIndexLeft + i];
							
								uint64_t rightIndex = childHowManyBeforeKClass + baseIndexRight + i - 1; //mirror
								if (rightIndex >= childTotalLength) {
									rightIndex = (childTotalLength << 1) - rightIndex - 1;
								}
								uint32_t right = td[child].values[rightIndex];
								if (right < left) {
									(*td[node].forget_bitset)[counter] = true;
									vals[counter++] = right;
								}
								else {
									vals[counter++] = left;
								}
							}
						}

						if (counter >= totalLength)
							goto finishedCalc;
						//calculate next perm
						uint64_t t = perm | (perm - 1);
						unsigned long temp;
						_BitScanForward64(&temp, perm); //compiler directive, temp is loaded with amount of trailing zeros of perm.
						perm = (t + 1) | (((~t & -~t) - 1) >> (temp + 1));
					}
					childHowManyBeforeKClassLate += binomial(bagSize + 1, k) * (degOfFreedom - 1);
				}
				delete[] indices;
			}
			else if (td[node].type == 3) {
				size_t leftChild = td[node].children[0];
				size_t rightChild = td[node].children[1];
				int degFreedomLeft = td[leftChild].inducedSubgraphSize - bagSize + 1;
				int degFreedomRight = td[rightChild].inducedSubgraphSize - bagSize + 1;
				uint64_t totalLengthLeft = ((uint64_t)1 << (bagSize - 1)) * degFreedomLeft;
				uint64_t totalLengthRight = ((uint64_t)1 << (bagSize - 1)) * degFreedomRight;
				uint64_t counter = 0;
				int* indices = new int[bagSize];
				
				td[node].join_j = new uint16_t[totalLength];

				uint64_t childHowManyBeforeKClassLeft = 0;
				uint64_t childHowManyBeforeKClassRight = 0;
				for (int k = 0; k < bagSize + 1; k++) {//for every k-class

					uint64_t perm = ((uint64_t)1 << k) - 1;
					
					uint64_t maxSubsets = binomial(bagSize, k);
					for (int s = 1; s <= maxSubsets; s++) {
						int cutWeight = computeWeightJoin(td, td[node].bag, perm);
						perm2Indices(indices, perm, k);
						uint64_t baseIndex = getIndexOfSubset(indices, k);
						uint64_t leftBaseIndex = (baseIndex - 1) * degFreedomLeft;
						uint64_t rightBaseIndex = (baseIndex - 1) * degFreedomRight;
						for (int i = 0; i < degOfFreedom; i++) {
							uint64_t min = UINT64_MAX;
							int j_lowerBound = i - degFreedomRight + 1;
							if (j_lowerBound < 0)
								j_lowerBound = 0;
							int j_upperBound = degFreedomLeft - 1;
							if (j_upperBound > i)
								j_upperBound = i;
							for (int j = j_lowerBound; j <= j_upperBound; j++) {
								int notJ = i - j;
								if (j < degFreedomLeft && notJ < degFreedomRight) {
									uint64_t candidate = (uint64_t)td[leftChild].values[childHowManyBeforeKClassLeft + leftBaseIndex + j] + (uint64_t)td[rightChild].values[childHowManyBeforeKClassRight + rightBaseIndex + notJ];
									//cout << "node " << node << " calculated td[leftChild].values[" << childHowManyBeforeKClassLeft << " + " << leftBaseIndex << " + " << j << "] + td[rightChild].values[" << childHowManyBeforeKClassRight << " + " << rightBaseIndex << " + " << notJ << "] which is candidate value " << candidate << " = " << (uint64_t)td[leftChild].values[childHowManyBeforeKClassLeft + leftBaseIndex + j] << " + " << (uint64_t)td[rightChild].values[childHowManyBeforeKClassRight + rightBaseIndex + notJ] << endl;
									if (candidate < min) {
										min = candidate;
										td[node].join_j[counter] = j;
									}
										
								}
							}
							vals[counter++] = min - cutWeight;
							//cout << "node " << node << " filled in minimum of " << min << " into td[" << node << "].values[" << counter - 1 << "]" << endl;
						}

						if (counter >= totalLength)
							goto finishedCalc;

						uint64_t t = perm | (perm - 1);
						unsigned long temp;
						_BitScanForward64(&temp, perm); //compiler directive, temp is loaded with amount of trailing zeros of perm.
						perm = (t + 1) | (((~t & -~t) - 1) >> (temp + 1));
					}
					uint64_t binom = binomial(bagSize, k);
					childHowManyBeforeKClassLeft += binom * (degFreedomLeft);
					childHowManyBeforeKClassRight += binom * (degFreedomRight);
				}
				delete[] indices;

			}
			else {//leaf
				vals[0] = 0; 
			}
finishedCalc:
			for (auto it = td[node].children.begin(); it != td[node].children.end(); it++) {
				delete[] td[*it].values;
			}
			if (node == root) {
				//cout << "thread with id " << std::this_thread::get_id() << " has reached root and will now cease activity!" << endl;
				retraceCut(td, root, calculateCutWeight(td, node));
				//retraceCut(td, 11, 22);
				return;
			}
			node = td[node].parent;
		}
	}
	//cout << "thread with id " << std::this_thread::get_id() << " has run out of leafs and will now cease activity!" << endl;
}
void Parser::debugAlgorithm(TreeDecomposition& td) {
	//expects fully prepared nice tree decomposition to work on. Then 
	TreeDecomposition::vertex_descriptor root = td[graph_bundle].root;
	vector<atomic_flag> a(num_vertices(td));
	joinWasVisited = std::move(a);

	thread myWorkers[THREADCOUNT];
	for (int i = 0; i < THREADCOUNT; i++) {
		myWorkers[i] = thread(&Parser::traverseUpThread, this, std::ref(td), root);
	}

	for (int i = 0; i < THREADCOUNT; i++) {
		myWorkers[i].join();
	}
	td[graph_bundle].tablesInitialized = true;
}



void Parser::makeNice(TreeDecomposition& td) { //ASSUMPTION: THERE ARE NO TWO EQUAL BAGS CONNECTED TO EACH OTHER, EXCEPT IF ITS A JOIN NODE THAT IS ALREADY NICE
	//make list of properties that need to be updated, such as height, leftChild, rightChild then in one central loop, first enforce the max 2 children with equal bag rule, then
	//subsequently the other rules such as forgetting and introducing instead of swapping. then add attachments to leafs, essentially strings of introduce nodes, already written in other code.
	// node is swap node if they have one child and the same bagsize, otherwise they are introduce nodes already(but check that bagsize grows by 1 at most!!). If a node has smaller bagsize (-1!!)than its child it must be a forget node already.
	if (td[graph_bundle].isNice) {
		return;
	}
	//ensure all unvisited
	for (auto it = vertices(td).first + 1; it != vertices(td).second; it++) {
		(td)[*it].visited = false;
	}

	TreeDecomposition::vertex_descriptor root = td[graph_bundle].root;
	std::queue<TreeDecomposition::vertex_descriptor> q;
	td[root].height = 0;
	td[root].visited = true;
	q.push(root);
	while (!q.empty()) {
		TreeDecomposition::vertex_descriptor v = q.front();
		for (auto it = out_edges(v, td).first; it != out_edges(v, td).second; it++) {
			TreeDecomposition::vertex_descriptor w = boost::target(*it, td);
			if (!td[w].visited) {
				td[w].visited = true;
				td[w].parent = v;
				td[v].children.push_back(w);
				q.push(w);
			}
		}
		if (v != root)
			td[v].height = td[td[v].parent].height + 1;
		//here do any manipulation on v that might be necessary to make it nice. Original children stored in td[v].children!
		size_t childrenAmount = td[v].children.size();
		if (childrenAmount == 0) {//attach introduce nodes to leaf
			vector<TreeDecomposition::vertex_descriptor> nextBag = td[v].bag;
			//handles artificial case of empty leaf in tree decomposition
			if (nextBag.size() > 1) {
				td[v].specialVertex = nextBag.back();
				nextBag.pop_back();
				td[v].type = INTRODUCE;
			}
			else {
				td[v].type = LEAF;
				leafs.push(v);
			}
			TreeDecomposition::vertex_descriptor currentParent = v;
			for (int i = 1; i < td[v].bag.size(); i++) {
				TreeDecomposition::vertex_descriptor newNode = add_vertex(td);
				if (i == td[v].bag.size() - 1) {
					td[newNode].type = LEAF;
					leafs.push(newNode);
				}
				else {
					td[newNode].type = INTRODUCE;
				}
				td[newNode].parent = currentParent;
				td[currentParent].children.push_back(newNode);
				td[newNode].bag = nextBag;
				td[newNode].specialVertex = nextBag.back();
				td[newNode].height = td[currentParent].height + 1;
				nextBag.pop_back();
				currentParent = newNode;
			}
		}
		else if (childrenAmount == 1) {//maxIntroduce = 1, maxForgotten = 1, attach chain between v and its children w
			insertChainBetweenUnequalNodes(td, v, td[v].children[0]);
		}
		else {// 1 split adds 1 additional exposedEndpoint
			if (childrenAmount == 2 && td[td[v].children[0]].bag == td[v].bag && td[td[v].children[1]].bag == td[v].bag) {//is nice join already
				td[v].type = JOIN;
			}
			else {
				vector<size_t> children = td[v].children;
				std::queue<size_t> exposedEndpoints = split(td, v, childrenAmount - 1);
				while (!exposedEndpoints.empty()) {
					TreeDecomposition::vertex_descriptor endPoint = exposedEndpoints.front();
					insertChainBetweenUnequalNodes(td, endPoint, children[children.size() - 1]);
					exposedEndpoints.pop();
					children.pop_back();
				}
			}
		}
		q.pop();
	}
	td[graph_bundle].isNice = true;
}
std::queue<size_t> Parser::split(TreeDecomposition& td, TreeDecomposition::vertex_descriptor v, int count) {
	std::queue<TreeDecomposition::vertex_descriptor> q;
	td[v].children.clear();
	q.push(v);
	for (int i = 0; i < count; i++) {
		TreeDecomposition::vertex_descriptor cloneMe = q.front();
		TreeDecomposition::vertex_descriptor newClone1 = add_vertex(td);
		TreeDecomposition::vertex_descriptor newClone2 = add_vertex(td);
		q.push(newClone1);
		q.push(newClone2);
		td[newClone1].bag = td[cloneMe].bag;
		td[newClone1].bagSize = td[cloneMe].bagSize;
		td[newClone1].height = td[cloneMe].height + 1;
		td[newClone1].parent = cloneMe;
		td[newClone2].bag = td[cloneMe].bag;
		td[newClone2].bagSize = td[cloneMe].bagSize;
		td[newClone2].height = td[cloneMe].height + 1;
		td[newClone2].parent = cloneMe;
		td[cloneMe].type = JOIN;

		td[cloneMe].children.push_back(newClone1);
		td[cloneMe].children.push_back(newClone2);
		q.pop();
	}
	return q;
}
void Parser::insertChainBetweenUnequalNodes(TreeDecomposition& td, TreeDecomposition::vertex_descriptor v, TreeDecomposition::vertex_descriptor w) {
	vector<TreeDecomposition::vertex_descriptor> introduceMe;
	vector<TreeDecomposition::vertex_descriptor> forgetMe;
	for (int i = 0; i < td[v].bag.size(); i++) {
		//if some vertex is not in w's bag, add it to introduceMe
		bool isContained = false;
		size_t check = td[v].bag[i];
		for (int k = 0; k < td[w].bag.size(); k++) {
			if (td[w].bag[k] > check) //make use of ordering of bags to end search early
				break;
			if (td[w].bag[k] == check) {
				isContained = true;
				break;
			}
		}
		if (!isContained)
			introduceMe.push_back(check);
	}


	for (int i = 0; i < td[w].bag.size(); i++) {
		//if some vertex is not in v's bag, add it to forgetMe
		bool isContained = false;
		size_t check = td[w].bag[i];
		for (int k = 0; k < td[v].bag.size(); k++) {
			if (td[v].bag[k] > check) //make use of ordering of bags to end search early
				break;
			if (td[v].bag[k] == check) {
				isContained = true;
				break;
			}
		}
		if (!isContained)
			forgetMe.push_back(check);
	}

	//handle special case that there is nothing to introduce, making v a forget node, and we need to add one less forget node to link them!
	vector<size_t> curBagIntroduce;
	int skipLast = 0;
	if (!introduceMe.empty()) {
		td[v].type = INTRODUCE;
		td[v].specialVertex = introduceMe[0];
		curBagIntroduce = td[v].bag;
		removeElementFromBag(curBagIntroduce, introduceMe[0]);
	}
	else {
		td[v].type = FORGET;
		td[v].specialVertex = forgetMe[0];
		skipLast = 1;
	}

	//add forget nodes from bottom up
	td[w].height = td[v].height + 1 + forgetMe.size() + introduceMe.size() - 1;
	vector<size_t> curBag = td[w].bag;
	TreeDecomposition::vertex_descriptor child = w;
	for (int i = 0; i < forgetMe.size() - skipLast; i++) {
		TreeDecomposition::vertex_descriptor newForgetNode = add_vertex(td);
		td[newForgetNode].type = FORGET;
		td[newForgetNode].bag = removeElementFromBag(curBag, forgetMe[i]); //curBag gets modified here!
		td[newForgetNode].specialVertex = forgetMe[i];
		td[newForgetNode].children.push_back(child);
		td[newForgetNode].height = td[v].height + introduceMe.size() - 1 + forgetMe.size() - i;
		td[child].parent = newForgetNode;
		child = newForgetNode;
	}


	td[v].children.clear(); //need to delete old children!!

	TreeDecomposition::vertex_descriptor curParent = v;
	for (int i = 1; i < introduceMe.size(); i++) { //highest introduce node already exists, is v!
		TreeDecomposition::vertex_descriptor newIntroduceNode = add_vertex(td);
		td[newIntroduceNode].type = INTRODUCE;
		td[newIntroduceNode].bag = curBagIntroduce;
		removeElementFromBag(curBagIntroduce, introduceMe[i]); //curBag gets modified here!
		td[newIntroduceNode].specialVertex = introduceMe[i];
		td[newIntroduceNode].parent = curParent;
		td[newIntroduceNode].height = td[curParent].height + 1;
		td[curParent].children.push_back(newIntroduceNode);
		curParent = newIntroduceNode;
	}
	//link two chains together
	td[curParent].children.push_back(child);
	td[child].parent = curParent;
}
vector<size_t>& Parser::removeElementFromBag(vector<size_t>& original, size_t elem) {
	for (auto it = original.begin(); it != original.end(); it++) {
		if (*it == elem) {
			original.erase(it);
			break;
		}
	}
	return original;
}
uint64_t Parser::binomial(int n, int k) {
	if (k < 0 || k > n)
		return 0;
	uint64_t t = 1;
	if (k < n - k) {
		for (int i = n; i > n - k; i--) {
			t *= i;
			t /= n - i + 1;
		}
	}
	else {
		for (int i = n; i > k; i--) {
			t *= i;
			t /= n - i + 1;
		}
	}
	return t;
}
uint64_t Parser::getIndexOfSubset(int* indices, int len) {
	uint64_t s = 1;
	for (int index = len - 1; index >= 0; index--) {
		//indices[index] accesses current index and index is also #leftToConsider!
		if (int movedBy = indices[index] - index) {//has moved!!
			s += binomial(index + movedBy, index + 1);
		}
	}
	return s;
}

void Parser::perm2Indices(int* indices, uint64_t perm, int k) {
	if (k == 0) {
		return;
	}
	int index = 0;
	for (int i = 0; i < 64; i++) {
		if (perm & 1) {
			indices[index++] = i;
			if (index == k)
				return;
		}
		perm = perm >> 1;
	}
}


uint64_t Parser::calculateCutWeight(TreeDecomposition& td, TreeDecomposition::vertex_descriptor root) {
	const int bagSize = td[root].bag.size();
	const int Y = td[root].inducedSubgraphSize;
	const int degOfFreedom = Y - bagSize + 1;
	uint64_t totalLength = ((uint64_t)1 << (bagSize - 1)) * degOfFreedom;
	uint32_t* vals = td[root].values;
	int minWeight = 0;
	uint64_t counter = 0;
	uint64_t counterOfMin = 0;
	double min = DBL_MAX;
	for (int k = 0; k < bagSize + 1; k++) {
		uint64_t maxSubsets = binomial(bagSize, k);
		for (int s = 1; s <= maxSubsets; s++) {
			for (int i = k; i < k + degOfFreedom; i++) {
				int notI = Y - i;
				if (!(notI && i)) {
					counter++;
					continue;
				}
				double candidate = (double)vals[counter++] / (i * (notI));
				//cout << candidate << endl;
				if (candidate < min) {
					min = candidate;
					minWeight = vals[counter - 1];
					counterOfMin = counter - 1;
				}
			}
			if (counter >= totalLength)
				goto finish;
		}
	}
finish:
	cout << "minWeight: " << minWeight << " minSparsestCutWeight: " << min << " selected index in top table: " << counterOfMin << endl;

	return counterOfMin;
}
uint32_t Parser::computeWeightIntroduceContained(TreeDecomposition& td, size_t introducedVertex, const vector<size_t>& bag, int* indices, int lengthOfSDash) {
	uint32_t weight = 0;

	size_t* set_difference = new size_t[bag.size()];
	vector<size_t> sDash;
	for (int i = 0; i < lengthOfSDash; i++) {
		sDash.push_back(bag[indices[i]]);
	}

	auto lastIndex = std::set_difference(bag.begin(), bag.end(), sDash.begin(), sDash.end(), set_difference);

	size_t row = introducedVertex - 1;
	for (auto it = set_difference; it != lastIndex; it++) {
		if (*it == introducedVertex)
			continue;
		weight += td[graph_bundle].adjacencyMatrix[row][*it - 1];
	}
	delete[] set_difference;
	return weight;
}
uint32_t Parser::computeWeightIntroduceNotContained(TreeDecomposition& td, size_t introducedVertex, const vector<size_t>& bag, int* indices, int amountNeighbors) {
	uint32_t weight = 0;
	size_t row = introducedVertex - 1;
	for (int i = 0; i < amountNeighbors; i++) {
		weight += td[graph_bundle].adjacencyMatrix[row][bag[indices[i]] - 1];
	}
	return weight;
}
uint32_t Parser::computeWeightJoin(TreeDecomposition& td, vector<size_t> copiedBag, uint64_t perm) {
	int weight = 0;
	uint64_t mask = 1;
	vector<size_t> sDash;
	const int leftOverSize = copiedBag.size();
	for (int i = 0; i < leftOverSize; i++) {
		if (mask & perm) {
			sDash.push_back(copiedBag[i]);
			copiedBag[i] = 0;
		}
		mask = mask << 1;
	}
	const int sDashSize = sDash.size();
	for (int i = 0; i < sDashSize; i++) {
		for (int j = 0; j < leftOverSize; j++) {
			if (copiedBag[j] == 0)
				continue;
			weight += td[graph_bundle].adjacencyMatrix[sDash[i] - 1][copiedBag[j] - 1];
		}
	}
	return weight;
}



size_t Parser::calculateOptimalRoot(TreeDecomposition& td) {
	uint64_t weight;
	uint64_t maxWeight = 0;
	uint64_t minWeight = UINT64_MAX;
	size_t worstRoot = 1;
	size_t bestRoot = 1;
	std::queue<TreeDecomposition::vertex_descriptor> q;
	for (auto it = vertices(td).first + 1; it != vertices(td).second; it++) {
		(td)[*it].visited = false;
	}
	
	for (auto rootCandidate = vertices(td).first + 1; rootCandidate != vertices(td).second; rootCandidate++) {
		weight = 0;
		TreeDecomposition c_td;
		boost::copy_graph(td, c_td);

		std::queue<size_t> leafsLocal = HOLLOWmakeNice(c_td, *rootCandidate);
		//print(*c_td, *rootCandidate);
		uint64_t weight = calculateNumberOfOperations(c_td, *rootCandidate, leafsLocal);
		//delete c_td;
		//cout << "root " << *rootCandidate << " has weight " << weight << endl;
		if (weight < minWeight) {
			bestRoot = *rootCandidate;
			minWeight = weight;
		}
		if (weight > maxWeight) {
			worstRoot = *rootCandidate;
			maxWeight = weight;
		}
	}
	cout << "Best root was node " << bestRoot << " with a weight of " << minWeight << endl;
	cout << "Worst root was node " << worstRoot << " with a weight of " << maxWeight << endl;
	return bestRoot;
}

std::queue<size_t> Parser::HOLLOWmakeNice(TreeDecomposition& td, size_t root){
	std::queue<size_t> leafsLocal;
	for (auto it = vertices(td).first + 1; it != vertices(td).second; it++) {
		td[*it].bagSize = td[*it].bag.size();
	}
	std::queue<TreeDecomposition::vertex_descriptor> q;
	td[root].height = 0;
	td[root].visited = true;
	
	q.push(root);
	while (!q.empty()) {
		TreeDecomposition::vertex_descriptor v = q.front();
		for (auto it = out_edges(v, td).first; it != out_edges(v, td).second; it++) {
			TreeDecomposition::vertex_descriptor w = boost::target(*it, td);
			if (!td[w].visited) {
				td[w].visited = true;
				td[w].parent = v;
				td[v].children.push_back(w);
				q.push(w);
			}
		}
		if (v != root) {
			td[v].height = td[td[v].parent].height + 1;
		}
		//here do any manipulation on v that might be necessary to make it nice. Original children stored in td[v].children!
		size_t childrenAmount = td[v].children.size();
		if (childrenAmount == 0) {//attach introduce nodes to leaf
			//vector<TreeDecomposition::vertex_descriptor> nextBag = td[v].bag;
			//handles artificial case of empty leaf in tree decomposition
			if (td[v].bagSize > 1) {
				//td[v].specialVertex = nextBag.back();
				//nextBag.pop_back();
				td[v].type = INTRODUCE;
			}
			else {
				td[v].type = LEAF;
				leafsLocal.push(v);
			}
			TreeDecomposition::vertex_descriptor currentParent = v;
			for (int i = 1; i < td[v].bagSize; i++) {
				TreeDecomposition::vertex_descriptor newNode = add_vertex(td);
				if (i == td[v].bagSize - 1) {
					td[newNode].type = LEAF;
					leafsLocal.push(newNode);
				}
				else {
					td[newNode].type = INTRODUCE;
				}
				td[newNode].parent = currentParent;
				td[currentParent].children.push_back(newNode);
				td[newNode].bagSize = td[currentParent].bagSize - 1;
				//td[newNode].specialVertex = nextBag.back();
				td[newNode].height = td[currentParent].height + 1;
				//nextBag.pop_back();
				currentParent = newNode;
			}
		}
		else if (childrenAmount == 1) {//maxIntroduce = 1, maxForgotten = 1, attach chain between v and its children w
			HOLLOWinsertChainBetweenUnequalNodes(td, v, td[v].children[0]);
		}
		else {// 1 split adds 1 additional exposedEndpoint
			if (childrenAmount == 2 && td[td[v].children[0]].bag == td[v].bag && td[td[v].children[1]].bag == td[v].bag) {//is nice join already
				td[v].type = JOIN;
			}
			else {
				vector<size_t> children = td[v].children;
				std::queue<size_t> exposedEndpoints = split(td, v, childrenAmount - 1);
				while (!exposedEndpoints.empty()) {
					TreeDecomposition::vertex_descriptor endPoint = exposedEndpoints.front();
					HOLLOWinsertChainBetweenUnequalNodes(td, endPoint, children[children.size() - 1]);
					exposedEndpoints.pop();
					children.pop_back();
				}
			}
		}
		q.pop();
	}
	return leafsLocal;
}
void Parser::HOLLOWinsertChainBetweenUnequalNodes(TreeDecomposition& td, TreeDecomposition::vertex_descriptor v, TreeDecomposition::vertex_descriptor w) {
	/*if (td[v].bag == td[w].bag) {
		throw std::exception("WARNING: HOLLOWinsertChainBetweenUnequalNodes was called for two nodes with equal bags!");
	}*/
	uint64_t introduceMeSize = 0;
	uint64_t forgetMeSize = 0;
	for (int i = 0; i < td[v].bag.size(); i++) {
		//if some vertex is not in w's bag, add it to introduceMe
		bool isContained = false;
		size_t check = td[v].bag[i];
		for (int k = 0; k < td[w].bag.size(); k++) {
			if (td[w].bag[k] > check) //make use of ordering of bags to end search early
				break;
			if (td[w].bag[k] == check) {
				isContained = true;
				break;
			}
		}
		if (!isContained)
			introduceMeSize++;
	}


	for (int i = 0; i < td[w].bag.size(); i++) {
		//if some vertex is not in v's bag, add it to forgetMe
		bool isContained = false;
		size_t check = td[w].bag[i];
		for (int k = 0; k < td[v].bag.size(); k++) {
			if (td[v].bag[k] > check) //make use of ordering of bags to end search early
				break;
			if (td[v].bag[k] == check) {
				isContained = true;
				break;
			}
		}
		if (!isContained)
			forgetMeSize++;
	}

	//handle special case that there is nothing to introduce, making v a forget node, and we need to add one less forget node to link them!
	//vector<size_t> curBagIntroduce;
	int skipLast = 0;
	if (introduceMeSize != 0) {
		td[v].type = INTRODUCE;
		//td[v].specialVertex = introduceMe[0];
		//curBagIntroduce = td[v].bag;
		//removeElementFromBag(curBagIntroduce, introduceMe[0]);
	}
	else {
		td[v].type = FORGET;
		//td[v].specialVertex = forgetMe[0];
		skipLast = 1;
	}

	//add forget nodes from bottom up
	td[w].height = td[v].height + forgetMeSize + introduceMeSize;
	//vector<size_t> curBag = td[w].bag;
	TreeDecomposition::vertex_descriptor child = w;
	for (int i = 0; i < forgetMeSize - skipLast; i++) {
		TreeDecomposition::vertex_descriptor newForgetNode = add_vertex(td);
		td[newForgetNode].type = FORGET;
		//td[newForgetNode].bag = removeElementFromBag(curBag, forgetMe[i]); //curBag gets modified here!
		td[newForgetNode].bagSize = td[child].bagSize - 1;
		//td[newForgetNode].specialVertex = forgetMe[i];
		td[newForgetNode].children.push_back(child);
		td[newForgetNode].height = td[v].height + introduceMeSize - 1 + forgetMeSize - i;
		td[child].parent = newForgetNode;
		child = newForgetNode;
	}


	td[v].children.clear(); //need to delete old children!!

	TreeDecomposition::vertex_descriptor curParent = v;
	for (int i = 1; i < introduceMeSize; i++) { //highest introduce node already exists, is v!
		TreeDecomposition::vertex_descriptor newIntroduceNode = add_vertex(td);
		td[newIntroduceNode].type = INTRODUCE;
		//td[newIntroduceNode].bag = curBagIntroduce;
		td[newIntroduceNode].bagSize = td[curParent].bagSize - 1;
		//removeElementFromBag(curBagIntroduce, introduceMe[i]); //curBag gets modified here!
		//td[newIntroduceNode].specialVertex = introduceMe[i];
		td[newIntroduceNode].parent = curParent;
		td[newIntroduceNode].height = td[curParent].height + 1;
		td[curParent].children.push_back(newIntroduceNode);
		curParent = newIntroduceNode;
	}
	//link two chains together
	td[curParent].children.push_back(child);
	td[child].parent = curParent;
}
uint64_t Parser::calculateNumberOfOperations(TreeDecomposition& td, TreeDecomposition::vertex_descriptor root, std::queue<size_t>& leafsLocal) {
	bool* visited = new bool[num_vertices(td)];
	for (int i = 1; i < num_vertices(td); i++) {
		visited[i] = false;
	}
	uint64_t weight = 0;
	while (!leafsLocal.empty()) {
		size_t node = leafsLocal.front();
		while (true) {
			if (td[node].type == 3) {//check if join
				if (!visited[node]) {//if 0th child(always use 0th child) was not checked before, the other thread has not checked it yet, meaning it will then receive a true -> false such that it will not break and keep on going
					visited[node] = true;
					break;
				}
				else {
					td[node].inducedSubgraphSize = td[td[node].children[0]].inducedSubgraphSize + td[td[node].children[1]].inducedSubgraphSize - td[node].bag.size();
				}
			}
			else if (td[node].type == 1) {
				td[node].inducedSubgraphSize = td[td[node].children[0]].inducedSubgraphSize + 1;
			}
			else if (td[node].type == 2) {
				td[node].inducedSubgraphSize = td[td[node].children[0]].inducedSubgraphSize;
			}
			else {
				td[node].inducedSubgraphSize = 1;
			}

			weight += ((uint64_t)1 << (td[node].bagSize - 1)) * (td[node].inducedSubgraphSize - td[node].bagSize + 1);
			//cout << "Node " << node << " contributed " << ((uint64_t)1 << td[node].bagSize) * (td[node].inducedSubgraphSize - td[node].bagSize + 1) << " to the total weight!" << endl;
			if (node == root) {
				break;
			}
			node = td[node].parent;
		}
		leafsLocal.pop();
	}
	delete[] visited;
	return weight;
}

void Parser::retraceCut(TreeDecomposition& td, TreeDecomposition::vertex_descriptor root, uint64_t index) {
	vector<size_t> cut;
	std::queue<NodeIndexPair> q;
	q.push(NodeIndexPair(root, index));
	while (!q.empty()) {
		NodeIndexPair nodeIndex = q.front();
		size_t node = nodeIndex.node;
		uint64_t index = nodeIndex.valueIndex;
		const int degOfFreedom = td[node].inducedSubgraphSize - td[node].bag.size() + 1;
		if (td[node].bag.size() == td[node].inducedSubgraphSize) {//dont add any more children, just add remaining to cut
			int* indices = new int[td[node].bag.size() + 1];
			vector<size_t> sDash = getSdash(td[node].bag, index, degOfFreedom, indices);
			addUncontainedElementsToCut(cut, sDash);
			delete[] indices;
		}
		else {
			int* indices = new int[td[node].bag.size() + 1]; //must suffice even for forget nodes
			vector<size_t> sDash = getSdash(td[node].bag, index, degOfFreedom, indices); //ALSO FILLS IN INDICES
			addUncontainedElementsToCut(cut, sDash);

			if (td[node].type == 1) {
				size_t child = td[node].children[0];
				bool contained = false;
				auto it = sDash.begin();
				while(it != sDash.end()) {
					if (*it == td[node].specialVertex) { 
						contained = true;
						break;
					}
					it++;
				}
				if (contained) {//need to find index where algorithm would have looked!
					int iTable = index % degOfFreedom;
					//auto end = sDash.erase(it); // make new indices
					for (int i = 0; i < sDash.size(); i++) {
						if (td[node].bag[indices[i]] > td[node].specialVertex) {
							indices[i - 1] = indices[i] - 1;
						}
					}
					uint64_t childIndexInKClass = (getIndexOfSubset(indices, sDash.size() - 1) - 1)*degOfFreedom;
					uint64_t subsetsBeforeKClass = 0;
					for (int k = 0; k < sDash.size() - 1; k++) { //sDash.size() == child.bagsize
						subsetsBeforeKClass += binomial(td[child].bag.size(), k);
					}
					subsetsBeforeKClass *= degOfFreedom;
					//cout << "INTRODUCE: Pushing child " << child << " with childIndex " << subsetsBeforeKClass + childIndexInKClass + iTable << endl;
					q.push(NodeIndexPair(child, subsetsBeforeKClass + childIndexInKClass + iTable)); //iTable in child stays the same
				}
				else {
					int iTable = index % degOfFreedom;
					for (int i = 0; i < sDash.size(); i++) {
						if (td[node].bag[indices[i]] > td[node].specialVertex) {
							indices[i]--;
						}
					}
					uint64_t childIndexInKClass = (getIndexOfSubset(indices, sDash.size()) - 1) * degOfFreedom;
					uint64_t subsetsBeforeKClass = 0;
					for (int k = 0; k < sDash.size(); k++) { //sDash.size() == child.bagsize
						subsetsBeforeKClass += binomial(td[child].bag.size(), k);
					}
					subsetsBeforeKClass *= degOfFreedom;
					//cout << "INTRODUCE: Pushing child " << child << " with childIndex " << subsetsBeforeKClass + childIndexInKClass + iTable << endl;
					q.push(NodeIndexPair(child, subsetsBeforeKClass + childIndexInKClass + iTable));
				}
				delete[] indices;
			}
			else if (td[node].type == 2) {
				//first: find out if left or right was taken, then calculate child index accordingly
				size_t child = td[node].children[0];
				uint64_t nodeLength = ((uint64_t)1 << (td[node].bag.size() - 1)) * degOfFreedom;
				size_t forgottenVertex = td[node].specialVertex;

				int forgottenVertexPos = 0;
				while (forgottenVertexPos < td[node].bag.size()) {
					if (forgottenVertex < td[node].bag[forgottenVertexPos])
						break;
					forgottenVertexPos++;
				}

				bool leftOrRight;
				if (index >= nodeLength) {
					leftOrRight = !(*td[node].forget_bitset)[(nodeLength << 1) - index - 1]; //mirror
				}
				else {
					leftOrRight = (*td[node].forget_bitset)[index];
				}

				if (leftOrRight) { // right was taken
					int iTable = (index % degOfFreedom) - 1; //is one less extra within degOfFreedom space!
					int i = sDash.size() - 1;
					for (; i >= 0; i--) {
						if (indices[i] >= forgottenVertexPos) {
							indices[i + 1] = indices[i] + 1;
						}
						else { //i so low that 

							break;
						}
					}
					indices[i + 1] = forgottenVertexPos;

					uint64_t childIndexInKClass = (getIndexOfSubset(indices, sDash.size() + 1) - 1) * (uint64_t)(degOfFreedom - 1);

					uint64_t subsetsBeforeKClass = 0;
					for (int k = 0; k < sDash.size() + 1; k++) {
						subsetsBeforeKClass += binomial(td[child].bag.size(), k);
					}
					subsetsBeforeKClass *= (uint64_t)(degOfFreedom - 1);
					//cout << "FORGET: Pushing child " << child << " with childIndex " << subsetsBeforeKClass << " + " << childIndexInKClass << " + " << iTable << endl;
					q.push(NodeIndexPair(child, subsetsBeforeKClass + childIndexInKClass + iTable));
					
				}
				else { //left was taken
					int iTable = index % degOfFreedom;
					//need to add 1 to every index that points to >= index of forgotten Vertex
					for (int i = 0; i < sDash.size(); i++) {
						if (indices[i] >= forgottenVertexPos) {
							indices[i] += 1;
						}
					}

					uint64_t childIndexInKClass = (getIndexOfSubset(indices, sDash.size()) - 1) * (uint64_t)(degOfFreedom - 1);
					uint64_t subsetsBeforeKClass = 0;
					for (int k = 0; k < sDash.size(); k++) { 
						subsetsBeforeKClass += binomial(td[child].bag.size(), k); 
					}
					subsetsBeforeKClass *= (uint64_t)(degOfFreedom - 1);
					//cout << "FORGET: Pushing child " << child << " with childIndex " << subsetsBeforeKClass + childIndexInKClass + iTable << endl;
					q.push(NodeIndexPair(child, subsetsBeforeKClass + childIndexInKClass + iTable));
				}
			}
			else if (td[node].type == 3) {
				uint64_t nodeLength = ((uint64_t)1 << (td[node].bag.size() - 1)) * degOfFreedom;
				size_t leftChild = td[node].children[0];
				size_t rightChild = td[node].children[1];
				int degFreedomLeft = td[leftChild].inducedSubgraphSize - td[node].bag.size() + 1;
				int degFreedomRight = td[rightChild].inducedSubgraphSize - td[node].bag.size() + 1;
				int i = index % degOfFreedom;
				int j;
				if (index >= nodeLength) {
					//1. calculate range of possible j's, then determine mirroredJ's position in it. then select its inverse respective to that range and use that as the real j.
					int mirroredIndex = (nodeLength << 1) - index - 1;
					int mirroredJ = td[node].join_j[(nodeLength << 1) - index - 1]; //mirror
					int mirroredI = mirroredIndex % degOfFreedom; // same i for mirrored and index position
					int minMirroredJ = max(0, mirroredI - degFreedomRight + 1);
					int distanceFromLeft = mirroredJ - minMirroredJ;
					int maxj = min(i, degFreedomLeft - 1);
					j = maxj - distanceFromLeft;
					
				}
				else {
					j = td[node].join_j[index];
				}
				//first find out j, then calculate both child indices
				
				int notJ = i - j;
				//indices already configured, only find out starting point and add respective j and notJ to push them.
				uint64_t leftIndexUpToSubset = (index / degOfFreedom) * degFreedomLeft;
				uint64_t rightIndexUpToSubset = (index / degOfFreedom) * degFreedomRight;
				//cout << "JOIN: Pushing leftChild " << leftChild << " with childIndex " << leftIndexUpToSubset << " + " << j << endl;
				//cout << "JOIN: Pushing rightChild " << rightChild << " with childIndex " << rightIndexUpToSubset << " + " << notJ << endl;
				q.push(NodeIndexPair(leftChild, leftIndexUpToSubset + j));
				q.push(NodeIndexPair(rightChild, rightIndexUpToSubset + notJ));
			}
			//leaves can't occur because they always will be the end and thus never reach this if else

			for (auto it = td[node].children.begin(); it != td[node].children.end(); it++) {
				//q.push(NodeIndexPair(*it, ))
			}
		}
		q.pop();
	}
	cout << "Finished computing cut [";
	for (auto it = cut.begin(); it != cut.end(); it++)
		cout << *it << ", ";
	cout << "]" << endl;
}

void Parser::addUncontainedElementsToCut(vector<size_t>& cut, const vector<size_t>& sDash) {
	for (auto it = sDash.begin(); it != sDash.end(); it++) {
		bool contained = false;
		for (auto cutIt = cut.begin(); cutIt != cut.end(); cutIt++) {
			if (*it == *cutIt) {
				contained = true;
				break;
			}
		}
		if (!contained) {
			cut.push_back(*it);
		}
	}
}

vector<size_t> Parser::getSdash(const vector<size_t>& bag, uint64_t indexInTotalLength, int degOfFreedom, int* indices) {
	vector<size_t> sDash;
	uint64_t elemsBeforeKClass = 0;
	int choose = 0;
	uint64_t totalIndexNormalized = indexInTotalLength / degOfFreedom; //make use of integer division
	while (elemsBeforeKClass <= totalIndexNormalized) {
		elemsBeforeKClass += binomial(bag.size(), choose++);
	}
	elemsBeforeKClass -= binomial(bag.size(), choose - 1);
	//elemsBeforeKClass *= degOfFreedom;

	uint64_t s = totalIndexNormalized - elemsBeforeKClass + 1;
	
	fillIndices(indices, choose - 1, s);
	for (int i = 0; i < choose - 1; i++) {
		sDash.push_back(bag[indices[i]]);
	}
	return sDash;
}

void Parser::fillIndices(int* indices, int l, int s) {
	int* temp = indices;
	int i = 0;
	while (i < l) { //prefill
		*(temp++) = i++;
	}
	for (int index = l - 1; index >= 0; index--) { //move index in indexes[k-1], then [k-2] .., index = #leftToConsider -> index+1 = k
		int k = index + 1;
		int n = indices[index] + 1; //0->1 shift
		uint64_t bin = 1; //always one because of initialization of indices[]
		while (true) {
			if (bin < s) { //need to skip more
				n++;
				bin *= n;
				bin /= n - k;
				//cout << "Increased n because bin was not sufficient to cover full distance to s. Newly calculated " << n << " / (" << n - k << " ) and gotten new bin=" << bin << endl;
			}
			else { //skipped just once too much, or direct hit. anyway, roll back, save and move on to correct with lower index
				indices[index] = n - 1; //0<-1 shift
				bin *= n - k;
				bin /= n;
				s -= bin;
				break;
			}
		}
	}
}
uint64_t Parser::getSubsetsBeforeKClass(uint64_t normalizedIndex, int bagSize) {
	uint64_t elemsBeforeKClass = 0;
	int choose = 0; //make use of integer division
	while (elemsBeforeKClass <= normalizedIndex) {
		elemsBeforeKClass += binomial(bagSize, choose++);
	}
	elemsBeforeKClass -= binomial(bagSize, choose - 1);
	return elemsBeforeKClass;
}