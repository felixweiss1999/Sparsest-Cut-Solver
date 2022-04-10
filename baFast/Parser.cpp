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
			//cout << "detected auto fill" << endl;
			istringstream lineStream(line);
			string trash;
			int width;
			if (!(lineStream >> trash >> trash >> width)) {
				throw std::exception("WARNING: some error occurred while parsing the width of the original graph!");
			}
			fillAdjacencyMatrix(*td, in, width);
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
	//cout << "setup adjacency matrix of width " << width << endl;
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
			const uint64_t degOfFreedom = td[node].inducedSubgraphSize - bagSize + 1;
			uint64_t s_max = (uint64_t)1 << (bagSize - 1);
			uint64_t totalLength = s_max * degOfFreedom;
			td[node].values = new uint32_t [totalLength];//want to be able to choose full bagSize as well
			uint32_t* vals = td[node].values;

			//cout << "thread with id " << std::this_thread::get_id() << " is calculating node " << node << " of type " << td[node].type << " at height " << td[node].height << " with " << ((uint64_t)1 << bagSize) << " * " << degOfFreedom << " values needing to be computed" << endl;


			if (td[node].type == 1) {//introduce
				uint64_t p = 0, c_1 = 0, c_2 = 0;
				size_t child = td[node].children[0];
				size_t intVertexInAdjacencyMatrix = td[node].specialVertex-1;
				int intVertexPos;
				for (int i = 0; i < bagSize; i++) {
					if (td[node].bag[i] == td[node].specialVertex) {
						intVertexPos = i;
						break;
					}
				}

				if (intVertexPos != bagSize - 1) { //no mirror necessary
					uint64_t checkContainednessMask = (uint64_t)1 << intVertexPos;
					for (uint64_t s = 0; s < s_max; s++) {
						if (s & checkContainednessMask) {

							//computeWeight when x in S'
							uint64_t setRep = (~s) & (((uint64_t)1 << bagSize) - 1);
							uint64_t weight = 0, i = 0;
							while (setRep > 0) {
								if (setRep & 1) {
									weight += td[graph_bundle].adjacencyMatrix[intVertexInAdjacencyMatrix][td[node].bag[i] - 1];
								}
								setRep = setRep >> 1;
								i++;
							}


							for (uint64_t i = 0; i < degOfFreedom; i++)
								vals[p++] = td[child].values[c_1++] + weight;
						}
						else {

							//compute weight for x not in S'
							uint64_t setRep = s;
							uint64_t weight = 0, i = 0;
							while (setRep > 0) {
								if (setRep & 1) {
									weight += td[graph_bundle].adjacencyMatrix[intVertexInAdjacencyMatrix][td[node].bag[i] - 1];
								}
								setRep = setRep >> 1;
								i++;
							}


							for (uint64_t i = 0; i < degOfFreedom; i++)
								vals[p++] = td[child].values[c_2++] + weight;
						}
					}
				}
				else { //mirror necessary, and also only x not in S' possible
					s_max = s_max >> 1;
					for (uint64_t s = 0; s < s_max; s++) {

						//compute weight for x not in S'
						uint64_t setRep = s;
						uint64_t weight = 0, i = 0;
						while (setRep > 0) {
							if (setRep & 1) {
								weight += td[graph_bundle].adjacencyMatrix[intVertexInAdjacencyMatrix][td[node].bag[i] - 1];
							}
							setRep = setRep >> 1;
							i++;
						}


						for (uint64_t i = 0; i < degOfFreedom; i++)
							vals[p++] = td[child].values[c_2++] + weight;
					}
					s_max = s_max << 1;
					for (uint64_t s = s_max >> 1; s < s_max; s++) { //second half of child indices needs to be mirrored

						//compute weight for x not in S'
						uint64_t setRep = s;
						uint64_t weight = 0, i = 0;
						while (setRep > 0) {
							if (setRep & 1) {
								weight += td[graph_bundle].adjacencyMatrix[intVertexInAdjacencyMatrix][td[node].bag[i] - 1];
							}
							setRep = setRep >> 1;
							i++;
						}


						for (uint64_t i = 0; i < degOfFreedom; i++)
							vals[p++] = td[child].values[--c_2] + weight;
					}
				}
			}
			else if (td[node].type == 2) {
				td[node].forget_bitset = new boost::dynamic_bitset<uint64_t>(totalLength);
				size_t child = td[node].children[0];
				int frgtVertexPos;
				for (int i = 0; i < bagSize + 1; i++) { //CHILD bag is being searched!
					if (td[child].bag[i] == td[node].specialVertex) {
						frgtVertexPos = i;
						break;
					}
				}
				uint64_t leftPartMask = UINT64_MAX << frgtVertexPos;
				uint64_t rightPartMask = ((uint64_t)1 << frgtVertexPos) - 1;
				uint64_t forgottenVertexBit = (uint64_t)1 << frgtVertexPos;

				uint64_t p = 0;
				if (frgtVertexPos != bagSize) { //bagsize = |X_u| = |X_v| - 1. In this case, no mirroring ever necessary
					for (uint64_t s = 0; s < s_max; s++) {
						uint64_t s_child = ((leftPartMask & s) << 1) | (rightPartMask & s); //inserted 0 at forgottenVertexPos
						uint64_t leftBaseIndex = (degOfFreedom - 1) * s_child;
						s_child = forgottenVertexBit | s_child;
						uint64_t rightBaseIndex = (degOfFreedom - 1) * s_child;

						//(*td[node].forget_bitset)[p] = false; unnecessary because always init to 0!
						vals[p++] = td[child].values[leftBaseIndex];
						for (uint64_t i = 1; i < degOfFreedom - 1; i++) {
							uint64_t left = td[child].values[leftBaseIndex + i];
							uint64_t right = td[child].values[rightBaseIndex + i - 1];
							if (right < left) {
								(*td[node].forget_bitset)[p] = true;
								vals[p++] = right;
							}
							else {
								vals[p++] = left;
							}
						}
						(*td[node].forget_bitset)[p] = true;
						vals[p++] = td[child].values[rightBaseIndex + degOfFreedom - 2];
					}
				}
				else { //mirroring may be necessary
					uint64_t childTotalLength = (s_max << 1) * (degOfFreedom - 1);
					for (uint64_t s = 0; s < s_max; s++) {
						uint64_t s_child = ((leftPartMask & s) << 1) | (rightPartMask & s); //inserted 0 at forgottenVertexPos
						uint64_t leftBaseIndex = (degOfFreedom - 1) * s_child;
						s_child = forgottenVertexBit | s_child;
						uint64_t rightBaseIndex = (degOfFreedom - 1) * s_child;

						//(*td[node].forget_bitset)[p] = false; unnecessary because always init to 0!
						vals[p++] = td[child].values[leftBaseIndex];
						for (uint64_t i = 1; i < degOfFreedom - 1; i++) {
							uint64_t left = td[child].values[leftBaseIndex + i];

							uint64_t rightIndex = rightBaseIndex + i - 1; //mirror
							if (rightIndex >= childTotalLength) {
								rightIndex = (childTotalLength << 1) - rightIndex - 1;
							}

							uint64_t right = td[child].values[rightIndex];
							if (right < left) {
								(*td[node].forget_bitset)[p] = true;
								vals[p++] = right;
							}
							else {
								vals[p++] = left;
							}
						}

						uint64_t rightIndex = rightBaseIndex + degOfFreedom - 2; //mirror
						if (rightIndex >= childTotalLength) {
							rightIndex = (childTotalLength << 1) - rightIndex - 1;
						}

						(*td[node].forget_bitset)[p] = true;
						vals[p++] = td[child].values[rightIndex];
					}
				}
			}
			else if (td[node].type == 3) {
				td[node].join_j = new uint16_t[totalLength];
				size_t leftChild = td[node].children[0];
				size_t rightChild = td[node].children[1];
				int degFreedomLeft = td[leftChild].inducedSubgraphSize - bagSize + 1;
				int degFreedomRight = td[rightChild].inducedSubgraphSize - bagSize + 1;

				uint64_t p = 0, leftBaseIndex = 0, rightBaseIndex = 0;
				for (uint64_t s = 0; s < s_max; s++) {
					
					//compute weight
					uint64_t weight = 0;
					uint64_t leftRep = s;
					uint64_t a = 0;
					while (leftRep > 0) {
						if (leftRep & 1) {
							uint64_t rightRep = (~s) & (((uint64_t)1 << bagSize) - 1);
							uint64_t b = 0;
							while (rightRep > 0) {
								if (rightRep & 1) {
									weight += td[graph_bundle].adjacencyMatrix[td[node].bag[a] - 1][td[node].bag[b] - 1];
								}
								rightRep = rightRep >> 1;
								b++;
							}
						}
						leftRep = leftRep >> 1;
						a++;
					}


					for (uint64_t i = 0; i < degOfFreedom; i++) {
						uint64_t min = UINT64_MAX;
						int j_lowerBound = i - degFreedomRight + 1;
						if (j_lowerBound < 0)
							j_lowerBound = 0;
						int j_upperBound = degFreedomLeft - 1;
						if (j_upperBound > i)
							j_upperBound = i;
						for (int j = j_lowerBound; j <= j_upperBound; j++) {
							uint64_t candidate = (uint64_t)td[leftChild].values[leftBaseIndex + j] + (uint64_t)td[rightChild].values[rightBaseIndex + i - j];
							if (candidate < min) {
								min = candidate;
								td[node].join_j[p] = j - j_lowerBound;
							}
						}
						vals[p++] = min - weight;
					}
					leftBaseIndex += degFreedomLeft;
					rightBaseIndex += degFreedomRight;
				}
			}
			else {//leaf
				vals[0] = 0; 
			}
			for (auto it = td[node].children.begin(); it != td[node].children.end(); it++) {
				delete[] td[*it].values;
			}
			if (node == root) {
				//cout << "thread with id " << std::this_thread::get_id() << " has reached root and will now cease activity!" << endl;
				uint64_t rootIndex = calculateCutWeight(td, root);
				delete[] td[root].values;
				retraceCut(td, root, rootIndex);
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



uint64_t Parser::calculateCutWeight(TreeDecomposition& td, TreeDecomposition::vertex_descriptor root) {
	const int bagSize = td[root].bag.size();
	const int degOfFreedom = td[root].inducedSubgraphSize - bagSize + 1;
	uint32_t* vals = td[root].values;
	//int minWeight = 0;
	uint64_t p = 0;
	uint64_t pOfMin = 0;
	//double lowestI, lowestNotI;
	double minSparsestCutDensity = DBL_MAX;
	uint64_t s_max = (uint64_t)1 << (bagSize - 1);
	for (int s = 0; s < s_max; s++) {
		int k = std::_Checked_x86_x64_popcount(s);
		for (int i = k; i < k + degOfFreedom; i++) {
			int notI = td[root].inducedSubgraphSize - i;
			if (!(notI && i)) {
				p++;
				continue;
			}
			double candidateDensity = (double)vals[p++] / (i * (notI));
			//cout << candidate << endl;
			if (candidateDensity < minSparsestCutDensity) {
				minSparsestCutDensity = candidateDensity;
				//minWeight = vals[p - 1];
				//lowestI = i;
				//lowestNotI = notI;
				pOfMin = p - 1;
			}

		}
	}
	//cout << "minWeight: " << minWeight << " minSparsestCutWeight: " << minSparsestCutDensity << " selected index in top table: " << pOfMin << " i: " << lowestI <<" notI: " << lowestNotI << endl;

	return pOfMin;
}


size_t Parser::calculateOptimalRoot(TreeDecomposition& td) {
	uint64_t weight;
	//uint64_t maxWeight = 0;
	uint64_t minWeight = UINT64_MAX;
	//size_t worstRoot = 1;
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
		//if (weight > maxWeight) {
		//	worstRoot = *rootCandidate;
		//	maxWeight = weight;
		//}
	}
	//cout << "Best root was node " << bestRoot << " with a weight of " << minWeight << endl;
	//cout << "Worst root was node " << worstRoot << " with a weight of " << maxWeight << endl;
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
			if (td[node].type == 3) {
				uint64_t jAmount = (td[td[node].children[0]].inducedSubgraphSize - td[td[node].children[0]].bagSize + 1) * (td[td[node].children[1]].inducedSubgraphSize - td[td[node].children[1]].bagSize + 1);
				weight += ((uint64_t)1 << (td[node].bagSize - 1)) * jAmount;
			}
			else {
				weight += ((uint64_t)1 << (td[node].bagSize - 1)) * (td[node].inducedSubgraphSize - td[node].bagSize + 1);
			}
			
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
	
	uint64_t rootRep = index / (td[root].inducedSubgraphSize - td[root].bag.size() + 1); //add root S' to cut
	int i = 0;
	while (rootRep > 0) {
		if (rootRep & 1) {
			cut.push_back(td[root].bag[i]);
		}
		rootRep = rootRep >> 1;
		i++;
	}
	std::queue<NodeIndexPair> q;
	q.push(NodeIndexPair(root, index));

	while (!q.empty()) {
		NodeIndexPair nodeIndex = q.front();
		size_t node = nodeIndex.node;
		uint64_t index = nodeIndex.valueIndex;
		const uint64_t degOfFreedom = td[node].inducedSubgraphSize - td[node].bag.size() + 1;
		
		if (node != root && td[td[node].parent].type == 2) {
			size_t x = td[td[node].parent].specialVertex;
			int intVertexPos;
			for (int i = 0; i < td[node].bag.size(); i++) {
				if (td[node].bag[i] == x) {
					intVertexPos = i;
					break;
				}
			}
			if (((uint64_t)index / degOfFreedom) & ((uint64_t)1 << intVertexPos)) {
				cut.push_back(x);
			}
		}

		if (td[node].type == 1) {
			int intVertexPos;
			for (int i = 0; i < td[node].bag.size(); i++) {
				if (td[node].bag[i] == td[node].specialVertex) {
					intVertexPos = i;
					break;
				}
			}
			uint64_t i = index % degOfFreedom;
			uint64_t s = index / degOfFreedom;
			uint64_t leftPart = s & (UINT64_MAX << (intVertexPos + 1));
			leftPart = leftPart >> 1;
			s = s & (((uint64_t)1 << intVertexPos) - 1);
			q.push(NodeIndexPair(td[node].children[0], degOfFreedom * (leftPart | s) + i));
			//cout << "INTRODUCE: Pushed child " << td[node].children[0] << " with childIndex " << degOfFreedom * (leftPart | s) + i << endl;
		}
		else if (td[node].type == 2) {
			int frgtVertexPos;
			size_t child = td[node].children[0];
			for (int i = 0; i < td[child].bag.size(); i++) {
				if (td[child].bag[i] == td[node].specialVertex) {
					frgtVertexPos = i;
					break;
				}
			}

			bool decision;
			if (index >= degOfFreedom * ((uint64_t)1 << (td[node].bag.size() - 1))) {
				//cout << "flag1" << endl;
				decision = !(*td[node].forget_bitset)[degOfFreedom * ((uint64_t)1 << (td[node].bag.size())) - index - 1];
			}
			else {
				//cout << "flag2" << endl;
				decision = (*td[node].forget_bitset)[index];
			}

			uint64_t i = index % degOfFreedom;
			uint64_t s = index / degOfFreedom;
			uint64_t leftPart = (s & (UINT64_MAX << frgtVertexPos)) << 1;
			s = s & (((uint64_t)1 << frgtVertexPos) - 1);
			uint64_t child_s = leftPart | s; //inserted 0
			if (decision) { // right value was taken
				child_s = child_s | ((uint64_t)1 << frgtVertexPos);
				q.push(NodeIndexPair(child, (degOfFreedom - 1) * child_s + i - 1));
				//cout << "ForgetRight: Pushed child " << child << " with childIndex " << (degOfFreedom - 1) * child_s + i - 1 << endl;
			}
			else {
				q.push(NodeIndexPair(child, (degOfFreedom - 1) * child_s + i));
				//cout << "ForgetLeft: Pushed child " << child << " with childIndex " << (degOfFreedom - 1) * child_s + i << endl;
			}
		}
		else if (td[node].type == 3) {
			size_t leftChild = td[node].children[0];
			size_t rightChild = td[node].children[1];
			uint64_t degFreedomLeft = td[leftChild].inducedSubgraphSize - td[node].bag.size() + 1;
			uint64_t degFreedomRight = td[rightChild].inducedSubgraphSize - td[node].bag.size() + 1;
			uint64_t i = index % degOfFreedom;
			uint64_t s = index / degOfFreedom;
			uint64_t j;
			if (index >= degOfFreedom * ((uint64_t)1 << (td[node].bag.size() - 1))) {
				//cout << "flag3" << endl;
				j = std::min(i, degFreedomLeft - 1) - td[node].join_j[(uint64_t)degOfFreedom * ((uint64_t)1 << td[node].bag.size()) - index - 1];
			}
			else {
				j = td[node].join_j[index] + std::max((long long int)0, (long long int)i - (long long int)degFreedomRight + 1); //careful with types!
				//cout << "flag4" << endl;
			}
			q.push(NodeIndexPair(leftChild, degFreedomLeft * s + j));
			//cout << "JOIN_1: Pushed child " << leftChild << " with childIndex " << degFreedomLeft * s + j << endl;
			q.push(NodeIndexPair(rightChild, degFreedomRight * s + i - j));
			//cout << "JOIN_2: Pushed child " << rightChild << " with childIndex " << degFreedomRight * s + i - j << endl;
		}

		q.pop();
	}
	//cout << "Finished computing cut [";
	for (auto it = cut.begin(); it != cut.end(); it++)
		cout << *it << endl;
	//cout << "]" << endl;
}


