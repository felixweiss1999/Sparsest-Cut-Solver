// baFast.cpp : Diese Datei enthält die Funktion "main". Hier beginnt und endet die Ausführung des Programms.
//

#include <boost/graph/adjacency_list.hpp>
#include <iostream>
#include "TreeDecomposition.h"
#include "Parser.h"

#include <string>
#include <sstream>
#include <fstream>
#include <chrono>
using namespace boost;

int main()
{
    Parser p;
    //std::string a = "s td 4 5 4\nb 1 1 2 3 4\nb 2 2 3 4\nb 3\nb 4\n1 2\n1 3\n2 4";
    //std::istringstream s(a);

    try {

        std::ifstream file;
        file.open("testInstances/ex001.td");
        auto start_time = std::chrono::high_resolution_clock::now();
        TreeDecomposition td = *p.parse(file);
        auto end_time = std::chrono::high_resolution_clock::now();
        std::cout << (end_time - start_time) / std::chrono::microseconds(1) << " microseconds to find best root " << std::endl;
        file.close();
        
        /*std::ofstream writeFile;
        writeFile.open("test.td");
        p.exportDimax(td, writeFile);
        writeFile.close();*/
        file.open("testInstances/ex001.gr");
        p.fillAdjacencyMatrix(td, file);
        file.close();
        
        std::cout << "Ready?" << std::endl;
        std::cin.get();
        start_time = std::chrono::high_resolution_clock::now();
        p.debugAlgorithm(td);
        end_time = std::chrono::high_resolution_clock::now();
        std::cout << (end_time - start_time) / std::chrono::microseconds(1) << " microseconds to run computation with root " << td[graph_bundle].root << std::endl;
        //p.print(td, td[graph_bundle].root);
    }
    catch (const std::exception& e) {
        std::cout << e.what();
    }
}