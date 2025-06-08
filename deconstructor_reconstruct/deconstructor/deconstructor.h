#ifndef DECONSTRUCTOR_H
#define DECONSTRUCTOR_H
#include <iomanip>
#include <iostream>

#include "bond_detector.h"
#include "graph.h"
using namespace std;

class deconstructor {
 private:
  map<int, vector<int>> adjacency_list;
  vector<atom> moleclular_structure;
  vector<node> graph;

 public:
  deconstructor(const string& filename) {
    BondDetector bonds(filename);
    adjacency_list = bonds.inferBonds();
    moleclular_structure = bonds.molecularStructure();
  }

  void showBonds() {
    for (const auto& pair : adjacency_list) {
      cout << pair.first << ": ";
      for (const auto& neighbor : pair.second) {
        cout << neighbor << ",";
      }
      cout << endl;
    }
  }
  void showMolecule() {
    cout << setw(5) << "ID" << setw(8) << "Elem" << setw(10) << "X" << setw(10)
         << "Y" << setw(10) << "Z" << endl;
    cout << string(43, '-') << endl;
    for (const auto& atom : moleclular_structure) {
      cout << setw(5) << atom.resID << setw(8) << atom.element << setw(10)
           << atom.pos.y << setw(10) << atom.pos.z << setw(10) << atom.pos.x
           << endl;
    }
  }
  void generateGraph(atom root) {
    graph = createGraph(adjacency_list, moleclular_structure, root);
  }
  void showGraph() { printGraph(graph); }

  void save() {
    ofstream outfile("graph.txt");
    if (outfile.is_open()) {
      streambuf* coutBuf = cout.rdbuf();

      cout.rdbuf(outfile.rdbuf());

      cout << "ORIGINAL\n";
      showMolecule();
      cout << "\nGRAPH\n";
      showGraph();

      cout.rdbuf(coutBuf);
    } else {
      cerr << "Unable to open file for writing graph." << endl;
    }
  }
};

#endif