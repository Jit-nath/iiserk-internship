#ifndef ADJ_MAT_GEN_H
#define ADJ_MAT_GEN_H

#include <cmath>
#include <fstream>
#include <map>
#include <vector>

#include "atom.h"
using namespace std;

class BondDetector {
 private:
  vector<atom> atoms;
  map<int, vector<int>> bond_list;

  bool areBonded(atom& a1, atom& a2) {
    double dx = a1.pos.x - a2.pos.x;
    double dy = a1.pos.y - a2.pos.y;
    double dz = a1.pos.z - a2.pos.z;
    return sqrt(dx * dx + dy * dy + dz * dz) < 2.0;
  }

 public:
  BondDetector(const string& filename) {
    ifstream file(filename);
    string line, atomName;
    int atomId;
    double x, y, z;
    atoms.clear();

    while (getline(file, line)) {
      if (line.substr(0, 4) == "ATOM") {
        atomId = stoi(line.substr(6, 5));
        atomName = line.substr(12, 4);
        x = stod(line.substr(30, 8));
        y = stod(line.substr(38, 8));
        z = stod(line.substr(46, 8));

        atoms.emplace_back(atomId, atomName, x, y, z);
      }
    }
    file.close();
  }

  vector<atom> molecularStructure() { return atoms; }

  map<int, vector<int>> inferBonds() {
    bond_list.clear();

    for (size_t i = 0; i < atoms.size(); i++) {
      for (size_t j = i + 1; j < atoms.size(); j++) {
        if (areBonded(atoms[i], atoms[j])) {
          bond_list[atoms[i].resID].push_back(atoms[j].resID);
          bond_list[atoms[j].resID].push_back(atoms[i].resID);
        }
      }
    }
    return bond_list;
  }
};

#endif