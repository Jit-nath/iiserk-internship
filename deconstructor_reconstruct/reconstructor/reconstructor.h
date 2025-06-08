#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H
#include "atom.h"

class reconstructor {
 private:
  string inp_molecule;
  pos inp_pos;
  pos inp_direction;
  map<int, vector<int>> structure;

 public:
  reconstructor(const string& name, const pos& molecule_pos,
                const pos& molecule_direction)
      : inp_molecule(name),
        inp_pos(molecule_pos),
        inp_direction(molecule_direction) {}

  void loadStructure(string filename) {
    ifstream file(filename);
    string line;
    while (getline(file, line)) {
    }
  }
  void reconstruct() {}
  void show() {}
  void saveToPdb() {}
};

void PrintAlogrithm() {
  cout << "                        ┌───────────────────────────┐" << endl;
  cout << "                        │        Algorithm          │" << endl;
  cout << "                        └───────────────────────────┘" << endl;
  cout << "┌───────────────────────────────────────────────────────────────────"
          "─────────┐"
       << endl;
  cout << "│ 1. Take the Atom name, position [x, y, z] and direction vector "
          "[x, y, z]   │"
       << endl;
  cout << "│ 2. Take the reference graph which we want to reconstruct          "
          "         │"
       << endl;
  cout << "│ 3. Read the reference graph                                       "
          "         │"
       << endl;
  cout << "│ 4. Initialize a map {Atom id : {bond length, bond angle}}         "
          "         │"
       << endl;
  cout << "│                                                                   "
          "         │"
       << endl;
  cout << "└───────────────────────────────────────────────────────────────────"
          "─────────┘"
       << endl;
}

#endif