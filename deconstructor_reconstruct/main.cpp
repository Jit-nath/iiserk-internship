#include "deconstructor/deconstructor.h"

int main() {
  deconstructor molecule("./lipid.pdb");

  molecule.generateGraph(atom(1, "N", -1.306, 0.669, 20.620));

  molecule.save();

  return 0;
}