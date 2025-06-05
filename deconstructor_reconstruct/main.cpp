#include "atom.h"
#include "bond_detector.h"
#include "include.h"
#include "lipid.h"
#include "reconstructor.h"

int main()
{
    lipid molecule("lipid.pdb");
    Atom atom = molecule.getAtomById(1);
    atom.show();
    atom.show_pdb_format();

    return 0;
}