#include "include.h"
#include "atom.h"
#include "bond_detector.h"
#include "lipid.h"
#include "reconstructor.h"

int main()
{
    lipid molecule("lipid.pdb");
    Atom atom = molecule.getAtomById(1);
    atom.show();
    BondDetector bonds("lipid.pdb");
    bonds.createBondInfo();
    bonds.show();
    return 0;
}