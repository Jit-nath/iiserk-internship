#include "atom.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

class pdb {
    fstream file;
    string file_line;
    vector<atom> atoms;

  public:
    pdb() {}
    pdb(string filename) {
        file.open(filename);
    }
    ~pdb() { file.close(); }

    void parse() {
        while (getline(file, file_line)) {
            string record_name = file_line.substr(0, 6);
            if (record_name != "ATOM  " && record_name != "HETATM")
                continue;

            try {
                atom single_atom;
                single_atom.atomID = stoi(file_line.substr(6, 5));
                single_atom.atomName = file_line.substr(12, 4);
                single_atom.resID = stoi(file_line.substr(22, 4));
                single_atom.recordName = record_name;
                single_atom.x = stod(file_line.substr(30, 8));
                single_atom.y = stod(file_line.substr(38, 8));
                single_atom.z = stod(file_line.substr(46, 8));
                single_atom.occupancy = stod(file_line.substr(54, 6));
                single_atom.tempFactor = stod(file_line.substr(60, 6));

                atoms.push_back(single_atom);
            } catch (const std::exception &e) {
                cerr << "Failed to parse line: " << file_line << endl;
                cerr << "Reason: " << e.what() << endl;
            }
        }
    }

    const vector<atom> &getAtoms() const {
        return atoms;
    }

    atom getByAtomID(unsigned num) {
        for (const auto &a : atoms) {
            if (a.atomID == num) {
                return a;
            }
        }
        throw out_of_range("Atom with given ID not found");
    }

    void printAtom(unsigned num) {

        const atom printable = getByAtomID(num);
        cout << left
             << setw(12) << printable.recordName
             << setw(8) << printable.atomID
             << setw(12) << printable.atomName
             << setw(8) << printable.resID
             << setw(10) << printable.x
             << setw(10) << printable.y
             << setw(10) << printable.z
             << setw(12) << printable.occupancy
             << setw(12) << printable.tempFactor
             << endl;
    }

    void show() {
        cout << left
             << setw(12) << "RecordName"
             << setw(8) << "AtomID"
             << setw(12) << "AtomName"
             << setw(8) << "ResID"
             << setw(10) << "x"
             << setw(10) << "y"
             << setw(10) << "z"
             << setw(12) << "Occupancy"
             << setw(12) << "TempFactor"
             << endl;

        for (const auto &a : atoms) {
            cout << left
                 << setw(12) << a.recordName
                 << setw(8) << a.atomID
                 << setw(12) << a.atomName
                 << setw(8) << a.resID
                 << setw(10) << a.x
                 << setw(10) << a.y
                 << setw(10) << a.z
                 << setw(12) << a.occupancy
                 << setw(12) << a.tempFactor
                 << endl;
        }
    }
};
