#ifndef ADJ_MAT_GEN_H
#define ADJ_MAT_GEN_H

#include "atom.h"
#include "include.h"

class BondDetector {
private:
    vector<Atom> atoms;
    map<int, vector<int>> adjacencyList;

    double calculateDistance(const Atom& a1, const Atom& a2)
    {
        double dx = a1.coords.x - a2.coords.x;
        double dy = a1.coords.y - a2.coords.y;
        double dz = a1.coords.z - a2.coords.z;
        return sqrt(dx * dx + dy * dy + dz * dz);
    }

    string getElementFromName(const string& atomName)
    {
        if (atomName.empty())
            return "";

        if (atomName[0] == 'H')
            return "H";
        if (atomName[0] == 'C')
            return "C";
        if (atomName[0] == 'N')
            return "N";
        if (atomName[0] == 'O')
            return "O";
        if (atomName[0] == 'P')
            return "P";
        if (atomName[0] == 'S')
            return "S";

        if (atomName.length() >= 2) {
            string twoChar = atomName.substr(0, 2);
            if (twoChar == "Cl" || twoChar == "Br" || twoChar == "Ca" || twoChar == "Mg" || twoChar == "Na" || twoChar == "Fe") {
                return twoChar;
            }
        }

        return string(1, atomName[0]);
    }

    bool areBonded(const Atom& a1, const Atom& a2)
    {
        double distance = calculateDistance(a1, a2);
        return distance < 2.0;
    }

    bool readPDBFile(const string& filename)
    {
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Error: Could not open file " << filename << endl;
            return false;
        }

        string line;
        atoms.clear();

        while (getline(file, line)) {
            if (line.length() < 54)
                continue; // Skip lines that are too short

            if (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM") {
                try {
                    int atomId = stoi(line.substr(6, 5));
                    string atomName = line.substr(12, 4);
                    atomName.erase(0, atomName.find_first_not_of(" \t"));
                    atomName.erase(atomName.find_last_not_of(" \t") + 1);

                    double x = stod(line.substr(30, 8));
                    double y = stod(line.substr(38, 8));
                    double z = stod(line.substr(46, 8));

                    string element = getElementFromName(atomName);

                    atoms.emplace_back(atomId, element, x, y, z);
                } catch (const exception& e) {
                    cerr << "Error parsing line: " << line << endl;
                    continue;
                }
            }
        }

        file.close();
        return !atoms.empty();
    }

public:
    BondDetector(const string& filename)
    {
        if (!readPDBFile(filename)) {
            cerr << "Failed to read PDB file: " << filename << endl;
        }
    }

    map<int, vector<int>> createBondInfo()
    {
        adjacencyList.clear();

        for (const auto& atom : atoms) {
            adjacencyList[atom.id] = vector<int>();
        }

        for (size_t i = 0; i < atoms.size(); i++) {
            for (size_t j = i + 1; j < atoms.size(); j++) {
                if (areBonded(atoms[i], atoms[j])) {
                    adjacencyList[atoms[i].id].push_back(atoms[j].id);
                    adjacencyList[atoms[j].id].push_back(atoms[i].id);
                }
            }
        }

        return adjacencyList;
    }

    void show() const
    {
        if (adjacencyList.empty()) {
            cout << "No connectivity data available." << endl;
            return;
        }

        cout << "\n=== Molecular Connectivity ===\n";
        cout << "Total atoms: " << adjacencyList.size() << "\n\n";

        for (const auto& [atomId, connections] : adjacencyList) {
            cout << "Atom " << setw(3) << atomId << " -> ";

            if (connections.empty()) {
                cout << "No bonds";
            } else {
                cout << "Connected to: ";
                for (size_t i = 0; i < connections.size(); ++i) {
                    cout << connections[i];
                    if (i < connections.size() - 1) {
                        cout << ", ";
                    }
                }
                cout << " (" << connections.size() << " bond"
                     << (connections.size() == 1 ? "" : "s") << ")";
            }
            cout << "\n";
        }

        cout << "\n"
             << string(30, '=') << "\n";
    }

    // Alternative version with more detailed output
    void showDetailed() const
    {
        if (adjacencyList.empty()) {
            cout << "No connectivity data available." << endl;
            return;
        }

        cout << "\n"
             << string(50, '=') << "\n";
        cout << "           MOLECULAR CONNECTIVITY REPORT\n";
        cout << string(50, '=') << "\n";

        int totalBonds = 0;
        int isolatedAtoms = 0;

        for (const auto& [atomId, connections] : adjacencyList) {
            totalBonds += connections.size();
            if (connections.empty())
                isolatedAtoms++;
        }

        cout << "Summary:\n";
        cout << "  Total atoms: " << adjacencyList.size() << "\n";
        cout << "  Total bonds: " << totalBonds / 2 << "\n"; // Divide by 2 since each bond is counted twice
        cout << "  Isolated atoms: " << isolatedAtoms << "\n\n";

        cout << "Detailed Connectivity:\n";
        cout << string(50, '-') << "\n";

        for (const auto& [atomId, connections] : adjacencyList) {
            cout << "Atom " << setw(4) << atomId << " | ";

            if (connections.empty()) {
                cout << "ISOLATED (no bonds)";
            } else {
                cout << "Bonds: ";
                for (size_t i = 0; i < connections.size(); ++i) {
                    cout << connections[i];
                    if (i < connections.size() - 1) {
                        cout << ", ";
                    }
                }
                cout << " [" << connections.size() << "]";
            }
            cout << "\n";
        }

        cout << string(50, '=') << "\n\n";
    }

    // Compact version for quick overview
    void showCompact() const
    {
        if (adjacencyList.empty()) {
            cout << "No data.\n";
            return;
        }

        cout << "Connectivity (" << adjacencyList.size() << " atoms):\n";
        for (const auto& [atomId, connections] : adjacencyList) {
            if (!connections.empty()) {
                cout << atomId << ": ";
                for (size_t i = 0; i < connections.size(); ++i) {
                    cout << connections[i] << (i < connections.size() - 1 ? "," : "");
                }
                cout << "\n";
            }
        }
    }

    // Modern C++ version using ranges (C++20)
    void showModern() const
    {
        if (adjacencyList.empty()) {
            cout << "No connectivity data available.\n";
            return;
        }

        cout << "\n=== Molecular Connectivity ===\n";

        for (const auto& [atomId, connections] : adjacencyList) {
            cout << "Atom " << setw(3) << atomId << " -> ";

            if (connections.empty()) {
                cout << "No bonds\n";
                continue;
            }

            // Join connections with comma separator
            ostringstream oss;
            for (auto it = connections.begin(); it != connections.end(); ++it) {
                if (it != connections.begin())
                    oss << ", ";
                oss << *it;
            }

            cout << oss.str() << " (" << connections.size() << " bonds)\n";
        }
        cout << string(30, '=') << "\n";
    }
};

#endif