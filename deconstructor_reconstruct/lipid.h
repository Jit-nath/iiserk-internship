#ifndef LIPID_H
#define LIPID_H

#include "include.h"
#include "atom.h"
#include "bond_detector.h"

class lipid
{
private:
    fstream file;
    map<int, vector<int>> adjacency_list;
    string file_contents;

public:
    lipid(const string &filename)
    {
        file.open(filename);
        if (!file.is_open())
        {
            cerr << "Error: Could not open file " << filename << endl;
            return;
        }

        file_contents.assign((istreambuf_iterator<char>(file)), istreambuf_iterator<char>());
        file.close();

        BondDetector generator(filename);
        adjacency_list = generator.createBondInfo();
    }

    Atom getAtomById(int id)
    {
        istringstream stream(file_contents);
        string line;

        while (getline(stream, line))
        {
            if (line.length() < 54)
                continue; // Skip lines that are too short

            if (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM")
            {
                try
                {
                    int atomId = stoi(line.substr(6, 5));
                    if (atomId != id)
                        continue;

                    string atomName = line.substr(12, 4);
                    atomName.erase(0, atomName.find_first_not_of(" \t"));
                    atomName.erase(atomName.find_last_not_of(" \t") + 1);

                    double x = stod(line.substr(30, 8));
                    double y = stod(line.substr(38, 8));
                    double z = stod(line.substr(46, 8));

                    return Atom(atomId, atomName, x, y, z);
                }
                catch (const exception &e)
                {
                    cerr << "Error parsing atom with ID " << id << endl;
                }
            }
        }

        // Return default atom if not found
        cerr << "Warning: Atom with ID " << id << " not found" << endl;
        return Atom();
    }

    const map<int, vector<int>> &getAdjacencyList() const
    {
        return adjacency_list;
    }
    /*
    TODO : generate graph, Show Graph
    */
};

#endif