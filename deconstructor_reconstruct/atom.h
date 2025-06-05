#ifndef ATOM_H
#define ATOM_H

#include "include.h"

struct pos {
    double x, y, z;

    // Constructor
    pos(double x = 0.0, double y = 0.0, double z = 0.0)
        : x(x)
        , y(y)
        , z(z)
    {
    }
};

class Atom {
public:
    int id;
    string element;
    pos coords;

    Atom(int id = 0, string element = "", double x = 0.0, double y = 0.0, double z = 0.0)

    {
        this->id = id;
        this->element = element;
        this->coords = { x, y, z };
    }

    inline void show()
    {
        cout << "Id\t" << "Element\t" << "X\t" << "Y\t" << "Z" << endl;
        cout << id << "\t" << element << "\t" << coords.x << "\t" << coords.y << "\t" << coords.z << endl;
    }

    inline void show_pdb_format()
    {
        cout << left << setw(6) << "ATOM"
             << right << setw(5) << id << "  "
             << left << setw(4) << element
             << "POPCX" << setw(2) << " " << "1"
             << fixed << setprecision(3)
             << setw(11) << coords.x
             << setw(8) << coords.y
             << setw(8) << coords.z
             << setw(6) << "1.00"
             << setw(6) << "0.00"
             << "      MEMB"
             << endl;
    }
};

#endif