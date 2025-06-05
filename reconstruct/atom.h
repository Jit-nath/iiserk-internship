#ifndef ATOM_H
#define ATOM_H

#include "include.h"

class Atom
{
public:
    int id;
    string element;
    double x, y, z;

    Atom(int id = 0, const string &element = "", double x = 0.0, double y = 0.0, double z = 0.0)
        : id(id), element(element), x(x), y(y), z(z) {}

    inline void show()
    {
        cout << "Id\t" << "Element\t" << "X\t" << "Y\t" << "Z" << endl;
        cout << id << "\t" << element << "\t" << x << "\t" << y << "\t" << z << endl;
    }
};

#endif