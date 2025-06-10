#pragma once
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

struct atom {
    unsigned atomID;
    string recordName;
    string atomName;
    unsigned resID;
    double x;
    double y;
    double z;
    double occupancy;
    double tempFactor;

    bool operator<(const atom &other) const {
        return atomID < other.atomID;
    }
    bool operator>(const atom &other) const {
        return atomID > other.atomID;
    }

    void show() {
        cout << "AtomID: " << atomID << ", Name: " << atomName
             << ", ResID: " << resID
             << ", Coords: (" << x << ", " << y << ", " << z << ")\n";
    }
    double distanceTo(const atom &other) const {
        double dx = x - other.x;
        double dy = y - other.y;
        double dz = z - other.z;
        return sqrt(dx * dx + dy * dy + dz * dz);
    }

    void translate(double dx, double dy, double dz) {
        x += dx;
        y += dy;
        z += dz;
    }

    vector<float> direction(const atom &root, const atom &child) const {
        vector<float> dir(3);
        dir[0] = child.x - root.x;
        dir[1] = child.y - root.y;
        dir[2] = child.z - root.z;
        return dir;
    }

    bool isSameResidue(const atom &other) const {
        return (resID == other.resID);
    }
};
