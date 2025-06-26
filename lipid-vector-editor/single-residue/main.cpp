#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

struct atom {
    int atomID;
    string atomName;
    string residueName;
    int resID;
    double x;
    double y;
    double z;
    double occupancy;
    double tempFactor;
    string chainID;

    bool operator<(const atom &other) const {
        return atomID < other.atomID;
    }
};

// Fixed distance calculation
float dist(atom &a1, atom &a2) {
    float dx = a1.x - a2.x;
    float dy = a1.y - a2.y;
    float dz = a1.z - a2.z;

    float distance = sqrt(dx * dx + dy * dy + dz * dz);
    return distance;
}

// Calculate direction vector from root to child
vector<float> direction(atom root, atom child) {
    vector<float> dir(3);
    dir[0] = child.x - root.x;
    dir[1] = child.y - root.y;
    dir[2] = child.z - root.z;
    return dir;
}

// Normalize a vector
vector<float> normalize(vector<float> vec) {
    float magnitude = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    if (magnitude == 0)
        return vec;

    vector<float> normalized(3);
    normalized[0] = vec[0] / magnitude;
    normalized[1] = vec[1] / magnitude;
    normalized[2] = vec[2] / magnitude;
    return normalized;
}

// Calculate dot product of two vectors
float dotProduct(vector<float> a, vector<float> b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Calculate cross product of two vectors
vector<float> crossProduct(vector<float> a, vector<float> b) {
    vector<float> result(3);
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    return result;
}

// Create rotation matrix to rotate from original direction to new direction
vector<vector<float>> createRotationMatrix(vector<float> from, vector<float> to) {
    // Normalize both vectors
    from = normalize(from);
    to = normalize(to);

    // Initialize 3x3 identity matrix
    vector<vector<float>> rotMatrix(3, vector<float>(3, 0));
    rotMatrix[0][0] = rotMatrix[1][1] = rotMatrix[2][2] = 1.0;

    float cosTheta = dotProduct(from, to);

    // If vectors are already aligned, return identity matrix
    if (abs(cosTheta - 1.0) < 1e-6) {
        return rotMatrix;
    }

    // If vectors are opposite, create 180-degree rotation
    if (abs(cosTheta + 1.0) < 1e-6) {
        // Find a perpendicular vector
        vector<float> perp(3);
        if (abs(from[0]) < 0.9) {
            perp = {1, 0, 0};
        } else {
            perp = {0, 1, 0};
        }
        vector<float> axis = normalize(crossProduct(from, perp));

        // Create 180-degree rotation matrix around axis
        rotMatrix[0][0] = 2 * axis[0] * axis[0] - 1;
        rotMatrix[0][1] = 2 * axis[0] * axis[1];
        rotMatrix[0][2] = 2 * axis[0] * axis[2];
        rotMatrix[1][0] = 2 * axis[1] * axis[0];
        rotMatrix[1][1] = 2 * axis[1] * axis[1] - 1;
        rotMatrix[1][2] = 2 * axis[1] * axis[2];
        rotMatrix[2][0] = 2 * axis[2] * axis[0];
        rotMatrix[2][1] = 2 * axis[2] * axis[1];
        rotMatrix[2][2] = 2 * axis[2] * axis[2] - 1;

        return rotMatrix;
    }

    // General case: use Rodrigues' rotation formula
    vector<float> axis = normalize(crossProduct(from, to));
    float sinTheta = sqrt(1 - cosTheta * cosTheta);

    // Rodrigues' rotation matrix
    float oneMinusCos = 1 - cosTheta;

    rotMatrix[0][0] = cosTheta + axis[0] * axis[0] * oneMinusCos;
    rotMatrix[0][1] = axis[0] * axis[1] * oneMinusCos - axis[2] * sinTheta;
    rotMatrix[0][2] = axis[0] * axis[2] * oneMinusCos + axis[1] * sinTheta;

    rotMatrix[1][0] = axis[1] * axis[0] * oneMinusCos + axis[2] * sinTheta;
    rotMatrix[1][1] = cosTheta + axis[1] * axis[1] * oneMinusCos;
    rotMatrix[1][2] = axis[1] * axis[2] * oneMinusCos - axis[0] * sinTheta;

    rotMatrix[2][0] = axis[2] * axis[0] * oneMinusCos - axis[1] * sinTheta;
    rotMatrix[2][1] = axis[2] * axis[1] * oneMinusCos + axis[0] * sinTheta;
    rotMatrix[2][2] = cosTheta + axis[2] * axis[2] * oneMinusCos;

    return rotMatrix;
}

// Apply rotation matrix to a point
atom rotateAtom(atom original, vector<vector<float>> rotMatrix, atom pivot) {
    atom rotated = original;

    // Translate to origin (relative to pivot)
    float tempX = original.x - pivot.x;
    float tempY = original.y - pivot.y;
    float tempZ = original.z - pivot.z;

    // Apply rotation
    rotated.x = rotMatrix[0][0] * tempX + rotMatrix[0][1] * tempY + rotMatrix[0][2] * tempZ;
    rotated.y = rotMatrix[1][0] * tempX + rotMatrix[1][1] * tempY + rotMatrix[1][2] * tempZ;
    rotated.z = rotMatrix[2][0] * tempX + rotMatrix[2][1] * tempY + rotMatrix[2][2] * tempZ;

    // Translate back
    rotated.x += pivot.x;
    rotated.y += pivot.y;
    rotated.z += pivot.z;

    return rotated;
}

// Parse PDB file and extract atoms
vector<atom> readPDBFile(string &filename) {
    vector<atom> atoms;
    ifstream file(filename);
    string line;

    if (!file.is_open()) {
        cout << "Error: Could not open file " << filename << endl;
        return atoms;
    }

    while (getline(file, line)) {
        if (line.substr(0, 4) == "ATOM") {
            atom singleAtom;

            // Parse according to PDB format
            singleAtom.atomID = stoi(line.substr(6, 5));
            singleAtom.atomName = line.substr(12, 4);
            singleAtom.residueName = line.substr(17, 3);
            singleAtom.resID = stoi(line.substr(22, 4));
            singleAtom.x = stod(line.substr(30, 8));
            singleAtom.y = stod(line.substr(38, 8));
            singleAtom.z = stod(line.substr(46, 8));
            singleAtom.occupancy = stod(line.substr(54, 6));
            singleAtom.tempFactor = stod(line.substr(60, 6));
            singleAtom.chainID = line.substr(72, 4);

            atoms.push_back(singleAtom);
        }
    }

    file.close();
    cout << "Read " << atoms.size() << " atoms from " << filename << endl;
    return atoms;
}

// Write atoms to PDB file
void writePDBFile(string filename, vector<atom> atoms, string title = "REPOSITIONED LIPID") {
    ofstream file(filename);

    if (!file.is_open()) {
        cout << "Error: Could not create file " << filename << endl;
        return;
    }

    // Write header
    file << "CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1" << endl;
    file << "REMARK Generated by lipid repositioning program" << endl;
    file << "REMARK " << title << endl;

    // Write atoms
    for (atom a : atoms) {
        file << "ATOM  "
             << setw(5) << a.atomID
             << " " << setw(4) << a.atomName
             << " " << setw(3) << a.residueName
             << " " << setw(1) << "1"
             << setw(4) << a.resID
             << "    "
             << setw(8) << fixed << setprecision(3) << a.x
             << setw(8) << fixed << setprecision(3) << a.y
             << setw(8) << fixed << setprecision(3) << a.z
             << setw(6) << fixed << setprecision(2) << a.occupancy
             << setw(6) << fixed << setprecision(2) << a.tempFactor
             << "      " << a.chainID << endl;
    }

    file << "END" << endl;
    file.close();
    cout << "Written " << atoms.size() << " atoms to " << filename << endl;
}

// Find root atom (usually the first nitrogen or phosphorus)
atom findRootAtom(vector<atom> atoms) {
    // Look for nitrogen first (polar head)
    for (atom a : atoms) {
        if (a.atomName.find("N") != string::npos) {
            return a;
        }
    }

    // If no nitrogen, look for phosphorus
    for (atom a : atoms) {
        if (a.atomName.find("P") != string::npos) {
            return a;
        }
    }

    // If neither found, return first atom
    if (!atoms.empty()) {
        return atoms[0];
    }

    // Return empty atom if no atoms
    atom empty = {0, "", "", 0, 0, 0, 0, 0, 0, ""};
    return empty;
}

// Find the atoms that are furthest from the root (tail ends)
vector<atom> findTailEnds(vector<atom> atoms, atom root, int numTails = 2) {
    vector<atom> tailEnds;

    // Calculate distances from root
    vector<pair<float, int>> distances; // Use index instead of atom object
    for (int i = 0; i < atoms.size(); i++) {
        float d = dist(atoms[i], root);
        distances.push_back({d, i});
    }

    // Sort by distance (furthest first) - using custom comparator
    sort(distances.begin(), distances.end(), [](const pair<float, int> &a, const pair<float, int> &b) {
        return a.first > b.first; // Sort by distance in descending order
    });

    // Take the furthest atoms
    for (int i = 0; i < min(numTails, (int)distances.size()); i++) {
        tailEnds.push_back(atoms[distances[i].second]);
    }

    return tailEnds;
}

// Calculate overall lipid direction
vector<float> calculateLipidDirection(vector<atom> atoms, atom root) {
    vector<atom> tailEnds = findTailEnds(atoms, root);

    if (tailEnds.empty()) {
        return {0, 0, 1}; // Default direction
    }

    // Calculate average direction to tail ends
    vector<float> avgDirection = {0, 0, 0};
    for (atom tail : tailEnds) {
        vector<float> dir = direction(root, tail);
        avgDirection[0] += dir[0];
        avgDirection[1] += dir[1];
        avgDirection[2] += dir[2];
    }

    avgDirection[0] /= tailEnds.size();
    avgDirection[1] /= tailEnds.size();
    avgDirection[2] /= tailEnds.size();

    return avgDirection;
}

// Reposition all atoms based on new direction vector
vector<atom> repositionLipid(vector<atom> atoms, atom root, vector<float> originalDirection, vector<float> newDirection) {
    vector<atom> repositionedAtoms;

    // Create rotation matrix
    vector<vector<float>> rotMatrix = createRotationMatrix(originalDirection, newDirection);

    // Apply rotation to each atom
    for (atom a : atoms) {
        atom rotatedAtom = rotateAtom(a, rotMatrix, root);
        repositionedAtoms.push_back(rotatedAtom);
    }

    return repositionedAtoms;
}

int main() {
    string inputFile = "./lipid.pdb";
    string outputFile = "./repositioned_lipid.pdb";

    // Read atoms from PDB file
    vector<atom> originalAtoms = readPDBFile(inputFile);

    if (originalAtoms.empty()) {
        cout << "No atoms found in input file. Exiting." << endl;
        return 1;
    }

    // Find root atom and calculate original direction
    atom root = findRootAtom(originalAtoms);
    vector<float> originalDirection = calculateLipidDirection(originalAtoms, root);

    cout << "Root atom: " << root.atomName << " at ("
         << root.x << ", " << root.y << ", " << root.z << ")" << endl;
    cout << "Original direction: [" << originalDirection[0] << ", "
         << originalDirection[1] << ", " << originalDirection[2] << "]" << endl;

    // Define new direction (example: pointing up in Z direction)
    vector<float> newDirection = {0.0, 0.0, 1.0};

    // Reposition the lipid
    vector<atom> repositionedAtoms = repositionLipid(originalAtoms, root, originalDirection, newDirection);

    // Write repositioned atoms to new PDB file
    writePDBFile(outputFile, repositionedAtoms, "LIPID POINTING UP [0,0,1]");

    // Create additional examples with different directions

    // Example 2: Point in +X direction
    vector<float> xDirection = {1.0, 0.0, 0.0};
    vector<atom> xRepositioned = repositionLipid(originalAtoms, root, originalDirection, xDirection);
    writePDBFile("lipid_x_direction.pdb", xRepositioned, "LIPID POINTING RIGHT [1,0,0]");

    // Example 3: Point in +Y direction
    vector<float> yDirection = {0.0, 1.0, 0.0};
    vector<atom> yRepositioned = repositionLipid(originalAtoms, root, originalDirection, yDirection);
    writePDBFile("lipid_y_direction.pdb", yRepositioned, "LIPID POINTING FORWARD [0,1,0]");

    // Example 4: Point diagonally
    vector<float> diagDirection = {1.0, 1.0, 0.0};
    vector<atom> diagRepositioned = repositionLipid(originalAtoms, root, originalDirection, diagDirection);
    writePDBFile("lipid_diagonal.pdb", diagRepositioned, "LIPID POINTING DIAGONALLY [1,1,0]");

    cout << "\nGenerated PDB files:" << endl;
    cout << "- " << outputFile << " (pointing up)" << endl;
    cout << "- lipid_x_direction.pdb (pointing right)" << endl;
    cout << "- lipid_y_direction.pdb (pointing forward)" << endl;
    cout << "- lipid_diagonal.pdb (pointing diagonally)" << endl;

    return 0;
}