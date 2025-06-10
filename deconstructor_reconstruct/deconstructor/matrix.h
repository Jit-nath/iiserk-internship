#include "atom.h"
#include <cmath>
#include <vector>
using namespace std;

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
