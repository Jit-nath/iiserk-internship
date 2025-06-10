#pragma once
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <vector>

#include "atom.h"
using namespace std;

struct node {
    int resID;
    int child_resID;
    double bond_length;
    double bond_angle;
};

double calculateDistance(const atom &a1, const atom &a2) {
    double dx = a1.pos.x - a2.pos.x;
    double dy = a1.pos.y - a2.pos.y;
    double dz = a1.pos.z - a2.pos.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

double calculateAngle(const atom &a1, const atom &middle, const atom &a2) {
    // Vector from middle to a1
    double v1x = a1.pos.x - middle.pos.x;
    double v1y = a1.pos.y - middle.pos.y;
    double v1z = a1.pos.z - middle.pos.z;

    // Vector from middle to a2
    double v2x = a2.pos.x - middle.pos.x;
    double v2y = a2.pos.y - middle.pos.y;
    double v2z = a2.pos.z - middle.pos.z;

    // Calculate dot product
    double dot = v1x * v2x + v1y * v2y + v1z * v2z;

    // Calculate magnitudes
    double mag1 = sqrt(v1x * v1x + v1y * v1y + v1z * v1z);
    double mag2 = sqrt(v2x * v2x + v2y * v2y + v2z * v2z);

    // Calculate angle in radians, then convert to degrees
    double cos_angle = dot / (mag1 * mag2);
    // Clamp to avoid numerical errors
    cos_angle = max(-1.0, min(1.0, cos_angle));
    double angle_rad = acos(cos_angle);
    return angle_rad * 180.0 / M_PI;
}

atom findAtom(const vector<atom> &molecular_structure, int resID) {
    for (const auto &a : molecular_structure) {
        if (a.resID == resID) {
            return a;
        }
    }
    return atom();
}

void dfs(int current,
         int parent,
         const map<int, vector<int>> &adj_list,
         vector<int> &path,
         vector<vector<int>> &all_paths,
         vector<bool> &visited) {

    visited[current] = true;
    path.push_back(current);

    bool hasUnvisitedChild = false;
    if (adj_list.find(current) != adj_list.end()) {
        for (int neighbor : adj_list.at(current)) {
            if (neighbor != parent && !visited[neighbor]) {
                hasUnvisitedChild = true;
                dfs(neighbor, current, adj_list, path, all_paths, visited);
            }
        }
    }

    // If this is a leaf node (no unvisited children), save the path
    if (!hasUnvisitedChild) {
        all_paths.push_back(path);
    }

    path.pop_back();
    visited[current] = false;
}

vector<node> createGraph(map<int,
                             vector<int>> adjacency_list,
                         vector<atom> molecular_structure,
                         atom root) {
    vector<node> nodes;

    int maxResID = 0;
    for (const auto &pair : adjacency_list) {
        maxResID = max(maxResID, pair.first);
        for (int neighbor : pair.second) {
            maxResID = max(maxResID, neighbor);
        }
    }

    // Find all paths from root to leaves using DFS
    vector<vector<int>> all_paths;
    vector<int> current_path;
    vector<bool> visited(maxResID + 1, false);

    dfs(root.resID, -1, adjacency_list, current_path, all_paths, visited);

    // Sort paths by length to get the deepest branches
    sort(all_paths.begin(), all_paths.end(), [](const vector<int> &a, const vector<int> &b) {
        return a.size() > b.size();
    });

    // Process each connection in the adjacency list
    for (const auto &pair : adjacency_list) {
        int parentID = pair.first;
        atom parentAtom = findAtom(molecular_structure, parentID);

        for (int childID : pair.second) {
            // Avoid duplicate edges (only process each edge once)
            if (parentID < childID) {
                atom childAtom = findAtom(molecular_structure, childID);

                // Calculate bond length
                double bondLength = calculateDistance(parentAtom, childAtom);

                // Calculate bond angle
                double bondAngle = 0.0;

                // Find a third atom to calculate angle
                // Look for another neighbor of the parent atom
                for (int thirdID : pair.second) {
                    if (thirdID != childID) {
                        atom thirdAtom = findAtom(molecular_structure, thirdID);
                        bondAngle = calculateAngle(childAtom, parentAtom, thirdAtom);
                        break;
                    }
                }

                // If no third atom found in parent's neighbors, look in child's
                // neighbors
                if (bondAngle == 0.0 &&
                    adjacency_list.find(childID) != adjacency_list.end()) {
                    for (int thirdID : adjacency_list.at(childID)) {
                        if (thirdID != parentID) {
                            atom thirdAtom = findAtom(molecular_structure, thirdID);
                            bondAngle = calculateAngle(parentAtom, childAtom, thirdAtom);
                            break;
                        }
                    }
                }

                // Create node for parent->child relationship
                node newNode;
                newNode.resID = parentID;
                newNode.child_resID = childID;
                newNode.bond_length = bondLength;
                newNode.bond_angle = bondAngle;
                nodes.push_back(newNode);

                // Create node for child->parent relationship (if needed for
                // bidirectional)
                node reverseNode;
                reverseNode.resID = childID;
                reverseNode.child_resID = parentID;
                reverseNode.bond_length = bondLength;
                reverseNode.bond_angle = bondAngle;
                nodes.push_back(reverseNode);
            }
        }
    }

    return nodes;
}

void printGraph(const vector<node> &graph) {
    cout << left << setw(8) << "Parent" << setw(8) << "Child" << setw(15)
         << "Bond Length" << setw(12) << "Bond Angle" << endl;
    cout << left << setw(8) << "------" << setw(8) << "-----" << setw(15)
         << "-----------" << setw(12) << "----------" << endl;
    for (const auto &n : graph) {
        cout << left << setw(8) << n.resID << setw(8) << n.child_resID << setw(15)
             << fixed << setprecision(4) << n.bond_length << setw(12) << fixed
             << setprecision(2) << n.bond_angle << endl;
    }
}