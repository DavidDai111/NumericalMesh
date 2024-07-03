#pragma once
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>

struct Vertex {
    float x, y; // Vertex coordinates
    float phi;  // phi value
    float f;    // Laplacian operator value at this vertex
};

struct Triangle {
    int v1, v2, v3; // Indices of vertices forming a triangle
};

struct Edge {
    int v1, v2; // Vertex indices

    // Ensure consistent ordering for comparison purposes
    Edge(int v1_, int v2_) : v1(std::min(v1_, v2_)), v2(std::max(v1_, v2_)) {}

    bool operator<(const Edge& other) const {
        return v1 < other.v1 || (v1 == other.v1 && v2 < other.v2);
    }
};

class Mesh {
public:
    std::vector<Vertex> vertices;
    std::vector<Triangle> triangles;
    
    Mesh() {}
    Mesh(const std::string& filename);
    ~Mesh() {}
    void loadFromOBJ(const std::string& filename);
    void WriteToOBJ(const std::string& filename) const;

private:
    
};
