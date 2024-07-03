#include "Mesh.h"

Mesh::Mesh(const std::string& filename){
    loadFromOBJ(filename);
}

void Mesh::loadFromOBJ(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    std::string line;
    while (getline(file, line)) {
        std::stringstream ss(line);
        std::string type;
        ss >> type;
        if (type == "v") {
            Vertex vertex;
            ss >> vertex.x >> vertex.y >> vertex.phi >> vertex.f;
            vertices.push_back(vertex);
        }
        else if (type == "f") {
            Triangle t;
            ss >> t.v1 >> t.v2 >> t.v3;
            // Adjust indices to be zero-based if they are not
            t.v1--; t.v2--; t.v3--;
            triangles.push_back(t);
        }
    }

    file.close();
}

void Mesh::WriteToOBJ(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }

    // Write vertices
    for (const auto& vertex : vertices) {
        file << "v " << vertex.x << " " << vertex.y << " " << vertex.phi << " " << vertex.f << std::endl;
    }

    // Write triangles
    for (const auto& triangle : triangles) {
        // OBJ file indices are 1-based
        file << "f " << triangle.v1 + 1 << " " << triangle.v2 + 1 << " " << triangle.v3 + 1 << std::endl;
    }

    file.close();
}
