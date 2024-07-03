/*****************************************************************
 *  @file main.cpp                                               *
 *  @brief ���ڹ����������������ϵĲ��ɷ��̵�����ϵͳ��        *
 *                                                               *
 *  �ó����OBJ�ļ���ȡ������Ϣ����������������ʷ֣�            *
 *  �����ʷ���Ϣ�����ʾ���񶥵㲴�ɷ��̵�����ϵͳ��             *
 *  ʹ���ſɱȷ�������ϵͳ�������������д��OBJ�ļ���          *
 *                                                               *   
 *  @author ������ 21312016                                      *
 ****************************************************************/ 

#include "Mesh.h"

/**
 *  @brief ������������֮���ŷ����þ��롣
 *
 *  @param v1 ��һ�����㡣
 *  @param v2 �ڶ������㡣
 *  @return ������������֮��ľ��롣
 */
double get_dist(const Vertex& v1, const Vertex& v2) {
    return sqrt((v1.x - v2.x) * (v1.x - v2.x) + (v1.y - v2.y) * (v1.y - v2.y));
}

/**
 *  @brief �ж�ĳ�������Ƿ�����ĳ�������Ρ�
 *
 *  @param T �����ζ���
 *  @param idx ����������
 *  @return ����������ڸ������Σ��򷵻�true�����򷵻�false��
 */
bool is_vertex(const Triangle& T, int idx) {
    return idx == T.v1 || idx == T.v2 || idx == T.v3;
}

/**
 *  @brief �ж϶����Ƿ�Ϊ�߽�㡣
 *
 *  @param v �������
 *  @return ��������Ǳ߽�㣬�򷵻�true�����򷵻�false��
 */
bool is_boundary_vertex(const Vertex& v) {
    return v.phi > -1e29;
}

/**
 *  @brief �ж��������Ƿ���ȡ�
 *
 *  @param e1 ��һ���ߡ�
 *  @param e2 �ڶ����ߡ�
 *  @return �����������ȣ��򷵻�true�����򷵻�false��
 */
bool is_edge_equal(std::pair<int, int> e1, std::pair<int, int> e2) {
    return ((e1.first == e2.second && e1.second == e2.first) || e1 == e2);
}

/**
 *  @brief �ж�һ�����Ƿ��������ε����Բ�ڡ�
 *
 *  @param p ���жϵĵ㡣
 *  @param v1 �����εĵ�һ�����㡣
 *  @param v2 �����εĵڶ������㡣
 *  @param v3 �����εĵ��������㡣
 *  @return ����������Բ�ڣ��򷵻�true�����򷵻�false��
 */
bool is_in_circumcircle(const Vertex& p, const Vertex& v1, const Vertex& v2, const Vertex& v3) {
    double ax = v1.x - p.x;
    double ay = v1.y - p.y;
    double bx = v2.x - p.x;
    double by = v2.y - p.y;
    double cx = v3.x - p.x;
    double cy = v3.y - p.y;

    double det_ab = ax * by - ay * bx;
    double det_bc = bx * cy - by * cx;
    double det_ca = cx * ay - cy * ax;

    double a_squared = ax * ax + ay * ay;
    double b_squared = bx * bx + by * by;
    double c_squared = cx * cx + cy * cy;

    double det = a_squared * det_bc + b_squared * det_ca + c_squared * det_ab;
    return det > 0;
}

/**
 *  @brief ����һ���������е����ڵĳ������Ρ�
 *
 *  @param vertices �������顣
 *  @param triangles ���������顣
 */
void generate_super_triangle(std::vector<Vertex>& vertices, std::vector<Triangle>& triangles) {
    double min_x = 1e9, min_y = 1e9, max_x = -1e9, max_y = -1e9;

    for (const auto& vertex : vertices) {
        if (vertex.x < min_x) min_x = vertex.x;
        if (vertex.y < min_y) min_y = vertex.y;
        if (vertex.x > max_x) max_x = vertex.x;
        if (vertex.y > max_y) max_y = vertex.y;
    }

    double dx = max_x - min_x;
    double dy = max_y - min_y;
    double mid_x = (min_x + max_x) / 2.0;
    double mid_y = (min_y + max_y) / 2.0;

    Vertex v1 = { mid_x - 2 * dx, mid_y - dy, 0, 0 };
    Vertex v2 = { mid_x + 2 * dx, mid_y - dy, 0, 0 };
    Vertex v3 = { mid_x, mid_y + 2 * dy, 0, 0 };
    vertices.push_back(v1);
    vertices.push_back(v2);
    vertices.push_back(v3);

    Triangle super_tri = { vertices.size() - 3, vertices.size() - 2, vertices.size() - 1 };
    triangles.push_back(super_tri);
}

/**
 *  @brief ��������������ʷ֡�
 *
 *  @param vertices �������顣
 *  @param triangles ���������顣
 */
void construct_delaunay(std::vector<Vertex>& vertices, std::vector<Triangle>& triangles) {
    generate_super_triangle(vertices, triangles);
    for (int i = 0; i < vertices.size() - 3; ++i) {
        std::vector<std::pair<int, int>> edge_buffer;
        for (auto iter = triangles.begin(); iter != triangles.end();) {
            if (is_in_circumcircle(vertices[i], vertices[iter->v1], vertices[iter->v2], vertices[iter->v3])) {
                edge_buffer.push_back({ iter->v1, iter->v2 });
                edge_buffer.push_back({ iter->v2, iter->v3 });
                edge_buffer.push_back({ iter->v3, iter->v1 });
                iter = triangles.erase(iter);   // �����µĵ�����
            }
            else {
                ++iter;
            }
        }
        for (int j = 0; j < edge_buffer.size(); ++j) {
            for (int k = j + 1; k < edge_buffer.size(); ++k) {
                if (is_edge_equal(edge_buffer[j], edge_buffer[k])) {
                    edge_buffer.erase(edge_buffer.begin() + k);
                    edge_buffer.erase(edge_buffer.begin() + j);
                    --j;
                    break;
                }
            }
        }
        for (const auto& edge : edge_buffer) {
            Triangle tri;
            tri.v1 = edge.first;
            tri.v2 = edge.second;
            tri.v3 = i;
            triangles.push_back(tri);
        }
    }
    
    std::vector<int> vertex_removed;
    for (auto iter = triangles.begin(); iter != triangles.end();) {
        if (iter->v1 >= vertices.size() - 3 || iter->v2 >= vertices.size() - 3 || iter->v3 >= vertices.size() - 3 ||
            (is_boundary_vertex(vertices[iter->v1]) && is_boundary_vertex(vertices[iter->v2]) && is_boundary_vertex(vertices[iter->v3]))) {
            vertex_removed.push_back(iter->v1);
            vertex_removed.push_back(iter->v2);
            vertex_removed.push_back(iter->v3);
            iter = triangles.erase(iter);
        }
        else {
            ++iter;
        }
    }
    sort(vertex_removed.begin(), vertex_removed.end());
    vertex_removed.erase(std::unique(vertex_removed.begin(), vertex_removed.end()), vertex_removed.end());

    for (int i = 0; i < vertex_removed.size();) {
        bool flag = false;
        for (auto iter = triangles.begin(); iter != triangles.end(); ++iter) {
            if (vertex_removed[i] == iter->v1 || vertex_removed[i] == iter->v2 || vertex_removed[i] == iter->v3) {
                flag = true;
                break;
            }
        }
        if (flag) {
            vertex_removed.erase(vertex_removed.begin() + i);
        }
        else {
            ++i;
        }
    }

    std::map<int, int> index_map;
    int origin_size = vertices.size();
    for (int i = 0, j = 0, rm_idx = 0; i < origin_size; ++i) {
        if (rm_idx < vertex_removed.size() && i == vertex_removed[rm_idx]) {
            vertices.erase(vertices.begin() + j);
            ++rm_idx;
        }
        else {
            index_map[i] = j;
            ++j;
        }
    }
    for (auto iter = triangles.begin(); iter != triangles.end(); ++iter) {
        iter->v1 = index_map[iter->v1];
        iter->v2 = index_map[iter->v2];
        iter->v3 = index_map[iter->v3];
    }
}

/**
 *  @brief ���체�ɷ��̵�����ϵͳ��
 *
 *  @param vertices �������顣
 *  @param triangles ���������顣
 *  @param mat_data ���ϡ��������з���Ԫ�ء�
 *  @param indices ��¼��ӦԪ����ԭʼ�����е�����Ϣ��
 *  @param ptr ÿ����mat_data�е���ʼλ�á�
 *  @param vector �Ҳ�������
 *  @param n ����������
 */
void construct_equation(const std::vector<Vertex>& vertices, const std::vector<Triangle>& triangles,
    std::vector<double>& mat_data, std::vector<int>& indices, std::vector<int>& ptr,
    std::vector<double>& vector, size_t n) {

    for (int i = 0; i < n; ++i) {
        ptr.push_back(mat_data.size());
        if (is_boundary_vertex(vertices[i])) {
            mat_data.push_back(1.0);
            indices.push_back(i);
            vector[i] = vertices[i].phi;
            continue;
        }

        std::vector<int> triangle_indices;
        for (int j = 0; j < triangles.size(); ++j) {
            if (is_vertex(triangles[j], i)) {
                triangle_indices.push_back(j);
            }
        }

        double area = 0;
        double* row_data = new double[n];
        for (int j = 0; j < n; ++j) {
            row_data[j] = 0;
        }

        for (int j = 0; j < triangle_indices.size(); ++j) {
            const Triangle& T = triangles[triangle_indices[j]];
            int v1 = T.v1, v2 = T.v2;
            if (v1 == i)
                v1 = T.v3;
            else if (v2 == i)
                v2 = T.v3;
            double a = get_dist(vertices[v1], vertices[v2]);
            double b = get_dist(vertices[i], vertices[v2]);
            double c = get_dist(vertices[i], vertices[v1]);
            double cosA = (b * b + c * c - a * a) / (2 * b * c);
            double sinA = sqrt(1 - cosA * cosA);
            double R = a / (2 * sinA);
            double hb = sqrt(R * R - b * b / 4);
            double hc = sqrt(R * R - c * c / 4);
            area += (hb * b + hc * c) / 4;
            row_data[v1] += hb / b;
            row_data[v2] += hc / c;
            row_data[i] -= (hb / b + hc / c);
        }
        vector[i] = area * vertices[i].f;
        for (int j = 0; j < n; ++j) {
            if (row_data[j] != 0) {
                mat_data.push_back(row_data[j]);
                indices.push_back(j);
            }
        }

        delete[] row_data;
    }
    ptr.push_back(mat_data.size());
}

/**
 *  @brief ʹ���ſɱȵ������������ϵͳ��
 *
 *  @param mat_data ���ϡ��������з���Ԫ�ء�
 *  @param indices ��¼��ӦԪ����ԭʼ�����е�����Ϣ��
 *  @param ptr ÿ����mat_data�е���ʼλ�á�
 *  @param vector �Ҳ�������
 *  @param x ��������
 *  @param n �����С��
 *  @param maxIter ������������
 *  @param tol ����ֹͣ���ݲ
 */
void jacobi_method(std::vector<double>& mat_data, std::vector<int>& indices, std::vector<int>& ptr,
    std::vector<double>& vector, std::vector<double>& x, int n, int maxIter, double tol) {
    std::vector<double> x_old(n, 0.0);

    for (int iter = 0; iter < maxIter; ++iter) {
        for (int i = 0; i < n; ++i) {
            x_old[i] = x[i];
        }
        for (int i = 0; i < n; ++i) {
            double sigma = 0.0, div;
            for (int j = ptr[i]; j < ptr[i + 1]; ++j) {
                if (indices[j] == i) {
                    div = mat_data[j];
                }
                else {
                    sigma += mat_data[j] * x_old[indices[j]];
                }
            }
            x[i] = (vector[i] - sigma) / div;
        }

        double norm = 0.0;
        for (int i = 0; i < n; ++i) {
            norm += (x[i] - x_old[i]) * (x[i] - x_old[i]);
        }
        norm = std::sqrt(norm);
        if (norm < tol) {
            break;
        }
    }
}

/**
 *  @brief ����⺯������������ϵͳ����⡣
 *
 *  @param vertices �������顣
 *  @param triangles ���������顣
 */
void Solution(std::vector<Vertex>& vertices, std::vector<Triangle>& triangles)
{
    size_t n = vertices.size();
    std::vector<double> mat_data;
    std::vector<int> indices, ptr;
    std::vector<double> vector(n, 0.0);
    std::vector<double> result(n, 0.0);

    construct_delaunay(vertices, triangles);
    construct_equation(vertices, triangles, mat_data, indices, ptr, vector, n);
    jacobi_method(mat_data, indices, ptr, vector, result, n, 1000, 1e-6);

    for (int i = 0; i < n; ++i) {
        vertices[i].phi = result[i];
    }
}

int main() {
    Mesh mesh;

    mesh.loadFromOBJ("C:/Users/86153/Desktop/��ʹ��ͷ�����еȼ۽���/��������Ĳ�����ֵ����/��ĩ����ҵ/NumericalMesh/models/map_vertices.obj");

    Solution(mesh.vertices, mesh.triangles);

    mesh.WriteToOBJ("C:/Users/86153/Desktop/��ʹ��ͷ�����еȼ۽���/��������Ĳ�����ֵ����/��ĩ����ҵ/NumericalMesh/models/map_solved_new.obj");

    return 0;
}