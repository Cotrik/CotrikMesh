#include "KDTree.h"
#include "ParallelFor.h"

KDTree::KDTree(const Mesh& mesh_) {
    mesh = Mesh(mesh_);
    std::vector<Point> points;
    points.resize(mesh.V.size());
    PARALLEL_FOR_BEGIN(0, mesh.V.size()) {
        auto& v = mesh.V[i];
        glm::dvec3 n(0.0, 0.0, 0.0);
        for (auto fid: v.N_Fids) {
            auto& f = mesh.F.at(fid);
            int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
            auto& AB = mesh.V.at(f.Vids.at((idx+1)%f.Vids.size())).xyz() - v.xyz();
            auto& AC = mesh.V.at(f.Vids.at((idx+2)%f.Vids.size())).xyz() - v.xyz();
            auto& AD = mesh.V.at(f.Vids.at((idx+3)%f.Vids.size())).xyz() - v.xyz();
            auto T_area_a = glm::cross(AB, AC);
            auto T_area_b = glm::cross(AC, AD);
            n += (0.5 * glm::length(T_area_a) * T_area_a);
            n += (0.5 * glm::length(T_area_b) * T_area_b); 
        }
        n = glm::normalize(n);
        points[i] = Point(i, v.xyz(), n, (v.isBoundary || v.type == FEATURE));
    } PARALLEL_FOR_END();
    root = BuildKDTree(points, 0);
}

Node* KDTree::BuildKDTree(std::vector<Point> points, int depth) {
    if (points.size() == 0) return NULL;

    int axis = depth % 3;
    int medianIndex = points.size() / 2;

    std::nth_element(points.begin(), points.begin() + medianIndex, points.end(), [&](const Point& a, const Point& b) {
        return Compare(a.p, b.p, axis);
    });

    Node* node = new Node(points[medianIndex]);
    node->left = BuildKDTree(std::vector<Point>(points.begin(), points.begin() + medianIndex), depth + 1);
    node->right = BuildKDTree(std::vector<Point>(points.begin() + medianIndex + 1, points.end()), depth + 1);

    return node;
}

bool KDTree::Compare(const glm::dvec3& a, const glm::dvec3& b, int axis) {
    return a[axis] < b[axis];
}

double KDTree::Distance(const Point& a, const Point& b) {
    auto x = a.p-b.p;
    return glm::length(x);
    // return glm::length(x) * (1.0-glm::dot(glm::normalize(a.n), glm::normalize(b.n)));
    // return std::sqrt(glm::dot(x,x) + (1.0-glm::dot(glm::normalize(a.n), glm::normalize(b.n))));
}

Point KDTree::SearchNearestPoint(Node* node, const Point& p, Node* best, double bestDist, int depth, std::vector<Point>& points) {
    if (node == NULL) return Point(-1, glm::dvec3(0.0, 0.0, 0.0), glm::dvec3(0.0, 0.0, 0.0), false);

    int axis = depth % 3;
    int goLeft = Compare(p.p, node->p.p, axis);

    Node* next = goLeft ? node->left : node->right;
    Node* secondary = goLeft ? node->right : node->left;

    Point nearest = SearchNearestPoint(next, p, best, bestDist, depth + 1, points);
    double dist_nearest = Distance(p, nearest);
    double dist_p = Distance(p, node->p);
    // if (dist_nearest > dist_p && normal_strength_nearest < normal_strength_p) {
    if (dist_nearest > dist_p) {
        nearest = node->p;
        std::rotate(points.rbegin(), points.rbegin() + 1, points.rend());
        points[0] = node->p;
    }

    if (fabs(p.p[axis]-node->p.p[axis]) < Distance(p, nearest)) {
        Point nearest2 = SearchNearestPoint(secondary, p, best, bestDist, depth + 1, points);
        double dist_nearest2 = Distance(p, nearest2);
        if (dist_nearest2 < dist_nearest) {
            nearest = nearest2;
            std::rotate(points.rbegin(), points.rbegin() + 1, points.rend());
            points[0] = nearest2;
        }
    }

    return nearest;

    /*double dist = Distance(node->p, p);
    if (dist < bestDist) {
        bestDist = dist;
        best = node;
    }

    if (Compare(p.p, node->p.p, axis)) {
        SearchNearestPoint(node->left, p, best, bestDist, depth + 1);
        if (bestDist > fabs(p.p[axis] - node->p.p[axis])) {
            SearchNearestPoint(node->right, p, best, bestDist, depth + 1);
        }
    } else {
        SearchNearestPoint(node->right, p, best, bestDist, depth + 1);
        if (bestDist > fabs(p.p[axis] - node->p.p[axis])) {
            SearchNearestPoint(node->left, p, best, bestDist, depth + 1);
        }
    }*/
}

void KDTree::SearchAndProject(glm::dvec3& p, glm::dvec3& n) {
    Node* best = NULL;
    double bestDist = std::numeric_limits<double>::max();
    // SearchNearestPoint(root, Point(-1, p, n), best, bestDist, 0);
    std::vector<Point> points(3);
    Point nearest = SearchNearestPoint(root, Point(-1, p, n, false), best, bestDist, 0, points);
    // auto nearest_point = points[0];
    // if (nearest_point.id == -1) return;
    // for (int i = 1; i < points.size(); i++) {
    //     auto point = points[i];
    //     if (point.id == -1) continue;

    //     if (glm::dot(glm::normalize(point.n), glm::normalize(n)) > glm::dot(glm::normalize(nearest_point.n), glm::normalize(n))) {
    //         nearest_point = point;
    //     }
    // }
    if (nearest.id == -1) return;
    // if (nearest_node == NULL) return;
    // if (bestDist == std::numeric_limits<double>::max()) return;
    // p = nearest.p;
    // p = best->p.p;
    // return;
    // int vid = nearest_point.id;
    int vid = nearest.id;
    // int vid = best->p.id;
    auto& v = mesh.V.at(vid);
    if (v.isBoundary || v.type == FEATURE) return;
    // double min_distance = std::numeric_limits<double>::max();
    glm::dvec3 intersection(0.0, 0.0, 0.0);
    for (auto fid: v.N_Fids) {
        auto& f = mesh.F.at(fid);
        int idx = std::distance(f.Vids.begin(), std::find(f.Vids.begin(), f.Vids.end(), v.id));
        auto projection = [&] (Vertex& v1, Vertex& v2, Vertex& v3) {
            glm::dvec3 AB = v2.xyz() - v1.xyz();
            glm::dvec3 AC = v3.xyz() - v1.xyz();
            glm::dvec3 BC = v3.xyz() - v2.xyz();
            glm::dvec3 CA = v1.xyz() - v3.xyz();

            glm::dvec3 n = glm::cross(AB, AC);
            // double t = glm::dot(glm::normalize(n), v1.xyz()-p);
            // glm::dvec3 projected_point = p + (t * n);
            double t = ((n.x*v1.x) - (n.x*p.x) + (n.y*v1.y - n.y*p.y) + (n.z*v1.z - n.z*p.z)) / ((n.x*n.x) + (n.y*n.y) + (n.z*n.z));
            glm::dvec3 projected_point = glm::dvec3(p.x + t*n.x, p.y + t*n.y, p.z + t*n.z);

            // glm::dvec3 dV = p - v1.xyz();
            // glm::dvec3 projected_point = v1.xyz() + (dV - ((glm::dot(dV, n) / glm::length(n)) * n));
            
            
            // glm::dvec3 intersect(0.0, 0.0, 0.0);
            // auto plane = Plane(v1, v2, v3);
            // double distance = plane.DistanseFromPoint(p, intersect);

            glm::dvec3 AP = projected_point - v1.xyz();
            glm::dvec3 BP = projected_point - v2.xyz();

            double T_area = 0.5 * glm::length(n);
            if (0.5 * (glm::length(glm::cross(AB, AP)) + glm::length(glm::cross(AC, AP)), glm::length(glm::cross(BP, BC))) > T_area) return;
            intersection = projected_point;
        };
        projection(v, mesh.V.at(f.Vids.at((idx+1)%f.Vids.size())), mesh.V.at(f.Vids.at((idx+2)%f.Vids.size())));
        projection(v, mesh.V.at(f.Vids.at((idx+2)%f.Vids.size())), mesh.V.at(f.Vids.at((idx+3)%f.Vids.size())));
    }
    if (glm::length(intersection) > 0.0) {
        p = intersection;
    }
}