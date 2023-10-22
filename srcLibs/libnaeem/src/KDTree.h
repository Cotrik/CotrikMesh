#ifndef KDTREE_H_
#define KDTREE_H_

#include "Mesh.h"
#include <algorithm>
#include <chrono>
#include <thread>
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_segment_primitive.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_triangle_primitive.h>

class SurfaceProjector {
    public:
        // typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
        typedef CGAL::Simple_cartesian<double> K;
        typedef K::Point_3 Point;
        typedef K::Vector_3 Vector;
        typedef K::Segment_3 Segment;
        typedef K::Triangle_3 Triangle;

        typedef std::list<Triangle>::iterator Triangle_Iterator;
        typedef CGAL::AABB_triangle_primitive<K, Triangle_Iterator> Triangle_Primitive;
        typedef CGAL::AABB_traits<K, Triangle_Primitive> Triangle_AABB_Traits;
        typedef CGAL::AABB_tree<Triangle_AABB_Traits> Triangle_AABB_Tree;

        typedef std::list<Segment>::iterator Edge_Iterator;
        typedef CGAL::AABB_segment_primitive<K, Edge_Iterator> Edge_Primitive;
        typedef CGAL::AABB_traits<K, Edge_Primitive> Edge_AABB_Traits;
        typedef CGAL::AABB_tree<Edge_AABB_Traits> Edge_AABB_Tree;
        
        // std::unique_ptr<AABB_tree> main_tree;
        // std::unique_ptr<Sharp_Edge_AABBTree> sharp_edge_tree;
        std::list<Triangle> triangles;
        Triangle_AABB_Tree main_tree;
        std::list<Segment> edges;
        Edge_AABB_Tree sharp_edge_tree;
        SurfaceProjector(const Mesh& mesh) {
            // for (int i = 0; i < mesh.F.size(); i++) {
            //     auto& f = mesh.F[i];
            //     for (int j = 0; j < f.Vids.size(); j++) {
            //         auto& v1 = mesh.V[f.Vids[j]];
            //         auto& v2 = mesh.V[f.Vids[(j + 1) % f.Vids.size()]];
            //         Point a(v1.x, v1.y, v1.z);
            //         Point b(v2.x, v2.y, v2.z);
            //         edges.push_back(Segment(a, b));
            //     }    
            // }
            // Point a(1.0, 0.0, 0.0);
            // Point b(0.0, 1.0, 0.0);
            // Point c(0.0, 0.0, 1.0);
            // Point d(0.0, 0.0, 0.0);

            // edges.push_back(Segment(a, b));
            // edges.push_back(Segment(b, c));
            // edges.push_back(Segment(c, d));
            // sharp_edge_tree = std::make_unique<Sharp_Edge_AABBTree>(edges.begin(), edges.end());
            // sharp_edge_tree.insert(edges.begin(), edges.end());
            // sharp_edge_tree.accelerate_distance_queries();
            // Sharp_Edge_AABBTree tree(edges.begin(), edges.end());

            // Point query(0.1, 0.5, 0.5);
            // Point closest_point = sharp_edge_tree.closest_point(query);
            // std::cout << "closest point: " << closest_point << std::endl;
            // std::vector<Triangle_3> triangles;
            std::set<std::pair<int, int>> sharp_edges;
            for (int i = 0; i < mesh.F.size(); i++) {
                std::vector<Point> points(4);
                auto& f = mesh.F[i];
                // std::cout << f.Vids.size() << " ";
                for (int j = 0; j < f.Vids.size(); j++) {
                    // Vertex& v1 = mesh.V[f.Vids[i]];
                    Vertex v1 = mesh.V.at(f.Vids.at(j));
                    Vertex v2 = mesh.V[f.Vids[(j + 1) % f.Vids.size()]];
                    points[j] = Point(v1.x, v1.y, v1.z);
                    std::pair<int, int> edge = std::make_pair(std::min(v1.id, v2.id), std::max(v1.id, v2.id));
                    if ((v1.isBoundary || v1.type == FEATURE) && (v2.isBoundary || v2.type == FEATURE) && sharp_edges.find(edge) == sharp_edges.end()) {
                        // std::cout << "v1: " << edge.first << " v2: " << edge.second << std::endl;
                        sharp_edges.insert(edge);
                    }
                }
                // std::cout << "*************************" << std::endl;
                // for (auto point: points) {
                //     std::cout << point << std::endl;
                // }
                // std::cout << "-------------------------" << std::endl;
                Triangle t1(points[0], points[1], points[2]);
                Triangle t2(points[0], points[2], points[3]);
                triangles.push_back(t1);
                triangles.push_back(t2);
                // triangles.push_back(Triangle(points[0], points[1], points[3]));
                // triangles.push_back(Triangle(points[1], points[2], points[3]));
            }
            // std::cout << std::endl;
            // Point a(4.60846, 15.166, 2e-06);
            // Point b(4.52087, 15.2188, 2e-06);
            // Point c(4.41227, 15.1825, 2e-06);
            // Point d(4.40176, 15.1379, 2e-06);
            // triangles.push_back(Triangle(a, b, c));
            // triangles.push_back(Triangle(a, c, d));
            
            // main_tree = std::make_unique<AABB_tree>();
            main_tree.insert(triangles.begin(), triangles.end());
            main_tree.accelerate_distance_queries();

            // CGAL::Cartesian_converter<Kernel, K> converter;
            // std::vector<Segment_3> edges;
            for (auto& e: sharp_edges) {
                Vertex v1 = mesh.V[e.first];
                Vertex v2 = mesh.V[e.second];
                Point a(v1.x, v1.y, v1.z);
                Point b(v2.x, v2.y, v2.z);
                Segment s(a, b);
                edges.push_back(s);
                // std::cout << "v1: " << v1.id << " v2: " << v2.id << std::endl;
                // Kernel::Segment_3 s(Kernel::Point_3(v1.x, v1.y, v1.z), Kernel::Point_3(v2.x, v2.y, v2.z));
                // edges.push_back(converter(s));
            }
            // sharp_edge_tree = std::make_unique<Sharp_Edge_AABBTree>();
            // std::cout << "edges: " << edges.size() << std::endl;
            // edges.clear();
            // Kernel::Segment_3 s(Kernel::Point_3(0, 0, 0), Kernel::Point_3(1, 0, 0));
            // edges = {converter(s)};
            sharp_edge_tree.insert(edges.begin(), edges.end());
            sharp_edge_tree.accelerate_distance_queries();
        }

        glm::dvec3 projectToBoundary(const glm::dvec3& p_) {
            // CGAL::Cartesian_converter<Kernel, K> converter;
            // CGAL::Cartesian_converter<K, Kernel> inverser;
            // Kernel::Point_3 p(0.5, 0, 0);
            Point query(p_.x, p_.y, p_.z);
            // std::cout << "query: " << query << std::endl;
            Point closest_point = sharp_edge_tree.closest_point(query);
            // std::cout << "closest point: " << closest_point << std::endl;
            // glm::dvec3 result(closest_point.x(), closest_point.y(), closest_point.z());
            // std::cout << "result: " << result.x << " " << result.y << " " << result.z << std::endl;
            // std::cout << "p: " << p.x() << " " << p.y() << " " << p.z() << std::endl;
            // Kernel::Point_3 closest_point = inverser(sharp_edge_tree->closest_point(converter(p)));
            // std::cout << "closest_point: " << closest_point.x() << " " << closest_point.y() << " " << closest_point.z() << std::endl;
            // return glm::dvec3(closest_point.x(), closest_point.y(), closest_point.z()) - p_;
            return glm::dvec3(0.0);
        }

        glm::dvec3 projectToSurface(const glm::dvec3& p_, const glm::dvec3& dir_) {
            // typedef AABB_tree::Intersection_and_primitive_id<Segment>::Type segment_intersection;
            // Point_3 p(p_.x, p_.y, p_.z);
            // Vector_3 dir(dir_.x, dir_.y, dir_.z);
            // Segment_3 segment(p-dir, p+dir);
            glm::dvec3 result(0.0);
            typedef Triangle_AABB_Tree::Intersection_and_primitive_id<Segment>::Type Seg_intersection;
            Point p(p_.x, p_.y, p_.z);
            // std::cout << "query: " << p << std::endl;
            Vector dir(dir_.x, dir_.y, dir_.z);
            // Vector dir(0.0, 0.0, 0.0);
            // if (glm::distance(p_-dir_,p_+dir_) == 0.) return result;
            Segment segment(p-dir, p+dir);
            std::list<Seg_intersection> intersections;
            // std::list<std::pair<Point, Triangle>> intersections;
            // auto closest_point = main_tree.closest_point(p);
            // std::list<segment_intersection> intersections;
            main_tree.all_intersections(segment, std::back_inserter(intersections));
            double min_dist = std::numeric_limits<double>::max();
            dir = dir / std::sqrt(dir.squared_length());
            // std::cout << "intersections: " << intersections.size() << std::endl;
            for (auto& intersection: intersections) {
                // auto intersection_point = boost::get<Point>(&(intersection.first));
                // const Point& pt = *intersection_point;
                // auto pt = intersection.first;
                const Point *pt = boost::get<Point>(&intersection.first);
                // std::cout << "intersection: " << pt << std::endl;
                auto id = intersection.second;
                Triangle triangle = *id;
                Vector normal = CGAL::normal(triangle.vertex(0), triangle.vertex(1), triangle.vertex(2));
                normal = normal / std::sqrt(normal.squared_length());
                if (normal * dir > 0.707) {
                    glm::dvec3 pt_(pt->x(), pt->y(), pt->z());
                    double dist = glm::length(pt_ - p_);
                    // double dist = std::sqrt(CGAL::squared_distance(p, pt));
                    if (dist < min_dist) {
                        min_dist = dist;
                        result = pt_ - p_;
                    }
                } 
            }
            // std::cout << "intersection: " << result.x << " " << result.y << " " << result.z << std::endl;
            return result;
        }
};


struct Point {
    int id = -1;
    glm::dvec3 p;
    glm::dvec3 n;
    bool feature;
    Point() {}
    Point(int id_, const glm::dvec3& p_, const glm::dvec3& n_, bool feature_) : id(id_), p(p_), n(n_), feature(feature_) {}
};

struct Node {
    Point p;
    Node* left;
    Node* right;
    Node(const Point& p_) : p(p_), left(NULL), right(NULL) {}
};

class KDTree {
    public:
        KDTree(const Mesh& mesh_);

        Node* BuildKDTree(std::vector<Point> points, int depth);
        bool Compare(const glm::dvec3& a, const glm::dvec3& b, int axis);
        double Distance(const Point& a, const Point& b);
        Point SearchNearestPoint(Node* node, const Point& p, Node* best, double bestDist, int depth, std::vector<Point>& points);
        void SearchAndProject(glm::dvec3& p, glm::dvec3& n);

        Mesh mesh;
        Node* root;
};

#endif