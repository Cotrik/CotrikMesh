/*
 * IsotropicRemesh.cpp
 *
 *  Created on: Dec 29, 2017
 *      Author: cotrik
 */

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/OFF_reader.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <boost/function_output_iterator.hpp>
#include <fstream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor edge_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;
struct halfedge2edge {
    halfedge2edge(const Mesh& m, std::vector<edge_descriptor>& edges)
            : m_mesh(m), m_edges(edges) {
    }
    void operator()(const halfedge_descriptor& h) const {
        m_edges.push_back(edge(h, m_mesh));
    }
    const Mesh& m_mesh;
    std::vector<edge_descriptor>& m_edges;
};

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: IsotropicRemesh input.off output.off\n";
        return -1;
    }
    const char* filename = argv[1];
    const char* outfilename = argv[2];
    // polygon soup : reorientation
    {
        std::ifstream input(filename);
        std::vector<K::Point_3> points;
        std::vector<std::vector<std::size_t> > polygons;
        if (!CGAL::read_OFF(input, points, polygons)) {
            std::cerr << "Error parsing the OFF file " << std::endl;
            return 1;
        }
        CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);
        Polyhedron mesh;
        CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, mesh);
        if (CGAL::is_closed(mesh) && (!CGAL::Polygon_mesh_processing::is_outward_oriented(mesh)))
            CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);
        std::ofstream out("oriented.off");
        out << mesh;
        out.close();
    }
    std::ifstream input("oriented.off");
    Mesh mesh;
    if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
        std::cerr << "Not a valid input file." << std::endl;
        return 1;
    }
    double target_edge_length = 0.04;
    unsigned int nb_iter = 3;
    std::cout << "Split border...";
    std::vector<edge_descriptor> border;
    PMP::border_halfedges(faces(mesh), mesh, boost::make_function_output_iterator(halfedge2edge(mesh, border)));
    PMP::split_long_edges(border, target_edge_length, mesh);
    std::cout << "done." << std::endl;
    std::cout << "Start remeshing of " << filename << " (" << num_faces(mesh) << " faces)..." << std::endl;
    PMP::isotropic_remeshing(faces(mesh), target_edge_length, mesh, PMP::parameters::number_of_iterations(nb_iter).protect_constraints(true));
    std::cout << "Remeshing done." << std::endl;

    std::ofstream output(outfilename);
    output << mesh;
    return 0;
}

