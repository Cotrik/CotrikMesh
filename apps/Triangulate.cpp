/*
 * Triangulate.cpp
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
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <boost/foreach.hpp>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Surface_mesh;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: Triangulate input.off output.off\n";
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
    Surface_mesh mesh;
    if (!input || !(input >> mesh) || mesh.is_empty()) {
        std::cerr << "Not a valid off file." << std::endl;
        return 1;
    }
    CGAL::Polygon_mesh_processing::triangulate_faces(mesh);
    // Confirm that all faces are triangles.
    BOOST_FOREACH(boost::graph_traits<Surface_mesh>::face_descriptor fit, faces(mesh))
        if (next(next(halfedge(fit, mesh), mesh), mesh) != prev(halfedge(fit, mesh), mesh)) std::cerr << "Error: non-triangular face left in mesh." << std::endl;
    std::ofstream cube_off(outfilename);
    cube_off << mesh;
    return 0;
}
