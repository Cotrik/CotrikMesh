/*
 * HexGen.cpp
 *
 *  Created on: May 25, 2017
 *      Author: cotrik
 */


//


#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Ray_3 Ray;
typedef CGAL::Surface_mesh<Point> CGALMesh;
typedef boost::graph_traits<CGALMesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<CGALMesh>::halfedge_descriptor halfedge_descriptor;
typedef CGAL::AABB_face_graph_triangle_primitive<CGALMesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;


#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "ArgumentManager.h"
#include <iostream>
#include <fstream>

bool IsPointInsideMesh(const Point& p, const Tree& tree)
{
    Vector v(0, 0, 1);
    Ray ray(p, v);
    return tree.number_of_intersected_primitives(ray) % 2 != 0;
}


void GetBoundingBox(const Mesh&mesh, glm::dvec3&min_coordinate, glm::dvec3&max_coordinate);
void GenerateCubeMeshInBoundingBox(const glm::dvec3& min_coordinate, const glm::dvec3& max_coordinate, const size_t numberOfGridsPerAxis, Mesh& cubeMesh);

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cout << "Usage: HexGen input.vtk output.vtk size=<50>" << std::endl;
        return -1;
    }

    ArgumentManager argumentManager(argc, argv);
    std::string input_filename = argv[1];
    std::string output_filename = argv[2];
    size_t size = 50;
    std::string strSize = argumentManager.get("size");
    if (!strSize.empty()) size = std::stoi(strSize);

    std::cout << "---------------------------------------" << std::endl;
    std::cout << "input_filename = " << input_filename << std::endl;
    std::cout << "output_filename = " << output_filename << std::endl;
    std::cout << "size = " << size << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    MeshFileReader reader(input_filename.c_str());
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    const char* triFilename = "surface.tri.off";
    MeshFileWriter surfaceWriter(mesh, triFilename);
    surfaceWriter.WriteSurfaceOff();

    std::ifstream input(triFilename);
    CGALMesh cgalMesh;
    input >> cgalMesh;
    Tree tree(faces(cgalMesh).first, faces(cgalMesh).second, cgalMesh);

    glm::dvec3 min_coordinate(1000000.0f, 1000000.0f, 1000000.0f);
    glm::dvec3 max_coordinate(-1000000.0f, -1000000.0f, -1000000.0f);
    GetBoundingBox(mesh, min_coordinate, max_coordinate);

    Mesh hexMesh;
    GenerateCubeMeshInBoundingBox(min_coordinate, max_coordinate, size, hexMesh);

//    MeshFileWriter writer(hexMesh, output_filename.c_str());
//    writer.WriteFile();

    std::vector<Cell> C;
    C.reserve(hexMesh.C.size());
//#pragma omp parallel for
    for (size_t i = 0; i < hexMesh.C.size(); i++) {
        const Cell& cell = hexMesh.C.at(i);
        glm::dvec3 p0 = hexMesh.V.at(cell.Vids.at(0)).xyz();
        glm::dvec3 p6 = hexMesh.V.at(cell.Vids.at(6)).xyz();
        glm::dvec3 mid = (p0 + p6) * 0.5f;
        Point p(mid.x, mid.y, mid.z);
        //if (mesh.IsPointInside(p)) C.push_back(cell);
        if (IsPointInsideMesh(p, tree)) C.push_back(cell);
    }
    MeshFileWriter writer(hexMesh.V, C, output_filename.c_str());
    writer.WriteFile();

    return 0;
}

void GetBoundingBox(const Mesh& mesh, glm::dvec3& min_coordinate, glm::dvec3& max_coordinate)
{
    for (auto const & v : mesh.V) {
        //if (!v.isBoundary) continue;

        if (v.x < min_coordinate.x) min_coordinate.x = v.x;
        if (v.y < min_coordinate.y) min_coordinate.y = v.y;
        if (v.z < min_coordinate.z) min_coordinate.z = v.z;

        if (v.x > max_coordinate.x) max_coordinate.x = v.x;
        if (v.y > max_coordinate.y) max_coordinate.y = v.y;
        if (v.z > max_coordinate.z) max_coordinate.z = v.z;
    }
}

void GenerateCubeMeshInBoundingBox(const glm::dvec3& min_coordinate, const glm::dvec3& max_coordinate, const size_t numberOfGridsPerAxis/* = 50*/, Mesh& cubeMesh)
{
    double xLength = max_coordinate.x - min_coordinate.x;
    double yLength = max_coordinate.y - min_coordinate.y;
    double zLength = max_coordinate.z - min_coordinate.z;

    double maxLengh = xLength;
    if (yLength > maxLengh) maxLengh = yLength;
    if (zLength > maxLengh) maxLengh = zLength;

    double cubeLength = maxLengh/numberOfGridsPerAxis;

    size_t xTotalStep = round(xLength / cubeLength) + 1;
    size_t yTotalStep = round(yLength / cubeLength) + 1;
    size_t zTotalStep = round(zLength / cubeLength) + 1;
    size_t yzPlaneNumber = yTotalStep * zTotalStep;

  std::cout << "*************************************************" << std::endl;
  std::cout << " xLength: " << xLength << std::endl;
  std::cout << " yLength: " << yLength << std::endl;
  std::cout << " zLength: " << zLength << std::endl;
  std::cout << "*************************************************" << std::endl;
  std::cout << " xTotalStep: " << xTotalStep << std::endl;
  std::cout << " yTotalStep: " << yTotalStep << std::endl;
  std::cout << " zTotalStep: " << zTotalStep << std::endl;

  std::cout << " cubeLength: " << cubeLength << std::endl;
  std::cout << "*************************************************" << std::endl;

    cubeMesh.V.clear();
    cubeMesh.C.clear();
    cubeMesh.V.resize(xTotalStep * yTotalStep * zTotalStep);
    cubeMesh.C.resize((xTotalStep - 1) * (yTotalStep - 1) * (zTotalStep - 1));
    double xPos = min_coordinate.x;
    double yPos = min_coordinate.y;
    double zPos = min_coordinate.z;
    size_t cid = 0;
    //size_t vid = 0;
    std::vector<size_t> v8(8, 0);
//#pragma omp parallel for
    for (size_t xStep = 0; xStep < xTotalStep; xStep++) {
        xPos = min_coordinate.x + cubeLength * xStep;
        for (size_t yStep = 0; yStep < yTotalStep; yStep++) {
            yPos = min_coordinate.y + cubeLength * yStep;
            for (size_t zStep = 0; zStep < zTotalStep; zStep++) {
                zPos = min_coordinate.z + cubeLength * zStep;
                glm::dvec3 v(xPos, yPos, zPos);
                size_t vid = zStep + zTotalStep * yStep + yzPlaneNumber * xStep;
                cubeMesh.V[vid] = v;
                cubeMesh.V[vid].id = vid;
                if (xStep > 0 && yStep > 0 && zStep > 0) {
                    //std::vector<size_t> v8(8, 0);
                    v8[0] = vid - yzPlaneNumber;
                    v8[1] = vid - yzPlaneNumber - zTotalStep;
                    v8[2] = vid - yzPlaneNumber - zTotalStep - 1;
                    v8[3] = vid - yzPlaneNumber - 1;
                    v8[4] = vid;
                    v8[5] = vid - zTotalStep;
                    v8[6] = vid - zTotalStep - 1;
                    v8[7] = vid - 1;
                    //size_t cid = zStep - 1 + (zTotalStep - 1) * (yStep - 1) + (yTotalStep - 1) * (zTotalStep - 1) * (xStep - 1);
                    cubeMesh.C[cid].Vids = v8;
                    cubeMesh.C[cid].id = cid++;
                }
                //vid++;
            }
        }
    }
}

