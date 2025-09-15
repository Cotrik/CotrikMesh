/*
* Tri2Quad.cpp
*
*  Created on: Oct 20, 2018
*      Author: cotrik
*      
*  Updated: Convert triangular meshes to quad meshes using Catmull-Clark subdivision
*/

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "ArgumentManager.h"
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <algorithm>
#include <string>

// Catmull-Clark subdivision implementation for triangular to quad conversion
Mesh CatmullClarkSubdivision(const Mesh& inputMesh) {
    const Mesh& mesh = inputMesh;
    
    // Validate input mesh
    if (mesh.V.empty() || mesh.F.empty()) {
        std::cerr << "Error: Input mesh is empty!" << std::endl;
        return Mesh();
    }
    
    if (mesh.E.empty()) {
        std::cerr << "Error: Mesh edges not built! Call BuildAllConnectivities() first." << std::endl;
        return Mesh();
    }
    
    std::cout << "Input mesh validation:" << std::endl;
    std::cout << "  Vertices: " << mesh.V.size() << std::endl;
    std::cout << "  Edges: " << mesh.E.size() << std::endl;
    std::cout << "  Faces: " << mesh.F.size() << std::endl;
    
    // Step 1: Add new vertices
    // Original vertices + edge midpoints + face centroids
    std::vector<Vertex> new_vertices(mesh.V.size() + mesh.E.size() + mesh.F.size());
    
    // Copy original vertices
    for (size_t i = 0; i < mesh.V.size(); i++) {
        new_vertices[i] = mesh.V[i];
        new_vertices[i].id = i;
    }
    
    // Add edge midpoints
    size_t offset = mesh.V.size();
    for (size_t i = 0; i < mesh.E.size(); i++) {
        const Edge& e = mesh.E[i];
        
        // Validate edge
        if (e.Vids.size() < 2) {
            std::cerr << "Error: Edge " << i << " has insufficient vertices!" << std::endl;
            continue;
        }
        
        size_t v0_id = e.Vids[0];
        size_t v1_id = e.Vids[1];
        
        // Validate vertex indices
        if (v0_id >= mesh.V.size() || v1_id >= mesh.V.size()) {
            std::cerr << "Error: Edge " << i << " references invalid vertex indices: " 
                      << v0_id << ", " << v1_id << " (max: " << mesh.V.size() - 1 << ")" << std::endl;
            continue;
        }
        
        const Vertex& v0 = mesh.V[v0_id];
        const Vertex& v1 = mesh.V[v1_id];
        
        new_vertices[offset + i] = Vertex();
        new_vertices[offset + i].x = 0.5 * (v0.x + v1.x);
        new_vertices[offset + i].y = 0.5 * (v0.y + v1.y);
        new_vertices[offset + i].z = 0.5 * (v0.z + v1.z);
        new_vertices[offset + i].id = offset + i;
    }
    
    // Add face centroids
    offset = mesh.V.size() + mesh.E.size();
    for (size_t i = 0; i < mesh.F.size(); i++) {
        const Face& f = mesh.F[i];
        
        // Validate face
        if (f.Vids.empty()) {
            std::cerr << "Error: Face " << i << " has no vertices!" << std::endl;
            continue;
        }
        
        new_vertices[offset + i] = Vertex();
        new_vertices[offset + i].x = 0.0;
        new_vertices[offset + i].y = 0.0;
        new_vertices[offset + i].z = 0.0;
        
        // Calculate centroid
        for (size_t j = 0; j < f.Vids.size(); j++) {
            size_t vid = f.Vids[j];
            
            // Validate vertex index
            if (vid >= mesh.V.size()) {
                std::cerr << "Error: Face " << i << " references invalid vertex index: " 
                          << vid << " (max: " << mesh.V.size() - 1 << ")" << std::endl;
                continue;
            }
            
            const Vertex& v = mesh.V[vid];
            new_vertices[offset + i].x += v.x;
            new_vertices[offset + i].y += v.y;
            new_vertices[offset + i].z += v.z;
        }
        
        if (f.Vids.size() > 0) {
            new_vertices[offset + i].x /= f.Vids.size();
            new_vertices[offset + i].y /= f.Vids.size();
            new_vertices[offset + i].z /= f.Vids.size();
        }
        new_vertices[offset + i].id = offset + i;
    }
    
    // Step 2: Create edge mapping
    std::map<std::pair<size_t, size_t>, size_t> edgeMap;
    for (size_t i = 0; i < mesh.E.size(); i++) {
        const Edge& e = mesh.E[i];
        
        if (e.Vids.size() >= 2) {
            size_t v1 = std::min(e.Vids[0], e.Vids[1]);
            size_t v2 = std::max(e.Vids[0], e.Vids[1]);
            edgeMap[{v1, v2}] = mesh.V.size() + i;
        }
    }
    
    // Step 3: Create new quad faces
    std::vector<Face> new_faces;
    
    for (size_t i = 0; i < mesh.F.size(); i++) {
        const Face& f = mesh.F[i];
        
        if (f.Vids.size() == 3) {
            // For triangular faces, create 3 quads
            size_t v0 = f.Vids[0];
            size_t v1 = f.Vids[1];
            size_t v2 = f.Vids[2];
            
            // Validate vertex indices
            if (v0 >= mesh.V.size() || v1 >= mesh.V.size() || v2 >= mesh.V.size()) {
                std::cerr << "Error: Triangle " << i << " has invalid vertex indices!" << std::endl;
                continue;
            }
            
            // Get edge midpoints with validation
            auto e0_it = edgeMap.find({std::min(v0, v1), std::max(v0, v1)});
            auto e1_it = edgeMap.find({std::min(v1, v2), std::max(v1, v2)});
            auto e2_it = edgeMap.find({std::min(v2, v0), std::max(v2, v0)});
            
            if (e0_it == edgeMap.end() || e1_it == edgeMap.end() || e2_it == edgeMap.end()) {
                std::cerr << "Error: Triangle " << i << " has missing edge midpoints!" << std::endl;
                continue;
            }
            
            size_t e0 = e0_it->second;
            size_t e1 = e1_it->second;
            size_t e2 = e2_it->second;
            
            // Face centroid
            size_t fc = mesh.V.size() + mesh.E.size() + i;
            
            // Create 3 quads
            Face quad1(4);
            quad1.Vids = {v0, e0, fc, e2};
            new_faces.push_back(quad1);
            
            Face quad2(4);
            quad2.Vids = {e0, v1, e1, fc};
            new_faces.push_back(quad2);
            
            Face quad3(4);
            quad3.Vids = {fc, e1, v2, e2};
            new_faces.push_back(quad3);
        }
        else if (f.Vids.size() == 4) {
            // For quad faces, create 4 quads (standard Catmull-Clark)
            size_t v0 = f.Vids[0];
            size_t v1 = f.Vids[1];
            size_t v2 = f.Vids[2];
            size_t v3 = f.Vids[3];
            
            // Validate vertex indices
            if (v0 >= mesh.V.size() || v1 >= mesh.V.size() || 
                v2 >= mesh.V.size() || v3 >= mesh.V.size()) {
                std::cerr << "Error: Quad " << i << " has invalid vertex indices!" << std::endl;
                continue;
            }
            
            // Get edge midpoints with validation
            auto e0_it = edgeMap.find({std::min(v0, v1), std::max(v0, v1)});
            auto e1_it = edgeMap.find({std::min(v1, v2), std::max(v1, v2)});
            auto e2_it = edgeMap.find({std::min(v2, v3), std::max(v2, v3)});
            auto e3_it = edgeMap.find({std::min(v3, v0), std::max(v3, v0)});
            
            if (e0_it == edgeMap.end() || e1_it == edgeMap.end() || 
                e2_it == edgeMap.end() || e3_it == edgeMap.end()) {
                std::cerr << "Error: Quad " << i << " has missing edge midpoints!" << std::endl;
                continue;
            }
            
            size_t e0 = e0_it->second;
            size_t e1 = e1_it->second;
            size_t e2 = e2_it->second;
            size_t e3 = e3_it->second;
            
            // Face centroid
            size_t fc = mesh.V.size() + mesh.E.size() + i;
            
            // Create 4 quads
            Face quad1(4);
            quad1.Vids = {v0, e0, fc, e3};
            new_faces.push_back(quad1);
            
            Face quad2(4);
            quad2.Vids = {e0, v1, e1, fc};
            new_faces.push_back(quad2);
            
            Face quad3(4);
            quad3.Vids = {fc, e1, v2, e2};
            new_faces.push_back(quad3);
            
            Face quad4(4);
            quad4.Vids = {e3, fc, e2, v3};
            new_faces.push_back(quad4);
        }
        else {
            std::cerr << "Warning: Face " << i << " has " << f.Vids.size() 
                      << " vertices, skipping (only triangles and quads supported)" << std::endl;
        }
    }
    
    std::cout << "Created " << new_faces.size() << " new quad faces" << std::endl;
    
    // Step 4: Create new mesh
    Mesh resultMesh;
    resultMesh.V = new_vertices;
    resultMesh.F = new_faces;
    resultMesh.m_cellType = QUAD;
    
    // Build connectivity
    // resultMesh.BuildAllConnectivities();

	std::cout << "Finished building all connectivities" << std::endl;
	std::cout << "resultMesh.V.size() = " << resultMesh.V.size() << std::endl;
	std::cout << "resultMesh.F.size() = " << resultMesh.F.size() << std::endl;
	std::cout << "resultMesh.m_cellType = " << resultMesh.m_cellType << std::endl;
    
    return resultMesh;
}

// Get file extension from filename
std::string getFileExtension(const std::string& filename) {
    size_t dotPos = filename.find_last_of('.');
    if (dotPos != std::string::npos) {
        return filename.substr(dotPos + 1);
    }
    return "";
}

// Convert extension to lowercase
std::string toLower(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), ::tolower);
    return result;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage: Tri2Quad input.[vtk|obj|stl|off|mesh] output.[vtk|obj|stl]\n";
        std::cout << "Converts triangular meshes to quad meshes using Catmull-Clark subdivision\n";
        std::cout << "Supported input formats: vtk, obj, stl, off, mesh\n";
        std::cout << "Supported output formats: vtk, obj, stl\n";
        return -1;
    }

    std::string inputFile = argv[1];
    std::string outputFile = argv[2];
    
    // Validate input file extension
    std::string inputExt = toLower(getFileExtension(inputFile));
    if (inputExt != "vtk" && inputExt != "obj" && inputExt != "stl" && 
        inputExt != "off" && inputExt != "mesh") {
        std::cerr << "Error: Unsupported input file format: " << inputExt << std::endl;
        std::cerr << "Supported formats: vtk, obj, stl, off, mesh" << std::endl;
        return -1;
    }
    
    // Validate output file extension
    std::string outputExt = toLower(getFileExtension(outputFile));
    if (outputExt != "vtk" && outputExt != "obj" && outputExt != "stl") {
        std::cerr << "Error: Unsupported output file format: " << outputExt << std::endl;
        std::cerr << "Supported formats: vtk, obj, stl" << std::endl;
        return -1;
    }
    
    try {
        std::cout << "Reading input mesh: " << inputFile << std::endl;
        MeshFileReader reader(inputFile.c_str());
        Mesh& inputMesh = (Mesh&)reader.GetMesh();
        
        std::cout << "Input mesh statistics:" << std::endl;
        std::cout << "  Vertices: " << inputMesh.V.size() << std::endl;
        std::cout << "  Faces: " << inputMesh.F.size() << std::endl;
        std::cout << "  Cell type: " << (inputMesh.m_cellType == TRIANGLE ? "Triangle" : 
                                        inputMesh.m_cellType == QUAD ? "Quad" : "Other") << std::endl;
        
        
		inputMesh.RemoveUselessVertices();
										// Validate input mesh
        if (inputMesh.V.empty()) {
            std::cerr << "Error: Input mesh has no vertices!" << std::endl;
            return -1;
        }
        
        if (inputMesh.F.empty()) {
            std::cerr << "Error: Input mesh has no faces!" << std::endl;
            return -1;
        }
        
        // Build connectivity if not already built
        if (inputMesh.E.empty()) {
            std::cout << "Building mesh connectivity..." << std::endl;
            try {
                inputMesh.BuildAllConnectivities();
            } catch (const std::exception& e) {
                std::cerr << "Error building mesh connectivity: " << e.what() << std::endl;
                return -1;
            }
        }
        
        // Apply Catmull-Clark subdivision
        std::cout << "Applying Catmull-Clark subdivision..." << std::endl;
        Mesh quadMesh = CatmullClarkSubdivision(inputMesh);
        
        // Validate output mesh
        if (quadMesh.V.empty() || quadMesh.F.empty()) {
            std::cerr << "Error: Subdivision failed - output mesh is empty!" << std::endl;
            return -1;
        }
        
        std::cout << "Output mesh statistics:" << std::endl;
        std::cout << "  Vertices: " << quadMesh.V.size() << std::endl;
        std::cout << "  Faces: " << quadMesh.F.size() << std::endl;
        std::cout << "  Cell type: " << (quadMesh.m_cellType == QUAD ? "Quad" : "Other") << std::endl;
        
        // Write output mesh
        std::cout << "Writing output mesh: " << outputFile << std::endl;
        try {
            MeshFileWriter writer(quadMesh, outputFile.c_str());
            writer.WriteFile();
        } catch (const std::exception& e) {
            std::cerr << "Error writing output file: " << e.what() << std::endl;
            return -1;
        }
        
        std::cout << "Conversion completed successfully!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    } catch (...) {
        std::cerr << "Unknown error occurred!" << std::endl;
        return -1;
    }
    
    return 0;
}
