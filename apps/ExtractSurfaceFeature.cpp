/*
 * ExtractPatches.cpp
 *
 *  Created on: May 17, 2017
 *      Author: cotrik
 */
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "FrameField.h"
#include "PolyLine.h"
#include "FrameOpt.h"
#include "FeatureLine.h"
#include "MeshQuality.h"
#include "Patches.h"
#include "ArgumentManager.h"
#include "Util.h"

const double PI = 3.1415926535;

void GetFaceTypes(const Mesh& mesh, std::vector<FACE_TYPE>& faceTypes) {
	faceTypes.resize(mesh.F.size());
	for (unsigned long i = 0; i < mesh.F.size(); i++) {
		const Face& face = mesh.F.at(i);
		faceTypes.at(i) = Util::GetFaceType(face.normal);
	}
}


void WriteSharpEdgesVtk(const char* filename, const Mesh& m_mesh, const std::vector<FeatureLine>& featureLines) {
	const std::vector<Vertex>& V = m_mesh.V;
	const std::vector<Edge>& E = m_mesh.E;

	std::ofstream ofs(filename);
	ofs << "# vtk DataFile Version 2.0" << std::endl
		<< filename << std::endl
		<< "ASCII" << std::endl << std::endl
		<< "DATASET UNSTRUCTURED_GRID" << std::endl;
	ofs << "POINTS " << V.size() << " float" << std::endl;
	for (size_t i = 0; i < V.size(); i++)
		ofs << V.at(i).x << " " << V.at(i).y << " " << V.at(i).z << std::endl;
	size_t numOfSharpVertices = 0;
	for (size_t i = 0; i < featureLines.size(); i++) {
		const FeatureLine& fl = featureLines.at(i);
		numOfSharpVertices += fl.Vids.size();
	}

	ofs << "CELLS " << featureLines.size() << " " << numOfSharpVertices + featureLines.size() << std::endl;
	for (size_t i = 0; i < featureLines.size(); i++) {
		const FeatureLine& fl = featureLines.at(i);
		ofs << fl.Vids.size();
		for (size_t j = 0; j < fl.Vids.size(); j++) {
			const size_t vid = fl.Vids.at(j);
			ofs << " " << vid;
		}
		ofs << std::endl;
	}

	ofs << "CELL_TYPES " << featureLines.size() << std::endl;
	for (size_t i = 0; i < featureLines.size(); i++) {
		ofs << 4 << std::endl;
	}

	ofs << "CELL_DATA " << featureLines.size() << std::endl
		<< "SCALARS " << " Feature" << " int 1" << std::endl
		<< "LOOKUP_TABLE default" << std::endl;
	for (size_t i = 0; i < featureLines.size(); i++) {
		ofs << i << std::endl;
	}
}

int main(int argc, char* argv[]) {
	if (argc < 3) {
		std::cout << "Usage: ExtractSurfaceFeature <input_tri_file> <output_tri_vtk_file> min_numOfFaces_per_patch=<10> localangle=<170> globalangle=<90>" << std::endl;
		return -1;
	}
	ArgumentManager argumentManager(argc, argv);
	std::string filename = argv[1];
	std::string output_filename = argv[2];
	int min_numOfFaces_per_patch = 10;
	double localangle = 170.0;
	double globalangle = 90.0;

	double coslocalangle = cos((180.0 - localangle) * PI / 180.0);
	const std::string strLocalAngle = argumentManager.get("localangle");
	if (!strLocalAngle.empty()) coslocalangle = cos((180.0 - std::stod(strLocalAngle)) * PI / 180.0);

	double cosglobalangle = cos((180.0 - globalangle) * PI / 180.0);
	const std::string strGlobalaAngle = argumentManager.get("globalangle");
	if (!strGlobalaAngle.empty()) cosglobalangle = cos((180.0 - std::stod(strGlobalaAngle)) * PI / 180.0);

	const std::string strmin_numOfFaces_per_patch = argumentManager.get("min_numOfFaces_per_patch");
	if (!strmin_numOfFaces_per_patch.empty()) min_numOfFaces_per_patch = std::stoi(strmin_numOfFaces_per_patch);

	std::cout << "---------------------------------------" << std::endl;
	std::cout << "filename = " << filename << std::endl;
	std::cout << "output_filename = " << output_filename << std::endl;
	std::cout << "min_numOfFaces_per_patch = " << min_numOfFaces_per_patch << std::endl;
	std::cout << "localangle = " << strLocalAngle << std::endl;
	std::cout << "globalangle = " << strGlobalaAngle << std::endl;
	std::cout << "---------------------------------------" << std::endl;

	MeshFileReader reader(filename.c_str());
	Mesh& mesh = (Mesh&)reader.GetMesh();
	mesh.RemoveUselessVertices();
	mesh.BuildAllConnectivities();
	mesh.ExtractBoundary();
	//    mesh.ExtractLayers();
	mesh.ExtractSingularities();
	mesh.SetCosAngleThreshold(coslocalangle);
	//    mesh.LabelSurface();
	mesh.GetNormalOfSurfaceFaces();
	mesh.GetNormalOfSurfaceVertices();

	Patches patches(mesh);
	patches.SetGlobalCosAngle(cosglobalangle);
	if (!strGlobalaAngle.empty()) patches.Extract2();
	else patches.Extract();
	{	
		MeshFileWriter sharpEdgesFileWriter(mesh, "SharpEdges.vtk");
		sharpEdgesFileWriter.WriteSharpEdgesVtk();

		for (auto& v : mesh.V) {
			int count = 0;
			for (auto eid : v.N_Eids)
				if (mesh.E.at(eid).isSharpFeature) ++count;
			if (count == 1) {
				v.type = CORNER;
				v.isCorner = true;
			}
		}
		std::vector<size_t> cornerIds;
		for (const auto& v : mesh.V)
			if (v.isCorner) cornerIds.push_back(v.id);
		MeshFileWriter writer(mesh, "Corners.vtk");
		writer.WriteVerticesVtk(cornerIds);
	}

	std::vector<size_t> faceids;
	for (auto& f : mesh.F)
		if (f.isBoundary) faceids.push_back(f.id);
	MeshFileWriter facesFileWriter(mesh, output_filename.c_str());
	facesFileWriter.WriteFacesVtk(faceids);
	//std::vector<bool> copy(mesh.E.size());
	//auto E_copy = mesh.E;
	//for (auto& e : mesh.E) { 
	//	copy[e.id] = e.isSharpFeature;
	//	e.isSharpFeature = false; 
	//}
	mesh.LabelSharpEdges(true);
	// for (auto& e : mesh.E) e.isSharpFeature = copy[e.id];
	std::vector<FeatureLine> featureLines(mesh.numOfSharpEdges, FeatureLine(mesh));
	for (size_t i = 0; i < mesh.numOfSharpEdges; i++)
		featureLines.at(i).Extract(i);
	WriteSharpEdgesVtk("FeatureLines.vtk", mesh, featureLines);
}




