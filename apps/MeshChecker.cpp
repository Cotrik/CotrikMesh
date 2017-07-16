/*
 * MeshChecker.cpp
 *
 *  Created on: Mar 11, 2017
 *      Author: cotrik
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include <iostream>
#include "ArgumentManager.h"

#include <string.h>
#include "frame/io.h"
#include "frame/frame_of_mesh.h"

char pathT[300] = "tri.off";
char path_IOH[300] = "hex.off";

//char Choices[300]="VOX";
char Choices[300] = "PO";
char Hex_NUM[300] = "8648";

h_io hio;
//int main(int argc, char* argv[])
//{
//    char pathioh[300] = "./";
//    sprintf(path_IOH, "%s", argv[1]);
//
//    printf("reading hex mesh...\n");
//    char path_temp[300];
//    //sprintf(path_temp, "%s%s", path_IOH, ".off");
//    hio.read_hex_mesh_off(hex_mesh.HVs, hex_mesh.HHs, path_IOH);
//    printf("constructing Es connectivities...\n");
//    construct_Es(hex_mesh.HVs, hex_mesh.HEs, hex_mesh.HHs);
//    printf("constructing Fs connectivities...\n");
//    construct_Fs(hex_mesh.HVs, hex_mesh.HEs, hex_mesh.HFs, hex_mesh.HHs);
//    determine_boundary_info(hex_mesh.HVs, hex_mesh.HEs, hex_mesh.HFs, hex_mesh.HHs);
//    hex_mesh.average_e_len = average_len(hex_mesh.HVs, hex_mesh.HEs);
//
//    //extract base-complex
//    frame_of_mesh fom;
//    printf("base-complex extraction...\n");
//    fom.base_complex_extraction();
//    hio.save_base_complex_VTK("base_complex.vtk");
//    hio.save_BasecomplexNode_VTK("base_complex_nodes.vtk");
//    hio.save_BasecomplexEdge_VTK("base_complex_edges.vtk");
//    hio.save_singularG_VTK("singularity.vtk");
//    hio.save_singularG_VTK_Color_Nodes("singularityNodesColor.vtk");
//    hio.save_singularG_VTK_Color_Edges("singularityEdgesColor.vtk");
//    hio.save_singularNode_VTK("singularityNodes.vtk");
//    hio.save_component_VTK("components.vtk");
//    component_class com_class;
//    com_class.extract_all_component_sheets();
//    //hio.save_sheets_VTK("sheets.vtk");
//    for (int i = 0; i < all_com_sheets.size(); i++)
//    {
//        char sheetname[300] = {0};
//        sprintf(sheetname, "%s%d%s", "sheet_", i, ".hex.vtk");
//        hio.save_sheet_hex_VTK(sheetname, i);
//        sprintf(sheetname, "%s%d%s", "sheet_", i, ".geo.vtk");
//        hio.save_sheet_geo_VTK(sheetname, i);
//    }
//    return 0;
//}
int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: MeshChecker input_tri_file\n";
        return -1;
    }

    MeshFileReader reader(argv[1]);
    //reader.GetPointsScalarFields();
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();

    int sum = 0;
    if (mesh.m_cellType == QUAD)
        for (size_t i = 0; i < mesh.V.size(); i++) {
            const Vertex& v = mesh.V.at(i);
            const size_t valence = v.N_Vids.size();
            if (v.isBoundary && valence != 4)
                sum += 4 - valence;
        }
    else if (mesh.m_cellType == HEXAHEDRA) {
        char pathioh[300] = "./";
        sprintf(path_IOH, "%s", argv[1]);

        printf("reading hex mesh...\n");
        char path_temp[300];
        //sprintf(path_temp, "%s%s", path_IOH, ".off");
        hio.read_hex_mesh_off(hex_mesh.HVs, hex_mesh.HHs, path_IOH);
        printf("constructing Es connectivities...\n");
        construct_Es(hex_mesh.HVs, hex_mesh.HEs, hex_mesh.HHs);
        printf("constructing Fs connectivities...\n");
        construct_Fs(hex_mesh.HVs, hex_mesh.HEs, hex_mesh.HFs, hex_mesh.HHs);
        determine_boundary_info(hex_mesh.HVs, hex_mesh.HEs, hex_mesh.HFs, hex_mesh.HHs);
        hex_mesh.average_e_len = average_len(hex_mesh.HVs, hex_mesh.HEs);

        //extract base-complex
        frame_of_mesh fom;
        printf("base-complex extraction...\n");
        fom.base_complex_extraction();
        hio.save_base_complex_VTK("base_complex.vtk");
        hio.save_BasecomplexNode_VTK("base_complex_nodes.vtk");
        hio.save_BasecomplexEdge_VTK("base_complex_edges.vtk");
        hio.save_singularG_VTK("singularity.vtk");
        hio.save_singularG_VTK_Color_Nodes("singularityNodesColor.vtk");
        hio.save_singularG_VTK_Color_Edges("singularityEdgesColor.vtk");
        hio.save_singularNode_VTK("singularityNodes.vtk");
        hio.save_component_VTK("components.vtk");
        for (size_t i = 0; i < mesh.E.size(); i++) {
            const Edge& e = mesh.E.at(i);
            const size_t valence = e.N_Cids.size();
            if (valence != 4)
                sum += 4 - valence;
        }
    }

    cout << "singularities sum = " << sum << std::endl;

    return 0;
}
