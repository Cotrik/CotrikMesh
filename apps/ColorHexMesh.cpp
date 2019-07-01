/*
 * ColorHexMesh.cpp
 *
 *  Created on: Sep 26, 2018
 *      Author: davim
 */

#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "EdgeLines.h"
#include "BaseComplex.h"
#include "BaseComplexQuad.h"
#include "BaseComplexSheet.h"
#include "BaseComplexSheetQuad.h"
#include "BaseComplexChord.h"
#include "BaseComplexEditor.h"
#include "SingularityGraph.h"
#include "ArgumentManager.h"

#include <iostream>

void ColorEdge(Mesh& mesh);
void ColorFace(Mesh& mesh, BaseComplexChord& baseComplexChords);
void ColorFace(BaseComplex& baseComplex);

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: ColorHexMesh hex.vtk" << std::endl;
        return -1;
    }
    ArgumentManager argumentManager(argc, argv);
    std::string filename = argv[1];
    std::cout << "---------------------------------------" << std::endl;
    std::cout << "filename = " << filename << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    MeshFileReader reader(filename.c_str());
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    mesh.ExtractSingularities();
    //"Info: cos10 = 0.984807753; cos15 = 0.965925826; cos20 = 0.939692621; cos25 = 0.906307787; cos30 = 0.866025404\n\n";
    const double cosangle = 0.866025404;
    mesh.SetCosAngleThreshold(cosangle);
    mesh.LabelSurface();
    mesh.LabelSharpEdges(true);
    // For extracting singularity Graph
    mesh.BuildParallelE();
    mesh.BuildConsecutiveE();
    mesh.BuildOrthogonalE();

    BaseComplex baseComplex(mesh);
    baseComplex.Build();

    BaseComplexChord baseComplexChords(baseComplex);
    baseComplexChords.Extract();

    //ColorFace(mesh, baseComplexChords);
    ColorFace(baseComplex);
    return 0;
}

void ColorParallelEdges(Mesh& mesh, const size_t edgeid, std::vector<std::set<int>>& edge_colors) {

}

void ColorEdge(Mesh& mesh) {
    std::vector<std::set<int>> edge_colors;
    std::set<size_t> conflict_eids;
    for (auto& c : mesh.C) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 4; ++j) {

            }
    }
}

void ColorFace(Mesh& mesh, BaseComplexChord& baseComplexChords) {
    std::vector<int> componentFaceColors(baseComplexChords.baseComplex.componentF.size(), -1);
    for (int i = 0; i < baseComplexChords.chords_componentFaceIds.size(); ++i) {
        auto& chord_componentFaceIds = baseComplexChords.chords_componentFaceIds.at(i);
        std::set<int> colors;
        for (auto chord_componentFaceId : chord_componentFaceIds)
            colors.insert(componentFaceColors[chord_componentFaceId]);
        if (colors.size() == 1) {
            if (*colors.begin() == -1) {
                std::set<int> otherCompoentFacecolors;
                auto& chord_componentCellIds = baseComplexChords.chords_componentCellIds.at(i);
                for (auto chord_componentCellId : chord_componentCellIds) {
                    auto& componentCell = baseComplexChords.baseComplex.componentC.at(chord_componentCellId);
                    for (auto componentFaceId : componentCell.Fids)
                        otherCompoentFacecolors.insert(componentFaceColors.at(componentFaceId));
                }
                if (otherCompoentFacecolors.size() == 1) {
                    for (auto chord_componentFaceId : chord_componentFaceIds)
                        componentFaceColors.at(chord_componentFaceId) = 0;
                } else {
                    std::vector<int> sortedColors(otherCompoentFacecolors.begin(), otherCompoentFacecolors.end());
                    if (sortedColors == std::vector<int>{-1, 0}) {
                        for (auto chord_componentFaceId : chord_componentFaceIds)
                            componentFaceColors.at(chord_componentFaceId) = 1;
                    } else if (sortedColors == std::vector<int>{-1, 1}) {
                        for (auto chord_componentFaceId : chord_componentFaceIds)
                            componentFaceColors.at(chord_componentFaceId) = 2;
                    } else if (sortedColors == std::vector<int>{-1, 2}) {
                        for (auto chord_componentFaceId : chord_componentFaceIds)
                            componentFaceColors.at(chord_componentFaceId) = 0;
                    } else if (sortedColors == std::vector<int>{-1, 0, 1}) {
                        for (auto chord_componentFaceId : chord_componentFaceIds)
                            componentFaceColors.at(chord_componentFaceId) = 2;
                    } else if (sortedColors == std::vector<int>{-1, 0, 2}) {
                        for (auto chord_componentFaceId : chord_componentFaceIds)
                            componentFaceColors.at(chord_componentFaceId) = 1;
                    } else if (sortedColors == std::vector<int>{-1, 1, 2}) {
                        for (auto chord_componentFaceId : chord_componentFaceIds)
                            componentFaceColors.at(chord_componentFaceId) = 0;
                    }
                }
            }
        } else if (colors.size() == 2 && (*colors.begin() == -1 && *std::next(colors.begin()) == -1)) {
            int color = *colors.begin() == -1 ? *std::next(colors.begin()) : *colors.begin();
            for (auto chord_componentFaceId : chord_componentFaceIds)
                componentFaceColors.at(chord_componentFaceId) = color;
        } else if (colors.size() == 3) {
            std::set<int> otherCompoentFacecolors;
            auto& chord_componentCellIds = baseComplexChords.chords_componentCellIds.at(i);
            for (auto chord_componentCellId : chord_componentCellIds) {
                auto& componentCell = baseComplexChords.baseComplex.componentC.at(chord_componentCellId);
                for (auto componentFaceId : componentCell.Fids)
                    otherCompoentFacecolors.insert(componentFaceColors.at(componentFaceId));
            }
            std::vector<int> sortedColors(otherCompoentFacecolors.begin(), otherCompoentFacecolors.end());
            if (sortedColors == std::vector<int>{-1, 0, 1}) {
                for (auto chord_componentFaceId : chord_componentFaceIds)
                    componentFaceColors.at(chord_componentFaceId) = 2;
            } else if (sortedColors == std::vector<int>{-1, 0, 2}) {
                for (auto chord_componentFaceId : chord_componentFaceIds)
                    componentFaceColors.at(chord_componentFaceId) = 1;
            } else if (sortedColors == std::vector<int>{-1, 1, 2}) {
                for (auto chord_componentFaceId : chord_componentFaceIds)
                    componentFaceColors.at(chord_componentFaceId) = 0;
            }
            } else for (auto chord_componentFaceId : chord_componentFaceIds)
                componentFaceColors.at(chord_componentFaceId) = 3;
    }

    std::vector<size_t> colors(baseComplexChords.baseComplex.componentF.size(), MAXID);
    std::ofstream ofs("componentFaceColors.txt");

    ofs << "SCALARS " << "face_colors" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (int i = 0; i < baseComplexChords.chords_componentFaceIds.size(); ++i) {
        auto& chord_componentFaceIds = baseComplexChords.chords_componentFaceIds.at(i);
        for (const auto chord_componentFaceId : chord_componentFaceIds)
            for (const auto face_id : baseComplexChords.baseComplex.componentF.at(chord_componentFaceId).fids_patch)
                ofs << componentFaceColors.at(chord_componentFaceId) << "\n";
    }
}

void ColorFace(BaseComplex& baseComplex) {
    std::vector<int> componentFaceColors(baseComplex.componentF.size(), -1);
    for (auto& se : baseComplex.SingularityI.E) {
        std::set<int> colors = {0, 1, 2, 3, 4};
        std::vector<std::set<size_t>> separated_ComponentFidsGroups(se.neighborComponentFidsGroups.size());
        int i = 0;
        for (auto & separatedFacePatchId : se.separatedFacePatchIds) {
            for (auto face_id : baseComplex.separatedFacePatches.at(separatedFacePatchId))
                separated_ComponentFidsGroups.at(i).insert(baseComplex.mesh.F.at(face_id).componentFid);
            ++i;
        }
        for (auto& neighborComponentFidsGroup : separated_ComponentFidsGroups) {
            auto color = componentFaceColors.at(*neighborComponentFidsGroup.begin());
            if (color != -1) {
                colors.erase(color);
                continue;
            }
            for (auto neighborComponentFid : neighborComponentFidsGroup)
                componentFaceColors.at(neighborComponentFid) = *colors.begin();
            colors.erase(*colors.begin());
        }
    }

    std::ofstream ofs("componentFaceColors.txt");

    ofs << "SCALARS " << "face_colors" << " int 1\n"
        << "LOOKUP_TABLE default\n";
    for (auto& componentFace : baseComplex.componentF) {
        for (const auto face_id : componentFace.fids_patch)
            ofs << componentFaceColors.at(componentFace.id) << "\n";
    }
}
