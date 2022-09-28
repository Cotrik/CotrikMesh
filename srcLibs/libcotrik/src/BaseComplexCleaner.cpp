/*
 * BaseComplexCleaner.cpp
 *
 *  Created on: Aug 28, 2018
 *      Author: cotrik
 */

#include "BaseComplexCleaner.h"
#include "MeshFileWriter.h"
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <queue>
#include <iterator>
#include <Eigen/Dense>

const std::string red_color = "\033[1;31m";
const std::string end_color = "\033[0m";

BaseComplexCleaner::BaseComplexCleaner(BaseComplex& baseComplex)
        : baseComplex(baseComplex) {
    // TODO Auto-generated constructor stub

}

BaseComplexCleaner::~BaseComplexCleaner() {
    // TODO Auto-generated destructor stub
}

void BaseComplexCleaner::Run() {
    baseComplex.BuildSingularityConnectivity();

    std::vector<int> celldata(baseComplex.mesh.C.size(), 0);
    int cubelabel = 1;
    std::vector<bool> labeled(baseComplex.SingularityI.E.size(), false);
    for (auto& sv : baseComplex.SingularityI.V) {
        if (baseComplex.mesh.V.at(sv.id_mesh).isBoundary) continue;
        if (/*sv.N_Vids.size() != 4 || */NeighborOnBoundary(sv)) continue;
        if (!AllNeighborsValence4(sv)) continue;
        if (CubeStructureLabeled(sv, labeled)) continue;
        if (!IsCubeStructureCorrect(sv)) continue;
        LabelCubeStructure(sv, labeled, cubelabel, celldata);
    }

    MeshFileWriter writer(baseComplex.mesh, "out.vtk");
    writer.WriteFile();
    writer.WriteCellData(celldata, "cubeid");
    std::cout << "Finished\n";
}

bool BaseComplexCleaner::NeighborOnBoundary(const SingularV& sv) const {
    for (auto& nsvid : sv.N_Vids) {
        auto& nsv = baseComplex.SingularityI.V.at(nsvid);
        if (baseComplex.mesh.V.at(nsv.id_mesh).isBoundary) return true;
    }
    return false;
}

bool BaseComplexCleaner::AllNeighborsValence4(const SingularV& sv) const {
    for (auto& nsvid : sv.N_Vids) {
        auto& nsv = baseComplex.SingularityI.V.at(nsvid);
        if (nsv.N_Vids.size() != 4) return false;
    }
    return true;
}

bool BaseComplexCleaner::CubeStructureLabeled(const SingularV& sv, const std::vector<bool>& labeled) const {
    bool has_been_labeled = false;
    for (auto& seid : sv.N_Eids) {
        if (labeled.at(seid)) return true;

        bool neighbor_has_been_labeled = false;
        for (auto& nsvid : sv.N_Vids) {
            auto& nsv = baseComplex.SingularityI.V.at(nsvid);
            for (auto& nseid : nsv.N_Eids)
                if (labeled.at(nseid)) return true;
        }
    }
    return false;
}

//bool BaseComplexCleaner::IsCubeStructureCorrect(const SingularV& sv) const {
//    std::set<size_t> svids;
//    for (auto& nsvid : sv.N_Vids) {
//        auto& nsv = baseComplex.SingularityI.V.at(nsvid);
//        if (nsv.N_Vids.size() != 4) return false;
//        svids.insert(nsv.N_Vids.begin(), nsv.N_Vids.end());
//    }
//    return svids.size() == (sv.N_Vids.size() * 5 - sv.N_Vids.size() + 1);
//}

bool BaseComplexCleaner::IsCubeStructureCorrect(const SingularV& sv) const {
    std::set<size_t> svids;
    int count = 0;
    for (auto& nsvid : sv.N_Vids) {
        auto& nsv = baseComplex.SingularityI.V.at(nsvid);
        if (nsv.N_Vids.size() != 4) return false;
        bool intersected = false;
        for (auto& nsv_nvid : nsv.N_Vids)
            if (nsv_nvid != sv.id && svids.find(nsv_nvid) != svids.end()) {
                intersected = true;
                break;
            }
        if (!intersected) ++count;
        svids.insert(nsv.N_Vids.begin(), nsv.N_Vids.end());
    }

    return count >= 4;
    //return svids.size() == (sv.N_Vids.size() * 5 - sv.N_Vids.size() + 1);
}

void BaseComplexCleaner::LabelCubeStructure(const SingularV& sv, std::vector<bool>& labeled, int& cubelabel, std::vector<int>& celldata) const {
    for (auto& seid : sv.N_Eids) {
        for (auto& nsvid : sv.N_Vids) {
            auto& nsv = baseComplex.SingularityI.V.at(nsvid);
            for (auto& nseid : nsv.N_Eids) {
                // labeled.at(nseid) = true;
                auto& nse = baseComplex.SingularityI.E.at(nseid);
                for (auto eid : nse.es_link) {
                    auto& e = baseComplex.mesh.E.at(eid);
                    for (auto cid : e.N_Cids)
                        celldata[cid] = cubelabel;
                }
            }
        }
    }
    ++cubelabel;
}
