/*
 * FeatureLine.cpp
 *
 *  Created on: Jan 3, 2017
 *      Author: cotrik
 */

#include "FeatureLine.h"
#include <algorithm>
#include <iostream>

FeatureLine::FeatureLine(const Mesh& mesh)
: mesh(mesh)
{
    // TODO Auto-generated constructor stub

}

FeatureLine::FeatureLine(const FeatureLine& r)
: mesh(r.mesh)
, Vids(r.Vids)
, Eids(r.Eids)
{
    // TODO Auto-generated constructor stub

}

FeatureLine::~FeatureLine()
{
    // TODO Auto-generated destructor stub
}

void FeatureLine::Extract(const size_t label)
{
    for (size_t i = 0; i < mesh.E.size(); i++) {
        const Edge& edge = mesh.E.at(i);
        if (edge.isSharpFeature && edge.label == label) {
            Eids.push_back(edge.id);
            Vids.push_back(edge.Vids[0]);
            Vids.push_back(edge.Vids[1]);
        }
    }
    {
        std::sort(Vids.begin(), Vids.end());
        std::vector<size_t>::iterator iter = std::unique(Vids.begin(), Vids.end());
        Vids.resize(std::distance(Vids.begin(), iter));
    }
    {
        std::sort(Eids.begin(), Eids.end());
        std::vector<size_t>::iterator iter = std::unique(Eids.begin(), Eids.end());
        Eids.resize(std::distance(Eids.begin(), iter));
    }
    Vertex* pCorner = NULL;
    for (size_t i = 0; i < Vids.size(); i++) {
        const Vertex& v = mesh.V.at(Vids.at(i));
        if (v.isCorner) {
            pCorner = (Vertex*)&v;
            break;
        }
    }
    std::vector<size_t> newVids;
    std::vector<size_t> newEids;
    std::vector<bool> visited(Vids.size(), false);
    //if (pCorner != NULL)
    {
        // begin with corner
        Vertex* pVertex = pCorner != NULL ? pCorner : (Vertex*)&mesh.V.at(Vids.at(0));
        size_t vid = pVertex->id;
        //if (pCorner != NULL)
        {
            for (size_t i = 0; i < Vids.size(); i++) {
                if (Vids[i] == vid) {
                    visited[i] = true;
                    break;
                }
            }
            newVids.push_back(vid);
        }
        bool allVisited = false;
        while (!allVisited) {
            size_t eid = MAXID;
            size_t oldvid = vid;
            vid = FindNextVid(vid, visited, eid);
            if (vid != MAXID) {
                newVids.push_back(vid);
            } 
            /*else {
                std::ofstream ofs1("Vids.vtk");
                ofs1 << "# vtk DataFile Version 3.0\n"
                	<< "output.vtk\n"
                	<< "ASCII\n\n"
                	<< "DATASET UNSTRUCTURED_GRID\n";
                ofs1 << "POINTS " << mesh.V.size() << " double\n";
                
                for (size_t i = 0; i < mesh.V.size(); i++) {
                	ofs1 << mesh.V.at(i).x << " " <<  mesh.V.at(i).y << " " <<  mesh.V.at(i).z << "\n";
                }
                ofs1 << "CELLS " << Vids.size() << " " << 2 * Vids.size() << std::endl;
                for (size_t i = 0; i < Vids.size(); i++) {
                	ofs1 << "1 " << Vids.at(i) << std::endl;
                }
                ofs1 << "CELL_TYPES " << Vids.size() << "\n";
                for (size_t i = 0; i < Vids.size(); i++) {
                	ofs1 << "1" << std::endl;
                }

                std::ofstream ofs2("Eids.vtk");
                ofs2 << "# vtk DataFile Version 3.0\n"
                	<< "output.vtk\n"
                	<< "ASCII\n\n"
                	<< "DATASET UNSTRUCTURED_GRID\n";
                ofs2 << "POINTS " << mesh.V.size() << " double\n";
                
                for (size_t i = 0; i < mesh.V.size(); i++) {
                	ofs2 << mesh.V.at(i).x << " " <<  mesh.V.at(i).y << " " <<  mesh.V.at(i).z << "\n";
                }
                ofs2 << "CELLS " << Eids.size() << " " << 3 * Eids.size() << std::endl;
                for (size_t i = 0; i < Eids.size(); i++) {
                	ofs2 << "2 " << mesh.E.at(Eids.at(i)).Vids[0] << mesh.E.at(Eids.at(i)).Vids[1] << std::endl;
                }
                ofs2 << "CELL_TYPES " << Eids.size() << "\n";
                for (size_t i = 0; i < Eids.size(); i++) {
                	ofs2 << "3" << std::endl;
                }
            }*/
            newEids.push_back(eid);
            allVisited = true;
            for (size_t i = 0; i < Vids.size(); i++) {
                if (!visited[i]) {
                    allVisited = false;
                    break;
                }
            }
        }
        if (pCorner == NULL)
            newVids.push_back(newVids[0]);
    }
    Vids = newVids;
    Eids = newEids;
}

size_t FeatureLine::FindNextVid(const size_t vid, std::vector<bool>& visited, size_t& nextEid)
{
    const Vertex& v = mesh.V.at(vid);
    for (size_t i = 0; i < v.N_Vids.size(); i++)
        for (size_t j = 0; j < Vids.size(); j++)
            if (Vids[j] == v.N_Vids[i] && !visited[j]){
                nextEid = v.N_Eids[i];
                visited[j] = true;
                return v.N_Vids[i];
            }

    std::cout << "Error in FindNextVid! vid = " << vid << std::endl;
    return MAXID;
}
