/*
 * tet2hex.cpp
 *
 *  Created on: Sep 11, 2017
 *      Author: cotrik
 */


#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "ArgumentManager.h"
#include <algorithm>

size_t GetEdgeId(Mesh& tetMesh, const size_t vid1, const size_t vid2) {
    for (const auto& e : tetMesh.E) {
        const auto& v1 = tetMesh.V.at(e.Vids[0]);
        const auto& v2 = tetMesh.V.at(e.Vids[1]);
        if ((v1.id == vid1 && v2.id == vid2) || (v1.id == vid2 && v2.id == vid1)) return e.id;
    }
    std::cerr << "error in GetEdgeId\n";
    return MAXID;
}
size_t GetFaceId(Mesh& tetMesh, std::vector<size_t> vids) {
    std::sort(vids.begin(), vids.end());
    for (const auto& f : tetMesh.F) {
        std::vector<size_t> vs = f.Vids;
        std::sort(vs.begin(), vs.end());
       //if (vids == vs) return f.id;
        if (vids[0] == vs[0] && vids[1] == vs[1] && vids[2] == vs[2]) return f.id;
    }
    std::cerr << "error in GetFaceId\n";
    return MAXID;
}
int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: tet2hex tet_file hex_file\n";
        return -1;
    }
    MeshFileReader reader(argv[1]);
    Mesh& tetMesh = (Mesh&)reader.GetMesh();
    tetMesh.BuildAllConnectivities();
    tetMesh.ExtractBoundary();

    Mesh hexMesh;
    hexMesh.V.resize(tetMesh.V.size() + tetMesh.E.size() + tetMesh.F.size() + tetMesh.C.size());
    hexMesh.C.resize(4 * tetMesh.C.size());
    // Construct V;
    size_t vid = 0;
    for (const auto& v : tetMesh.V) {
        hexMesh.V.at(vid) = v.xyz();
        hexMesh.V.at(vid).isBoundary = v.isBoundary;
        hexMesh.V.at(vid).id = vid++;
    }
    for (const auto& e : tetMesh.E) {
        const auto& v1 = tetMesh.V.at(e.Vids[0]);
        const auto& v2 = tetMesh.V.at(e.Vids[1]);
        hexMesh.V.at(vid) = 0.5 * (v1 + v2);
        hexMesh.V.at(vid).isBoundary = e.isBoundary;
        hexMesh.V.at(vid).id = vid++;
    }
    for (const auto& f : tetMesh.F) {
        const auto& v1 = tetMesh.V.at(f.Vids[0]);
        const auto& v2 = tetMesh.V.at(f.Vids[1]);
        const auto& v3 = tetMesh.V.at(f.Vids[2]);
        hexMesh.V.at(vid) = 0.333333 * (v1 + v2 + v3);
        hexMesh.V.at(vid).isBoundary = f.isBoundary;
        hexMesh.V.at(vid).id = vid++;
    }
    for (const auto& c : tetMesh.C) {
        const auto& v1 = tetMesh.V.at(c.Vids[0]);
        const auto& v2 = tetMesh.V.at(c.Vids[1]);
        const auto& v3 = tetMesh.V.at(c.Vids[2]);
        const auto& v4 = tetMesh.V.at(c.Vids[3]);
        hexMesh.V.at(vid) = 0.25 * (v1 + v2 + v3 + v4);
        hexMesh.V.at(vid).isBoundary = false;
        hexMesh.V.at(vid).id = vid++;
    }
    // Construct C;
    size_t cid = 0;
    for (const auto& c : tetMesh.C) {
        // V
        const auto& v0 = hexMesh.V.at(c.Vids[0]);
        const auto& v1 = hexMesh.V.at(c.Vids[1]);
        const auto& v2 = hexMesh.V.at(c.Vids[2]);
        const auto& v3 = hexMesh.V.at(c.Vids[3]);
        // E
        const auto& v4 = hexMesh.V.at(GetEdgeId(tetMesh, c.Vids[0], c.Vids[1]) + tetMesh.V.size());
        const auto& v5 = hexMesh.V.at(GetEdgeId(tetMesh, c.Vids[1], c.Vids[2]) + tetMesh.V.size());
        const auto& v6 = hexMesh.V.at(GetEdgeId(tetMesh, c.Vids[0], c.Vids[2]) + tetMesh.V.size());
        const auto& v7 = hexMesh.V.at(GetEdgeId(tetMesh, c.Vids[2], c.Vids[3]) + tetMesh.V.size());
        const auto& v8 = hexMesh.V.at(GetEdgeId(tetMesh, c.Vids[0], c.Vids[3]) + tetMesh.V.size());
        const auto& v9 = hexMesh.V.at(GetEdgeId(tetMesh, c.Vids[1], c.Vids[3]) + tetMesh.V.size());
        // F
        const auto& v10 = hexMesh.V.at(GetFaceId(tetMesh, std::vector<size_t>{c.Vids[0], c.Vids[1], c.Vids[2]}) + tetMesh.V.size() + tetMesh.E.size());
        const auto& v11 = hexMesh.V.at(GetFaceId(tetMesh, std::vector<size_t>{c.Vids[0], c.Vids[2], c.Vids[3]}) + tetMesh.V.size() + tetMesh.E.size());
        const auto& v12 = hexMesh.V.at(GetFaceId(tetMesh, std::vector<size_t>{c.Vids[0], c.Vids[1], c.Vids[3]}) + tetMesh.V.size() + tetMesh.E.size());
        const auto& v13 = hexMesh.V.at(GetFaceId(tetMesh, std::vector<size_t>{c.Vids[1], c.Vids[2], c.Vids[3]}) + tetMesh.V.size() + tetMesh.E.size());
        // C
        const auto& v14 = hexMesh.V.at(c.id + tetMesh.V.size() + tetMesh.E.size() + tetMesh.F.size());

        std::vector<size_t> c1 = {v0.id, v8.id, v12.id, v4.id, v6.id, v11.id, v14.id, v10.id};
        std::vector<size_t> c2 = {v1.id, v5.id, v10.id, v4.id, v9.id, v13.id, v14.id, v12.id};
        std::vector<size_t> c3 = {v2.id, v5.id, v13.id, v7.id, v6.id, v10.id, v14.id, v11.id};
        std::vector<size_t> c4 = {v3.id, v8.id, v11.id, v7.id, v9.id, v12.id, v14.id, v13.id};

        hexMesh.C.at(cid).Vids = c1;
        hexMesh.C.at(cid).id = cid++;
        hexMesh.C.at(cid).Vids = c2;
        hexMesh.C.at(cid).id = cid++;
        hexMesh.C.at(cid).Vids = c3;
        hexMesh.C.at(cid).id = cid++;
        hexMesh.C.at(cid).Vids = c4;
        hexMesh.C.at(cid).id = cid++;
    }
    MeshFileWriter writer(hexMesh, argv[2]);
    writer.WriteFile();
    return 0;
}
