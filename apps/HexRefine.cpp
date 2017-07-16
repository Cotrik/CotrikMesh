/*
 * HexRefine.cpp
 *
 *  Created on: Dec 5, 2016
 *      Author: cotrik
 */

const int HexRefine[8][8] =
{
    11, 20, 10, 3, 22, 26, 25, 19,
    20, 9, 2, 10, 26, 23, 18, 25,
    22, 26, 25, 19, 15, 21, 14, 7,
    26, 23, 18, 25, 21, 13, 6, 14,
    0, 8, 20, 11, 16, 24, 26, 22,
    8, 1, 9, 20, 24, 17, 23, 26,
    16, 24, 26, 22, 4, 12, 21, 15,
    24, 17, 23, 26, 12, 5, 13, 21
};
int main(int argc, char* argv[])
{
    int clockwise = 0;
    int iternum = 1;
    if (argc == 4)
    {
        iternum = atoi(argv[3]);
    }
    else if (argc == 5)
    {
        clockwise = atoi(argv[4]);
        std::cout << "clockwise oriented clockwise = " << clockwise << std::endl;
    }
    else if (argc != 3)
    {
        std::cout << "Usage: HexRefine <input_hex_file> <output_hex_file> <iternum> " << std::endl;
        return -1;
    }

    if (iternum == 1)
    {
        MeshFileReader reader(argv[1]);
        Mesh& hex_mesh = reader.GetMesh();
        hex_mesh.BuildAllConnectivities();

        Mesh new_mesh(hex_mesh);
        ////////////////////////////////////////////////////////////////////////////
        if (!new_mesh.C.empty())
        {
            const std::vector<Cell>::const_iterator iterCellBegin = new_mesh.C.begin();
            const std::vector<Cell>::const_iterator iterCellEnd = new_mesh.C.end();
            std::vector<Cell>::const_iterator iterCell = iterCellBegin;

            //if (new_mesh.m_meshType == MESH_TYPE_HEXAHEDRON || new_mesh.m_meshType == MESH_TYPE_HEXAHEDRON_VTK)
            {
                for (iterCell = iterCellBegin; iterCell != iterCellEnd; ++iterCell)
                {
                    for (int i = 0; i < 12; i++)
                    {
                        const Edge e(iterCell->at(HexEdge[i][0]), iterCell->at(HexEdge[i][1]));
                        new_mesh.E.push_back(e);
                        const DirectedEdge de1(iterCell->at(HexEdge[i][0]), iterCell->at(HexEdge[i][1]));
                        const DirectedEdge de2(iterCell->at(HexEdge[i][1]), iterCell->at(HexEdge[i][0]));
                        new_mesh.DE.push_back(de1);
                        new_mesh.DE.push_back(de2);
                    }
                }
            }
        }
        if (new_mesh.uniqueE.empty())
        for (int i = 0; i < new_mesh.E.size(); i++)
        {
            const Edge& e = new_mesh.E.at(i);
            std::vector<Edge>::iterator iterE = std::find(new_mesh.uniqueE.begin(), new_mesh.uniqueE.end(), e);
            if (iterE == new_mesh.uniqueE.end())
            {
                new_mesh.uniqueE.push_back(e);
            }
        }

        if (new_mesh.uniqueF.empty())
        {
            new_mesh.SF.resize(new_mesh.F.size());
            std::copy(new_mesh.F.begin(), new_mesh.F.end(), new_mesh.SF.begin());
            for (int i = 0; i < new_mesh.SF.size(); i++)
            {
                const SFace& f = new_mesh.SF.at(i);
                std::vector<SFace>::iterator iterF = std::find(new_mesh.uniqueF.begin(), new_mesh.uniqueF.end(), f);
                if (iterF == new_mesh.uniqueF.end())
                {
                    new_mesh.uniqueF.push_back(f);
                }
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        // add vertices
        std::vector<Vertex> new_vertex(new_mesh.V.size() + new_mesh.uniqueE.size() + new_mesh.uniqueF.size() + new_mesh.C.size());
        for (size_t i = 0; i < new_mesh.V.size(); i++)
        {
            new_vertex.at(i) = new_mesh.V.at(i);
        }
        size_t offset = new_mesh.V.size();
        for (size_t i = 0; i < new_mesh.uniqueE.size(); i++)
        {
            const Edge& e = new_mesh.uniqueE.at(i);
            const Vertex& v0 = new_mesh.V.at(e.p0);
            const Vertex& v1 = new_mesh.V.at(e.p1);
            new_vertex.at(offset + i).x = 0.5 * (v0.x + v1.x);
            new_vertex.at(offset + i).y = 0.5 * (v0.y + v1.y);
            new_vertex.at(offset + i).z = 0.5 * (v0.z + v1.z);
        }
        offset = new_mesh.V.size() + new_mesh.uniqueE.size();
        for (size_t i = 0; i < new_mesh.uniqueF.size(); i++)
        {
            const Face& e = new_mesh.uniqueF.at(i);
            const Vertex& v0 = new_mesh.V.at(e.at(0));
            const Vertex& v1 = new_mesh.V.at(e.at(2));
            new_vertex.at(offset + i).x = 0.5 * (v0.x + v1.x);
            new_vertex.at(offset + i).y = 0.5 * (v0.y + v1.y);
            new_vertex.at(offset + i).z = 0.5 * (v0.z + v1.z);
        }
        offset = new_mesh.V.size() + new_mesh.uniqueE.size() + new_mesh.uniqueF.size();
        for (size_t i = 0; i < new_mesh.C.size(); i++)
        {
            const Cell& c = new_mesh.C.at(i);
            const Vertex& v0 = new_mesh.V.at(c.at(0));
            const Vertex& v1 = new_mesh.V.at(c.at(6));
            new_vertex.at(offset + i).x = 0.5 * (v0.x + v1.x);
            new_vertex.at(offset + i).y = 0.5 * (v0.y + v1.y);
            new_vertex.at(offset + i).z = 0.5 * (v0.z + v1.z);
        }
        //new_mesh.V = new_vertex;
        /////////////////////////////////////////////////////////////////
        // add cells
        Cell cell(8, 0);
        std::vector<Cell> new_cells(8 * new_mesh.C.size(), cell);
        int count = 0;
        for (size_t i = 0; i < new_mesh.C.size(); i++)
        {
            unsigned long v_index[27];
            const Cell& c = new_mesh.C.at(i);
            for (unsigned long j = 0; j < 8; j++)
            {
                v_index[j] = c.at(j);
            }
            if (clockwise != 0)
            {
                swap(v_index[1], v_index[3]);
                swap(v_index[5], v_index[7]);
            }
            for (unsigned long j = 0; j < 12; j++)
            {
                const Edge e(c.at(HexEdge[j][0]), c.at(HexEdge[j][1]));
                std::vector<Edge>::iterator it = std::find(new_mesh.uniqueE.begin(), new_mesh.uniqueE.end(), e);
                if (it == new_mesh.uniqueE.end())
                {
                    std::cout << "Edge search Error !" << std::endl;
                }
                const unsigned long e_index = std::distance(new_mesh.uniqueE.begin(), it);
                v_index[8 + j] = new_mesh.V.size() + e_index;
            }
            for (unsigned long j = 0; j < 6; j++)
            {
                Face f(4, 0);
                f.at(0) = c.at(HexFaces[j][0]);
                f.at(1) = c.at(HexFaces[j][1]);
                f.at(2) = c.at(HexFaces[j][2]);
                f.at(3) = c.at(HexFaces[j][3]);
                const SFace sf(f);
                std::vector<SFace>::iterator it = std::find(new_mesh.uniqueF.begin(), new_mesh.uniqueF.end(), sf);
                if (it == new_mesh.uniqueF.end())
                {
                    std::cout << "Face search Error !" << std::endl;
                }
                const unsigned long f_index = std::distance(new_mesh.uniqueF.begin(), it);
                v_index[20 + j] = new_mesh.V.size() + new_mesh.uniqueE.size() + f_index;
            }
            v_index[26] = new_mesh.V.size() + new_mesh.uniqueE.size() + new_mesh.uniqueF.size() + i;
            //Cell new_cell(8, 0);
            for (int k = 0; k < 8; k++)
            {
                for (int j = 0; j < 8; j++)
                {
                    new_cells.at(count).at(j) = v_index[HexRefine[k][j]];
                }
                count++;
            }
        }
        MeshFileWriter mesh(new_vertex, new_cells, argv[2], HEXAHEDRA);
        mesh.WriteFile();
    }
}


