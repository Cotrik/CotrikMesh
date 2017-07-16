#include "io.h"

h_io::h_io(void)
{
    counter = -1;
}

void h_io::read_triangle_mesh_off(vector<Vertex> &Vs, vector<Triangle> &Ts, const char * fname)
{
    char file[300];
    std::fstream f(fname, std::ios::in);
    char s[1024], sread[1024];
    int vnum, tnum;
    float x, y, z;
    f.getline(s, 1023);
    sscanf(s, "%s", &sread);
    f.getline(s, 1023);
    sscanf(s, "%d%d%f", &vnum, &tnum, &x);
    for (int i = 0; i < vnum; i++)
    {
        f.getline(s, 1023);
        sscanf(s, "%f %f %f", &x, &y, &z);
        Vertex v;
        v.v[0] = x;
        v.v[1] = y;
        v.v[2] = z;
        v.index = i;
        Vs.push_back(v);
    }
    for (int i = 0; i < tnum; i++)
    {
        f.getline(s, 1023);
        int num, a, b, c;
        sscanf(s, "%d %d %d %d", &num, &a, &b, &c);
        Triangle t;

        t.triangle_v[0] = a;
        t.triangle_v[1] = b;
        t.triangle_v[2] = c;
        t.index = i;

        Ts.push_back(t);

        Vs[t.triangle_v[0]].neighbort.push_back(t.index);
        Vs[t.triangle_v[1]].neighbort.push_back(t.index);
        Vs[t.triangle_v[2]].neighbort.push_back(t.index);

        Vs[t.triangle_v[0]].neighborv.push_back(t.triangle_v[1]);
        Vs[t.triangle_v[0]].neighborv.push_back(t.triangle_v[2]);
        Vs[t.triangle_v[1]].neighborv.push_back(t.triangle_v[0]);
        Vs[t.triangle_v[1]].neighborv.push_back(t.triangle_v[2]);
        Vs[t.triangle_v[2]].neighborv.push_back(t.triangle_v[0]);
        Vs[t.triangle_v[2]].neighborv.push_back(t.triangle_v[1]);
    }
    for (int i = 0; i < Vs.size(); i++)
    {
        std::vector<int> neighborvs = Vs[i].neighborv;
        Vs[i].neighborv.clear();
        for (int j = 0; j < neighborvs.size(); j++)
        {
            bool havesame = false;
            for (int k = j + 1; k < neighborvs.size(); k++)
                if (neighborvs[k] == neighborvs[j])
                    havesame = true;
            if (!havesame)
                Vs[i].neighborv.push_back(neighborvs[j]);
        }
    }

    f.close();
}
void h_io::read_hex_mesh_off(vector<Hex_V> &Vs, vector<Hex> &Hexs, const char * fname)
{
    Vs.clear();
    Hexs.clear();

    char file[300];
// 	sprintf(file,"%s",HEX_MESH_PATH_READ_OFF);
    std::fstream f(fname, std::ios::in);
    char s[1024], sread[1024];
    int vnum, tnum;
    float x, y, z;
    f.getline(s, 1023);
    sscanf(s, "%s", &sread);
    f.getline(s, 1023);
    sscanf(s, "%d%d%f", &vnum, &tnum, &x);
    for (int i = 0; i < vnum; i++)
    {
        f.getline(s, 1023);
        int temp = -1;
        //sscanf(s,"%f %f %f",&x,&y,&z);
        sscanf(s, "%f %f %f %d", &x, &y, &z, &temp);
        Hex_V v;
        v.v[0] = x;
        v.v[1] = y;
        v.v[2] = z;
        v.index = i;
        v.where_location = -1;
        v.fixed = false;
        v.slice_id = temp;
        Vs.push_back(v);
    }
    int hid = 0;
    for (int i = 0; i < tnum; i++)
    {
        f.getline(s, 1023);
        int a, b, c, d, e, f, g, m;
        sscanf(s, "%d %d %d %d %d %d %d %d %d", &vnum, &a, &b, &c, &d, &e, &f, &g, &m);
        Hex h;

        //a--;b--;c--;d--;e--;f--;g--;m--;

        h.V_Ids[0] = a;
        h.V_Ids[1] = b;
        h.V_Ids[2] = c;
        h.V_Ids[3] = d;
        h.V_Ids[4] = e;
        h.V_Ids[5] = f;
        h.V_Ids[6] = g;
        h.V_Ids[7] = m;
// 		h.v_ids.push_back(vnum-1);h.v_ids.push_back((int)x-1);h.v_ids.push_back((int)y-1);h.v_ids.push_back((int)z-1);
// 		h.v_ids.push_back(a-1);h.v_ids.push_back(b-1);h.v_ids.push_back(c-1);h.v_ids.push_back(d-1);

//  		bool all=true;//deal with deformed torus base-complex
//  		for(int j = 0; j < 8; j++)
//  			if(!(h.V_Ids[j] >= 0 && h.V_Ids[j] <= 961))
//  				all = false;
//  		if(all)
//  			continue;

        h.index = hid++;
        h.frame_component = -1;
        Hexs.push_back(h);

        Vs[a].neighbor_Hs.push_back(h.index);
        Vs[b].neighbor_Hs.push_back(h.index);
        Vs[c].neighbor_Hs.push_back(h.index);
        Vs[d].neighbor_Hs.push_back(h.index);
        Vs[e].neighbor_Hs.push_back(h.index);
        Vs[f].neighbor_Hs.push_back(h.index);
        Vs[g].neighbor_Hs.push_back(h.index);
        Vs[m].neighbor_Hs.push_back(h.index);
    }

    f.close();
}

void h_io::write_hex_mesh_off(vector<Hex_V> &Vs, vector<Hex> &Hexs, const char * fname)
{
    char file[300];
    fstream f(fname, ios::out);

    f << "OFF" << endl;
    f << Vs.size() << " " << Hexs.size() << " " << 0 << endl;
    for (int i = 0; i < Vs.size(); i++)
        //f<<setprecision(10)<<Vs[i].v[0]<<" "<<setprecision(10)<<Vs[i].v[1]<<" "<<setprecision(10)<<Vs[i].v[2]<<endl;
        f << setprecision(10) << Vs[i].v[0] << " " << setprecision(10) << Vs[i].v[1] << " " << setprecision(10) << Vs[i].v[2] << " " << Vs[i].slice_id << endl;

    for (int i = 0; i < Hexs.size(); i++)
    {
        f << "10 " << Hexs[i].V_Ids[0] << " " << Hexs[i].V_Ids[1] << " " << Hexs[i].V_Ids[2] << " " << Hexs[i].V_Ids[3];
        f << " " << Hexs[i].V_Ids[4] << " " << Hexs[i].V_Ids[5] << " " << Hexs[i].V_Ids[6] << " " << Hexs[i].V_Ids[7] << " 0 0" << endl;
    }

    f.close();
}

void h_io::write_VTK(const vector<Hex_V> &Vs, const vector<Hex> &Hexs, const char * fname)
{
    char file_meshlab_[300];
    std::fstream f_out_meshlab_(fname, std::ios::out);

    f_out_meshlab_ << "# vtk DataFile Version 2.0" << std::endl << "hex mesh vtk data - converted from .off" << std::endl;
    f_out_meshlab_ << "ASCII" << std::endl;
    f_out_meshlab_ << "DATASET UNSTRUCTURED_GRID" << std::endl;

    f_out_meshlab_ << "POINTS " << Vs.size() << " double" << std::endl;

    for (int i = 0; i < Vs.size(); i++)
        f_out_meshlab_ << Vs[i].v[0] << " " << Vs[i].v[1] << " " << Vs[i].v[2] << std::endl;

    f_out_meshlab_ << "CELLS " << Hexs.size() << " " << Hexs.size() * 9 << std::endl;
    for (int i = 0; i < Hexs.size(); i++)
    {
        f_out_meshlab_ << " " << 8 << " " << Hexs[i].V_Ids[0] << " " << Hexs[i].V_Ids[1] << " " << Hexs[i].V_Ids[2] << " " << Hexs[i].V_Ids[3];
        f_out_meshlab_ << " " << Hexs[i].V_Ids[4] << " " << Hexs[i].V_Ids[5] << " " << Hexs[i].V_Ids[6] << " " << Hexs[i].V_Ids[7] << std::endl;
    }

    f_out_meshlab_ << "CELL_TYPES " << Hexs.size() << std::endl;
    for (int i = 0; i < Hexs.size(); i++)
        f_out_meshlab_ << 12 << std::endl;

    f_out_meshlab_ << "POINT_DATA " << Vs.size() << std::endl;
    f_out_meshlab_ << "SCALARS fixed int" << std::endl;
    f_out_meshlab_ << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < Vs.size(); i++)
    {
        if (Vs[i].where_location == 1)
        {
            f_out_meshlab_ << "1" << std::endl;
        }
        else
        {
            f_out_meshlab_ << "0" << std::endl;
        }
    }
    f_out_meshlab_ << "CELL_DATA " << Hexs.size() << std::endl;
    f_out_meshlab_ << "SCALARS Hfixed int" << std::endl;
    f_out_meshlab_ << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < Hexs.size(); i++)
    {
        f_out_meshlab_ << Hexs[i].Color_ID << std::endl;
    }

    f_out_meshlab_ << "SCALARS Efixed int" << std::endl;
    f_out_meshlab_ << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < Hexs.size(); i++)
    {
        f_out_meshlab_ << Hexs[i].Color_E_ID << std::endl;
    }

    f_out_meshlab_ << "SCALARS EFfixed int" << std::endl;
    f_out_meshlab_ << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < Hexs.size(); i++)
    {
        f_out_meshlab_ << Hexs[i].Color_EF_ID << std::endl;
    }
    f_out_meshlab_.close();
}
void h_io::read_VTK(vector<Hex_V> &Vs, vector<Hex> &Hexs, const char * fname)
{
    Vs.clear();
    Hexs.clear();

    char file[300];
    std::fstream ff(fname, std::ios::in);
    char s[1024], sread[1024], sread2[1024];
    int vnum, hnum;
    float x, y, z;
    int is_boundary;

    int find = false;
    while (!find)
    {
        ff.getline(s, 1023);
        if (sscanf(s, "%s %d %s", &sread, &vnum, &sread2) == 3 && (strcmp(sread, "POINTS") == 0))
            find = true;
    }

    for (int i = 0; i < vnum; i++)
    {
        ff.getline(s, 1023);
        int framid;
        sscanf(s, "%f %f %f", &x, &y, &z);

        Hex_V v;
        //Vs[i].v[0]=x;Vs[i].v[1]=y;Vs[i].v[2]=z;
        v.v[0] = x;
        v.v[1] = y;
        v.v[2] = z;
        Vs.push_back(v);
    }

    find = false;
    while (!find)
    {
        int temp_int;
        ff.getline(s, 1023);
        if (sscanf(s, "%s %d %d", &sread, &hnum, &temp_int) == 3 && (strcmp(sread, "CELLS") == 0))
            find = true;
    }

    for (int i = 0; i < hnum; i++)
    {
        int a, b, c, d, e, f, g, m, nn, pp;
        ff.getline(s, 1023);
        sscanf(s, "%d %d %d %d %d %d %d %d %d", &vnum, &a, &b, &c, &d, &e, &f, &g, &m);
        Hex h;
//      Hexs[i].V_Ids[0]=(int)a;Hexs[i].V_Ids[1]=(int)b;Hexs[i].V_Ids[2]=(int)c;
//      Hexs[i].V_Ids[3]=(int)d;Hexs[i].V_Ids[4]=(int)e;Hexs[i].V_Ids[5]=(int)f;
//      Hexs[i].V_Ids[6]=(int)g;Hexs[i].V_Ids[7]=(int)m;

        h.V_Ids[0] = (int) a;
        h.V_Ids[1] = (int) b;
        h.V_Ids[2] = (int) c;
        h.V_Ids[3] = (int) d;
        h.V_Ids[4] = (int) e;
        h.V_Ids[5] = (int) f;
        h.V_Ids[6] = (int) g;
        h.V_Ids[7] = (int) m;
        Hexs.push_back(h);
    }
    ff.close();
}

void h_io::save_base_complex_VTK(const char *fname)
{
    char file_meshlab_[300];
    std::fstream ff(fname, std::ios::out);

    ff << "# vtk DataFile Version 2.0" << std::endl << "hex mesh vtk data - converted from .off" << std::endl;
    ff << "ASCII" << std::endl;
    ff << "DATASET POLYDATA" << std::endl;

    ff << "POINTS " << hex_mesh.HVs.size() << " double" << std::endl;

    for (int i = 0; i < hex_mesh.HVs.size(); i++)
        ff << hex_mesh.HVs[i].v[0] << " " << hex_mesh.HVs[i].v[1] << " " << hex_mesh.HVs[i].v[2] << std::endl;

    ff << "VERTICES " << FrameI.FVs.size() << " " << FrameI.FVs.size() * 2 << endl;
    for (int i = 0; i < FrameI.FVs.size(); i++)
    {
        ff << "1 " << FrameI.FVs[i].index_hex << endl;
    }

    int line_num = 0;
    for (int i = 0; i < FrameI.FEs.size(); i++)
    {
        line_num += FrameI.FEs[i].vs_link.size();
    }

    ff << "LINES " << FrameI.FEs.size() << " " << line_num + FrameI.FEs.size() << endl;
    for (int i = 0; i < FrameI.FEs.size(); i++)
    {
        ff << FrameI.FEs[i].vs_link.size() << " ";
        for (int j = 0; j < FrameI.FEs[i].vs_link.size(); j++)
            ff << FrameI.FEs[i].vs_link[j] << " ";
        ff << endl;
    }

    int polygon_num = 0;
    for (int i = 0; i < FrameI.FFs.size(); i++)
    {
        polygon_num += FrameI.FFs[i].hfs_net_another.size();
    }

    ff << "POLYGONS " << polygon_num << " " << polygon_num * 5 << endl;
    for (int i = 0; i < FrameI.FFs.size(); i++)
    {
        for (int j = 0; j < FrameI.FFs[i].hfs_net_another.size(); j++)
        {
            ff << "4 " << hex_mesh.HFs[FrameI.FFs[i].hfs_net_another[j]].cv_Ids[0] << " ";
            ff << hex_mesh.HFs[FrameI.FFs[i].hfs_net_another[j]].cv_Ids[1] << " " << hex_mesh.HFs[FrameI.FFs[i].hfs_net_another[j]].cv_Ids[2] << " ";
            ff << hex_mesh.HFs[FrameI.FFs[i].hfs_net_another[j]].cv_Ids[3] << endl;
        }
    }

    ff << "POINT_DATA " << hex_mesh.HVs.size() << std::endl;
    ff << "SCALARS V_Scalars int" << std::endl;
    ff << "LOOKUP_TABLE V_Table" << std::endl;
    for (int i = 0; i < hex_mesh.HVs.size(); i++)
    {
        if (hex_mesh.HVs[i].where_location == 1)
        {
            bool found = false;
            for (int j = 0; j < FrameI.FVs.size(); j++)
            {
                if (FrameI.FVs[j].index_hex == i)
                    found = true;
            }

            bool found2 = false;
            for (int j = 0; j < FrameI.FEs.size(); j++)
            {
                for (int k = 0; k < FrameI.FEs[j].vs_link.size(); k++)
                    if (FrameI.FEs[j].vs_link[k] == i)
                        found2 = true;
            }

            if (found)
                ff << "3" << std::endl; //Frame_V
            else if (found2)
                ff << "2" << std::endl; //Frame_Edge
            else
                ff << "1" << std::endl;
        }
        else
        {
            bool found = false;
            for (int j = 0; j < FrameI.FVs.size(); j++)
            {
                if (FrameI.FVs[j].index_hex == i)
                    found = true;
            }

            bool found2 = false;
            for (int j = 0; j < FrameI.FEs.size(); j++)
            {
                for (int k = 0; k < FrameI.FEs[j].vs_link.size(); k++)
                    if (FrameI.FEs[j].vs_link[k] == i)
                        found2 = true;
            }

            if (found)
                ff << "3" << std::endl; //Frame_V
            else if (found2)
                ff << "2" << std::endl; //Frame_Edge
            else
                ff << "0" << std::endl;
        }
    }

    /*ff<< "POINT_DATA " << hex_mesh.HVs.size() << std::endl;
     ff<< "KNOTS fixed int" <<" "<<1<< std::endl;
     ff << "LOOKUP_TABLE default" << std::endl;
     for(int i=0;i<hex_mesh.HVs.size();i++)
     {
     bool found=false;
     for(int j=0;j<FrameI.FVs.size();j++)
     {
     if(i==FrameI.FVs[j].index_hex)
     found=true;
     }
     if (found)
     {
     ff<< "1" << std::endl;
     }
     else
     {
     ff<< "0" << std::endl;
     }
     }*/
    ff.close();
}

void h_io::save_BasecomplexNode_VTK(const char *fname)
{
    char file_meshlab_[300];
    std::fstream ff(fname, std::ios::out);

    ff << "# vtk DataFile Version 2.0" << std::endl << "hex mesh vtk data - converted from .off" << std::endl;
    ff << "ASCII" << std::endl;
    ff << "DATASET POLYDATA" << std::endl;
    ff << "POINTS " << FrameI.FVs.size() << " double" << std::endl;

    for (int i = 0; i < FrameI.FVs.size(); i++)
    {
        int hvid = FrameI.FVs[i].index_hex;
        ff << hex_mesh.HVs[hvid].v[0] << " " << hex_mesh.HVs[hvid].v[1] << " " << hex_mesh.HVs[hvid].v[2] << std::endl;
    }

    ff << "POINT_DATA " << FrameI.FVs.size() << std::endl;
    ff << "SCALARS V_Scalars int" << std::endl;
    ff << "LOOKUP_TABLE V_Table" << std::endl;
    for (int i = 0; i < FrameI.FVs.size(); i++)
    {
        ff << "1" << std::endl;
    }

    ff.close();
}

void h_io::save_BasecomplexEdge_VTK(const char *fname)
{
    char file_meshlab_[300];
    std::fstream ff(fname, std::ios::out);

    ff << "# vtk DataFile Version 2.0" << std::endl << "hex mesh vtk data - converted from .off" << std::endl;
    ff << "ASCII" << std::endl;
    ff << "DATASET POLYDATA" << std::endl;

    ff << "POINTS " << hex_mesh.HVs.size() << " double" << std::endl;

    for (int i = 0; i < hex_mesh.HVs.size(); i++)
        ff << hex_mesh.HVs[i].v[0] << " " << hex_mesh.HVs[i].v[1] << " " << hex_mesh.HVs[i].v[2] << std::endl;

    int line_num = 0;
    for (int i = 0; i < FrameI.FEs.size(); i++)
    {
        line_num += FrameI.FEs[i].vs_link.size();
    }

    ff << "LINES " << FrameI.FEs.size() << " " << line_num + FrameI.FEs.size() << endl;
    for (int i = 0; i < FrameI.FEs.size(); i++)
    {
        ff << FrameI.FEs[i].vs_link.size() << " ";
        for (int j = 0; j < FrameI.FEs[i].vs_link.size(); j++)
            ff << FrameI.FEs[i].vs_link[j] << " ";
        ff << endl;
    }

    ff << "POINT_DATA " << hex_mesh.HVs.size() << std::endl;
    ff << "SCALARS V_Scalars int" << std::endl;
    ff << "LOOKUP_TABLE V_Table" << std::endl;
    for (int i = 0; i < hex_mesh.HVs.size(); i++)
    {
        if (hex_mesh.HVs[i].where_location == 1)
        {
            bool found = false;
            for (int j = 0; j < SingularityI.SVs.size(); j++)
            {
                if (SingularityI.SVs[j].index_hex == i)
                    found = true;
            }

            bool found2 = false;
            for (int j = 0; j < SingularityI.SEs.size(); j++)
            {
                for (int k = 0; k < SingularityI.SEs[j].vs_link.size(); k++)
                    if (SingularityI.SEs[j].vs_link[k] == i)
                        found2 = true;
            }

            if (found)
                ff << "3" << std::endl; //SV
            else if (found2)
                ff << "2" << std::endl; //SE
            else
                ff << "1" << std::endl;
        }
        else
        {
            bool found = false;
            for (int j = 0; j < SingularityI.SVs.size(); j++)
            {
                if (SingularityI.SVs[j].index_hex == i)
                    found = true;
            }

            bool found2 = false;
            for (int j = 0; j < SingularityI.SEs.size(); j++)
            {
                for (int k = 0; k < SingularityI.SEs[j].vs_link.size(); k++)
                    if (SingularityI.SEs[j].vs_link[k] == i)
                        found2 = true;
            }

            if (found)
                ff << "3" << std::endl; //SV
            else if (found2)
                ff << "2" << std::endl; //SE
            else
                ff << "0" << std::endl;
        }
    }

    ff.close();
}

void h_io::save_singularG_VTK(const char *fname)
{
    char file_meshlab_[300];
    std::fstream ff(fname, std::ios::out);

    ff << "# vtk DataFile Version 2.0" << std::endl << "hex mesh vtk data - converted from .off" << std::endl;
    ff << "ASCII" << std::endl;
    ff << "DATASET POLYDATA" << std::endl;

    ff << "POINTS " << hex_mesh.HVs.size() << " double" << std::endl;

    for (int i = 0; i < hex_mesh.HVs.size(); i++)
        ff << hex_mesh.HVs[i].v[0] << " " << hex_mesh.HVs[i].v[1] << " " << hex_mesh.HVs[i].v[2] << std::endl;

    ff << "VERTICES " << SingularityI.SVs.size() << " " << SingularityI.SVs.size() * 2 << endl;
    for (int i = 0; i < SingularityI.SVs.size(); i++)
    {
        ff << "1 " << SingularityI.SVs[i].index_hex << endl;
    }

    int line_num = 0;
    for (int i = 0; i < SingularityI.SEs.size(); i++)
    {
        line_num += SingularityI.SEs[i].vs_link.size();
    }

    ff << "LINES " << SingularityI.SEs.size() << " " << line_num + SingularityI.SEs.size() << endl;
    for (int i = 0; i < SingularityI.SEs.size(); i++)
    {
        ff << SingularityI.SEs[i].vs_link.size() << " ";
        for (int j = 0; j < SingularityI.SEs[i].vs_link.size(); j++)
            ff << SingularityI.SEs[i].vs_link[j] << " ";
        ff << endl;
    }

    ff << "POINT_DATA " << hex_mesh.HVs.size() << std::endl;
    ff << "SCALARS V_Scalars int" << std::endl;
    ff << "LOOKUP_TABLE V_Table" << std::endl;
    for (int i = 0; i < hex_mesh.HVs.size(); i++)
    {
        if (hex_mesh.HVs[i].where_location == 1)
        {
            bool found = false;
            for (int j = 0; j < SingularityI.SVs.size(); j++)
            {
                if (SingularityI.SVs[j].index_hex == i)
                    found = true;
            }

            bool found2 = false;
            for (int j = 0; j < SingularityI.SEs.size(); j++)
            {
                for (int k = 0; k < SingularityI.SEs[j].vs_link.size(); k++)
                    if (SingularityI.SEs[j].vs_link[k] == i)
                        found2 = true;
            }

            if (found)
                ff << "3" << std::endl; //SV
            else if (found2)
                ff << "2" << std::endl; //SE
            else
                ff << "1" << std::endl;
        }
        else
        {
            bool found = false;
            for (int j = 0; j < SingularityI.SVs.size(); j++)
            {
                if (SingularityI.SVs[j].index_hex == i)
                    found = true;
            }

            bool found2 = false;
            for (int j = 0; j < SingularityI.SEs.size(); j++)
            {
                for (int k = 0; k < SingularityI.SEs[j].vs_link.size(); k++)
                    if (SingularityI.SEs[j].vs_link[k] == i)
                        found2 = true;
            }

            if (found)
                ff << "3" << std::endl; //SV
            else if (found2)
                ff << "2" << std::endl; //SE
            else
                ff << "0" << std::endl;
        }
    }

    ff.close();
}

void h_io::save_singularG_VTK_Color_Nodes(const char *fname)
{
    std::ofstream ff(fname);

    ff << "# vtk DataFile Version 2.0" << std::endl << "Singularities Color" << std::endl;
    ff << "ASCII" << std::endl;
    ff << "DATASET POLYDATA" << std::endl;
    ff << "POINTS " << hex_mesh.HVs.size() << " double" << std::endl;

    for (int i = 0; i < hex_mesh.HVs.size(); i++)
        ff << hex_mesh.HVs[i].v[0] << " " << hex_mesh.HVs[i].v[1] << " " << hex_mesh.HVs[i].v[2] << std::endl;
    ff << "VERTICES " << SingularityI.SVs.size() << " " << SingularityI.SVs.size() * 2 << endl;
    for (int i = 0; i < SingularityI.SVs.size(); i++)
        ff << "1 " << SingularityI.SVs[i].index_hex << endl;
}

void h_io::save_singularG_VTK_Color_Edges(const char *fname)
{
    std::ofstream ff(fname);

    ff << "# vtk DataFile Version 2.0" << std::endl << "Singularities Color" << std::endl;
    ff << "ASCII" << std::endl;
    ff << "DATASET POLYDATA" << std::endl;

    ff << "POINTS " << hex_mesh.HVs.size() << " double" << std::endl;

    for (int i = 0; i < hex_mesh.HVs.size(); i++)
        ff << hex_mesh.HVs[i].v[0] << " " << hex_mesh.HVs[i].v[1] << " " << hex_mesh.HVs[i].v[2] << std::endl;

    int line_num = 0;
    for (int i = 0; i < SingularityI.SEs.size(); i++)
        line_num += SingularityI.SEs[i].vs_link.size();

    ff << "LINES " << SingularityI.SEs.size() << " " << line_num + SingularityI.SEs.size() << endl;
    for (int i = 0; i < SingularityI.SEs.size(); i++) {
        ff << SingularityI.SEs[i].vs_link.size() << " ";
        for (int j = 0; j < SingularityI.SEs[i].vs_link.size(); j++)
            ff << SingularityI.SEs[i].vs_link[j] << " ";
        ff << endl;
    }

    ff << "CELL_DATA " << line_num << std::endl;
    ff << "SCALARS SingularityEdge int" << std::endl;
    ff << "LOOKUP_TABLE SingularityEdge" << std::endl;
    for (int i = 0; i < SingularityI.SEs.size(); i++) ff << i << endl;
}

void h_io::save_singularNode_VTK(const char *fname)
{
    char file_meshlab_[300];
    std::fstream ff(fname, std::ios::out);

    ff << "# vtk DataFile Version 2.0" << std::endl << "hex mesh vtk data - converted from .off" << std::endl;
    ff << "ASCII" << std::endl;
    ff << "DATASET POLYDATA" << std::endl;

    ff << "POINTS " << hex_mesh.HVs.size() << " double" << std::endl;

    for (int i = 0; i < hex_mesh.HVs.size(); i++)
        ff << hex_mesh.HVs[i].v[0] << " " << hex_mesh.HVs[i].v[1] << " " << hex_mesh.HVs[i].v[2] << std::endl;

    int numOfSingularityNodes = 0;
    for (int i = 0; i < hex_mesh.HVs.size(); i++)
        for (int j = 0; j < SingularityI.SVs.size(); j++)
            if (SingularityI.SVs[j].index_hex == i)
                numOfSingularityNodes++;

    ff << "VERTICES " << numOfSingularityNodes << " " << numOfSingularityNodes * 2 << endl;
    for (int i = 0; i < hex_mesh.HVs.size(); i++)
        for (int j = 0; j < SingularityI.SVs.size(); j++)
            if (SingularityI.SVs[j].index_hex == i)
                ff << "1 " << SingularityI.SVs[j].index_hex << std::endl;

    ff.close();
}

void h_io::save_component_VTK(const char *fname)
{
    char file_meshlab_[300];
    std::fstream ff(fname, std::ios::out);

    ff << "# vtk DataFile Version 2.0" << std::endl << fname << std::endl;
    ff << "ASCII" << std::endl;
    ff << "DATASET UNSTRUCTURED_GRID" << std::endl;
    ff << "POINTS " << hex_mesh.HVs.size() << " float" << std::endl;

    for (int i = 0; i < hex_mesh.HVs.size(); i++)
        ff << hex_mesh.HVs[i].v[0] << " " << hex_mesh.HVs[i].v[1] << " " << hex_mesh.HVs[i].v[2] << std::endl;
    int cells_num = hex_mesh.HHs.size();
    ff << "CELLS " << cells_num << " " << cells_num *9 << endl;
    for (int i = 0; i < hex_mesh.HHs.size(); i++)
    {
        const Hex& cell = hex_mesh.HHs.at(i);
        ff << "8 ";
        for (int k = 0; k < 8; k++)
            ff << cell.V_Ids[k] << " ";
        ff << endl;
    }
    ff << "CELL_TYPES " << cells_num << endl;
    for (int i = 0; i < cells_num; i++)
        ff << "12" << endl;

    ff << "CELL_DATA " << cells_num << std::endl;
    ff << "SCALARS component int" << std::endl;
    ff << "LOOKUP_TABLE component_Table" << std::endl;
    for (int i = 0; i < cells_num; i++)
        ff << hex_mesh.HHs.at(i).Color_ID << endl;
    //////////////////////////////////////////////
    ff << "SCALARS raw_component_id int" << std::endl;
    ff << "LOOKUP_TABLE Default" << std::endl;
    std::vector<int> cells_label(cells_num);
    for (int i = 0; i < FrameI.FHs.size(); i++)
    {
        for (int j = 0; j < FrameI.FHs[i].hs_net.size(); j++)
        {
            const int cell_index = FrameI.FHs[i].hs_net.at(j);
            cells_label[cell_index] = i;
        }
    }
    for (int i = 0; i < cells_num; i++)
        ff << cells_label.at(i) << endl;
    ff.close();
}

void h_io::save_sheet_VTK(const char *fname, const unsigned long sheet_index)
{
    char file_meshlab_[300];
    std::fstream ff(fname, std::ios::out);

    ff << "# vtk DataFile Version 2.0" << std::endl << "hex mesh vtk data - converted from .off" << std::endl;
    ff << "ASCII" << std::endl;
    ff << "DATASET POLYDATA" << std::endl;

    ff << "POINTS " << hex_mesh.HVs.size() << " double" << std::endl;

    for (int i = 0; i < hex_mesh.HVs.size(); i++)
        ff << hex_mesh.HVs[i].v[0] << " " << hex_mesh.HVs[i].v[1] << " " << hex_mesh.HVs[i].v[2] << std::endl;

//    ff << "VERTICES " << SingularityI.SVs.size() << " " << SingularityI.SVs.size() * 2 << endl;
//    for (int i = 0; i < SingularityI.SVs.size(); i++)
//    {
//        ff << "1 " << SingularityI.SVs[i].index_hex << endl;
//    }
    ///////////////////////////////////////////////////////////////
    const Component_Sheet& cs = all_com_sheets.at(sheet_index);
    ff << "VERTICES " << cs.all_nodes.size() << " " << cs.all_nodes.size() * 2 << endl;
    for (int i = 0; i < cs.all_nodes.size(); i++)
    {
        ff << "1 " << FrameI.FVs.at(cs.all_nodes.at(i)).index_hex << std::endl;
    }
    int line_num = 0;
    for (int i = 0; i < cs.all_edges.size(); i++)
    {
        line_num += FrameI.FEs[cs.all_edges.at(i)].vs_link.size();
    }

    ff << "LINES " << cs.all_edges.size() << " " << line_num + cs.all_edges.size() << endl;
    for (int i = 0; i < cs.all_edges.size(); i++)
    {
        ff << FrameI.FEs[cs.all_edges.at(i)].vs_link.size() << " ";
        for (int j = 0; j < FrameI.FEs[cs.all_edges.at(i)].vs_link.size(); j++)
            ff << FrameI.FEs[cs.all_edges.at(i)].vs_link[j] << " ";
        ff << endl;
    }

    int cells_num = 0;
    for (int i = 0; i < cs.components.size(); i++)
    {
        cells_num += FrameI.FHs[cs.components.at(i)].hs_net.size();
    }
    ff << "CELLS " << cells_num << " " << cells_num *9 << endl;
    for (int i = 0; i < cs.components.size(); i++)
    {
        const int components_index = cs.components.at(i);
        for (int j = 0; j < FrameI.FHs[components_index].hs_net.size(); j++)
        {
            const int cell_index = FrameI.FHs[components_index].hs_net.at(j);
            const Hex& cell = hex_mesh.HHs.at(cell_index);
            ff << "8 ";
            for (int k = 0; k < 8; k++)
                ff << cell.V_Ids[k] << " ";
            ff << endl;
        }
    }

//    ff << "CELLS " << cs.components.size() << " " << 9 * cs.components.size() << endl;
//    for (int i = 0; i < cs.components.size(); i++)
//    {
//        ff << "8 ";
//        for (int j = 0; j < 8; j++)
//            ff << FrameI.FHs[cs.components.at(i)].FV_Ids[j] << " ";
//        ff << endl;
//    }
    ff << "CELL_TYPES " << cells_num << endl;
    for (int i = 0; i < cells_num; i++)
    {
        ff << "12" << endl;
    }
//    ff << "LINES " << SingularityI.SEs.size() << " " << 1 + SingularityI.SEs[sheet_index].vs_link.size() << endl;
//    for (int i = 0; i < 1; i++)
//    {
//        ff << SingularityI.SEs[i].vs_link.size() << " ";
//        for (int j = 0; j < SingularityI.SEs[i].vs_link.size(); j++)
//            ff << SingularityI.SEs[i].vs_link[j] << " ";
//        ff << endl;
//    }

    ///////////////////////////////////////////////////////////////
//    int line_num = 0;
//    for (int i = 0; i < SingularityI.SEs.size(); i++)
//    {
//        line_num += SingularityI.SEs[i].vs_link.size();
//    }
//
//    ff << "LINES " << SingularityI.SEs.size() << " " << line_num + SingularityI.SEs.size() << endl;
//    for (int i = 0; i < SingularityI.SEs.size(); i++)
//    {
//        ff << SingularityI.SEs[i].vs_link.size() << " ";
//        for (int j = 0; j < SingularityI.SEs[i].vs_link.size(); j++)
//            ff << SingularityI.SEs[i].vs_link[j] << " ";
//        ff << endl;
//    }
//
//    ff << "POINT_DATA " << hex_mesh.HVs.size() << std::endl;
//    ff << "SCALARS V_Scalars int" << std::endl;
//    ff << "LOOKUP_TABLE V_Table" << std::endl;
//    for (int i = 0; i < hex_mesh.HVs.size(); i++)
//    {
//        if (hex_mesh.HVs[i].where_location == 1)
//        {
//            bool found = false;
//            for (int j = 0; j < SingularityI.SVs.size(); j++)
//            {
//                if (SingularityI.SVs[j].index_hex == i)
//                    found = true;
//            }
//
//            bool found2 = false;
//            for (int j = 0; j < SingularityI.SEs.size(); j++)
//            {
//                for (int k = 0; k < SingularityI.SEs[j].vs_link.size(); k++)
//                    if (SingularityI.SEs[j].vs_link[k] == i)
//                        found2 = true;
//            }
//
//            if (found)
//                ff << "3" << std::endl; //SV
//            else if (found2)
//                ff << "2" << std::endl; //SE
//            else
//                ff << "1" << std::endl;
//        }
//        else
//        {
//            bool found = false;
//            for (int j = 0; j < SingularityI.SVs.size(); j++)
//            {
//                if (SingularityI.SVs[j].index_hex == i)
//                    found = true;
//            }
//
//            bool found2 = false;
//            for (int j = 0; j < SingularityI.SEs.size(); j++)
//            {
//                for (int k = 0; k < SingularityI.SEs[j].vs_link.size(); k++)
//                    if (SingularityI.SEs[j].vs_link[k] == i)
//                        found2 = true;
//            }
//
//            if (found)
//                ff << "3" << std::endl; //SV
//            else if (found2)
//                ff << "2" << std::endl; //SE
//            else
//                ff << "0" << std::endl;
//        }
//    }

    ff.close();
}

void h_io::save_sheet_hex_VTK(const char *fname, const unsigned long sheet_index)
{
    char file_meshlab_[300];
    std::fstream ff(fname, std::ios::out);

    ff << "# vtk DataFile Version 2.0" << std::endl << fname << std::endl;
    ff << "ASCII" << std::endl;
    ff << "DATASET UNSTRUCTURED_GRID" << std::endl;
    ff << "POINTS " << hex_mesh.HVs.size() << " float" << std::endl;

    for (int i = 0; i < hex_mesh.HVs.size(); i++)
        ff << hex_mesh.HVs[i].v[0] << " " << hex_mesh.HVs[i].v[1] << " " << hex_mesh.HVs[i].v[2] << std::endl;
    const Component_Sheet& cs = all_com_sheets.at(sheet_index);
    int cells_num = 0;
    for (int i = 0; i < cs.components.size(); i++)
        cells_num += FrameI.FHs[cs.components.at(i)].hs_net.size();
    ff << "CELLS " << cells_num << " " << cells_num *9 << endl;
    for (int i = 0; i < cs.components.size(); i++)
    {
        const int components_index = cs.components.at(i);
        for (int j = 0; j < FrameI.FHs[components_index].hs_net.size(); j++)
        {
            const int cell_index = FrameI.FHs[components_index].hs_net.at(j);
            const Hex& cell = hex_mesh.HHs.at(cell_index);
            ff << "8 ";
            for (int k = 0; k < 8; k++)
                ff << cell.V_Ids[k] << " ";
            ff << endl;
        }
    }
    ff << "CELL_TYPES " << cells_num << endl;
    for (int i = 0; i < cells_num; i++)
        ff << "12" << endl;
    //////////////////////////////////////////////
    ff << "CELL_DATA " << cells_num << std::endl;
    //////////////////////////////////////////////
    ff << "SCALARS component int" << std::endl;
    ff << "LOOKUP_TABLE component" << std::endl;
    for (int i = 0; i < cs.components.size(); i++)
    {
        const int components_index = cs.components.at(i);
        for (int j = 0; j < FrameI.FHs[components_index].hs_net.size(); j++)
        {
            const int cell_index = FrameI.FHs[components_index].hs_net.at(j);
            ff << hex_mesh.HHs.at(cell_index).Color_ID << endl;
        }
    }
    //////////////////////////////////////////////
    ff << "SCALARS raw_component_id int" << std::endl;
    ff << "LOOKUP_TABLE Default" << std::endl;
    for (int i = 0; i < cs.components.size(); i++)
    {
        const int components_index = cs.components.at(i);
        for (int j = 0; j < FrameI.FHs[components_index].hs_net.size(); j++)
        {
            ff << components_index << endl;
        }
    }
    //////////////////////////////////////////////
    ff << "SCALARS component_id int" << std::endl;
    ff << "LOOKUP_TABLE Default" << std::endl;
    int component_id = 0;
    for (int i = 0; i < cs.components.size(); i++)
    {
        const int components_index = cs.components.at(i);
        for (int j = 0; j < FrameI.FHs[components_index].hs_net.size(); j++)
        {
            ff << component_id << endl;
        }
        component_id++;
    }
    //////////////////////////////////////////////
    ff << "SCALARS sheet int" << std::endl;
    ff << "LOOKUP_TABLE Default" << std::endl;
    for (int i = 0; i < cells_num; i++)
        ff << sheet_index << endl;
    //////////////////////////////////////////////
    ff.close();
}

void h_io::save_sheet_geo_VTK(const char *fname, const unsigned long sheet_index)
{
    char file_meshlab_[300];
    std::fstream ff(fname, std::ios::out);

    ff << "# vtk DataFile Version 2.0" << std::endl << fname << std::endl;
    ff << "ASCII" << std::endl;
    ff << "DATASET POLYDATA" << std::endl;
    ff << "POINTS " << hex_mesh.HVs.size() << " double" << std::endl;
    for (int i = 0; i < hex_mesh.HVs.size(); i++)
        ff << hex_mesh.HVs[i].v[0] << " " << hex_mesh.HVs[i].v[1] << " " << hex_mesh.HVs[i].v[2] << std::endl;
    const Component_Sheet& cs = all_com_sheets.at(sheet_index);
    int vertices_num = cs.all_nodes.size();
    ff << "VERTICES " << cs.all_nodes.size() << " " << cs.all_nodes.size() * 2 << endl;
    for (int i = 0; i < cs.all_nodes.size(); i++)
        ff << "1 " << FrameI.FVs.at(cs.all_nodes.at(i)).index_hex << std::endl;
    int line_num = 0;
    for (int i = 0; i < cs.all_edges.size(); i++)
        line_num += FrameI.FEs[cs.all_edges.at(i)].vs_link.size();
    ff << "LINES " << cs.all_edges.size() << " " << line_num + cs.all_edges.size() << endl;
    for (int i = 0; i < cs.all_edges.size(); i++)
    {
        ff << FrameI.FEs[cs.all_edges.at(i)].vs_link.size() << " ";
        for (int j = 0; j < FrameI.FEs[cs.all_edges.at(i)].vs_link.size(); j++)
            ff << FrameI.FEs[cs.all_edges.at(i)].vs_link[j] << " ";
        ff << endl;
    }

    int polygon_num = FrameI.FFs[sheet_index].hfs_net_another.size();
//    polygon_num += cs.all_faces.size();
//    for (int i = 0; i < cs.all_faces.size(); i++)
//        polygon_num += FrameI.FFs[cs.all_faces.at(i)].hfs_net_another.size();
    ff << "POLYGONS " << polygon_num << " " << polygon_num * 5 << endl;
//    cout << "--------sheet " << sheet_index << endl;
//    cout << "-------- F size " << hex_mesh.HFs.size() << endl;
//    //for (int i = 0; i < all_com_sheets.size(); i++)
//    {
//        int i = sheet_index;
//        for (int j = 0; j < all_com_sheets[i].all_faces.size(); j++)
//            cout << all_com_sheets[i].all_faces[j] << " " ;
//        cout << endl;
//    }
//    cout << endl;
    //for (int i = 0; i < FrameI.FFs.size(); i++)
    {
        int i = sheet_index;
        for (int j = 0; j < FrameI.FFs[i].hfs_net_another.size(); j++)
        {
            ff << "4 " << hex_mesh.HFs[FrameI.FFs[i].hfs_net_another[j]].cv_Ids[0] << " ";
            ff << hex_mesh.HFs[FrameI.FFs[i].hfs_net_another[j]].cv_Ids[1] << " " << hex_mesh.HFs[FrameI.FFs[i].hfs_net_another[j]].cv_Ids[2] << " ";
            ff << hex_mesh.HFs[FrameI.FFs[i].hfs_net_another[j]].cv_Ids[3] << endl;
        }

//        int n = FrameI.FFs[cs.faces_right.at(i)].hfs_net_another.size();
        for (int j = 0; j < cs.all_faces.size(); j++)
        {

//            int id = FrameI.FFs[cs.faces_right[j]].index_own;
//            ff << "4 ";
//            ff << hex_mesh.HFs[id].cv_Ids[0] << " " << hex_mesh.HFs[id].cv_Ids[1] << " ";
//            ff << hex_mesh.HFs[id].cv_Ids[2] << " " << hex_mesh.HFs[id].cv_Ids[3] << endl;

//            int* fid = FrameI.FFs[cs.faces_right[j]].fv_Ids;
//            ff << "4 " << fid[0] << " " << fid[1] << " " << fid[2] << " " << fid[3] << endl;

//            int* fid = FrameI.FFs[cs.all_faces[j]].fv_Ids;
//            ff << "4 " << fid[0] << " " << fid[1] << " " << fid[2] << " " << fid[3] << endl;

//            ff << hex_mesh.HFs[id].cv_Ids[2] << " " << hex_mesh.HFs[id].cv_Ids[3] << endl;

//            int id = FrameI.FFs[cs.faces_right[j]].hfs_net_another;
//            ff << "4 ";
//            ff << hex_mesh.HFs[id].cv_Ids[0] << " " << hex_mesh.HFs[id].cv_Ids[1] << " ";
//            ff << hex_mesh.HFs[id].cv_Ids[2] << " " << hex_mesh.HFs[id].cv_Ids[3] << endl;
        }
    }
//    for (int i = 0; i < cs.all_faces.size(); i++)
//    {
//        for (int j = 0; j < FrameI.FFs[cs.all_faces.at(i)].hfs_net_another.size(); j++)
//        {
//            ff << "4 " << hex_mesh.HFs[FrameI.FFs[cs.all_faces.at(i)].hfs_net_another[j]].cv_Ids[0] << " ";
//            ff << hex_mesh.HFs[FrameI.FFs[cs.all_faces.at(i)].hfs_net_another[j]].cv_Ids[1] << " ";
//            ff << hex_mesh.HFs[FrameI.FFs[cs.all_faces.at(i)].hfs_net_another[j]].cv_Ids[2] << " ";
//            ff << hex_mesh.HFs[FrameI.FFs[cs.all_faces.at(i)].hfs_net_another[j]].cv_Ids[3] << endl;
//        }
//    }
//    ff << "CELL_DATA " << polygon_num << std::endl;
//    ff << "SCALARS sheet int 1" << std::endl;
//    ff << "LOOKUP_TABLE Default" << std::endl;
//    for (int j = 0; j < polygon_num; j++)
//        ff << sheet_index << std::endl;
    ff.close();
}
void h_io::save_sheets_VTK(const char *fname)
{
    char file_meshlab_[300];
    std::fstream ff(fname, std::ios::out);

    ff << "# vtk DataFile Version 2.0" << std::endl << fname << std::endl;
    ff << "ASCII" << std::endl;
    ff << "DATASET POLYDATA" << std::endl;
    ff << "POINTS " << hex_mesh.HVs.size() << " double" << std::endl;
    for (int i = 0; i < hex_mesh.HVs.size(); i++)
        ff << hex_mesh.HVs[i].v[0] << " " << hex_mesh.HVs[i].v[1] << " " << hex_mesh.HVs[i].v[2] << std::endl;

    int polygon_num = 0;
    for (int i = 0; i < all_com_sheets.size(); i++)
    {
        polygon_num += all_com_sheets[i].all_faces.size();
    }
    std::vector<int> cells_data(polygon_num);
    int count = 0;
    for (int i = 0; i < all_com_sheets.size(); i++)
    {
        for (int j = 0; j < all_com_sheets[i].all_faces.size(); j++)
        {
            cells_data.at(count++) = i;
        }
    }
    ff << "POLYGONS " << polygon_num << " " << polygon_num * 5 << endl;
    for (int i = 0; i < all_com_sheets.size(); i++)
    {
        //int i = sheet_index;
        for (int j = 0; j < all_com_sheets[i].all_faces.size(); j++)
        {
            int id = all_com_sheets[i].all_faces.at(j);
            ff << "4 " << hex_mesh.HFs[id].cv_Ids[0] << " ";
            ff << hex_mesh.HFs[id].cv_Ids[1] << " " << hex_mesh.HFs[id].cv_Ids[2] << " ";
            ff << hex_mesh.HFs[id].cv_Ids[3] << endl;
        }
    }

    ff << "POINT_DATA " << hex_mesh.HVs.size() << std::endl;
    ff << "SCALARS V_Scalars int" << std::endl;
    ff << "LOOKUP_TABLE V_Table" << std::endl;
    for (int i = 0; i < hex_mesh.HVs.size(); i++)
    {
        if (hex_mesh.HVs[i].where_location == 1)
        {
            bool found = false;
            for (int j = 0; j < FrameI.FVs.size(); j++)
            {
                if (FrameI.FVs[j].index_hex == i)
                    found = true;
            }

            bool found2 = false;
            for (int j = 0; j < FrameI.FEs.size(); j++)
            {
                for (int k = 0; k < FrameI.FEs[j].vs_link.size(); k++)
                    if (FrameI.FEs[j].vs_link[k] == i)
                        found2 = true;
            }

            if (found)
                ff << "3" << std::endl; //Frame_V
            else if (found2)
                ff << "2" << std::endl; //Frame_Edge
            else
                ff << "1" << std::endl;
        }
        else
        {
            bool found = false;
            for (int j = 0; j < FrameI.FVs.size(); j++)
            {
                if (FrameI.FVs[j].index_hex == i)
                    found = true;
            }

            bool found2 = false;
            for (int j = 0; j < FrameI.FEs.size(); j++)
            {
                for (int k = 0; k < FrameI.FEs[j].vs_link.size(); k++)
                    if (FrameI.FEs[j].vs_link[k] == i)
                        found2 = true;
            }

            if (found)
                ff << "3" << std::endl; //Frame_V
            else if (found2)
                ff << "2" << std::endl; //Frame_Edge
            else
                ff << "0" << std::endl;
        }
    }

    ff << "CELL_DATA " << polygon_num << std::endl;
    ff << "SCALARS sheet int 1" << std::endl;
    ff << "LOOKUP_TABLE Default" << std::endl;
    for (int j = 0; j < polygon_num; j++)
        ff << cells_data.at(j) << std::endl;
    ff.close();
}
h_io::~h_io(void)
{
}
