#include "component_class.h"

component_class::component_class(void)
{

}
//component chords
void component_class::extract_all_component_sheets()
{
    all_com_sheets.clear();

    vector<vector<int> > All_es_areas;
    vector<bool> arrayE_test;
    for (int i = 0; i < FrameI.FEs.size(); i++)
    {
        arrayE_test.push_back(false);
    }

    while (true)
    {
        int start_eid = -1;
        for (int i = 0; i < arrayE_test.size(); i++)
        {
            if (!arrayE_test[i])
            {
                start_eid = i;
                break;
            }
        }
        if (start_eid == -1)
            break;

        All_es_areas.push_back(an_edge_area(start_eid, arrayE_test));
    }

    for (int i = 0; i < All_es_areas.size(); i++)
    {
        Component_Sheet cs;
        cs.id = i;
        cs.middle_edges = All_es_areas[i];

        build_component_information(cs.middle_edges, cs.edges_left, cs.edges_right, cs.middle_faces, cs.faces_left, cs.faces_right, cs.components,
                cs.nodes_left, cs.nodes_right);

        cs.all_edges = cs.edges_left;
        append_vector(cs.all_edges, cs.edges_right);
        append_vector(cs.all_edges, cs.middle_edges);

        cs.all_faces = cs.faces_left;
        append_vector(cs.all_faces, cs.faces_right);
        append_vector(cs.all_faces, cs.middle_faces);

        cs.all_nodes = cs.nodes_left;
        append_vector(cs.all_nodes, cs.nodes_right);

        all_com_sheets.push_back(cs);
    }
}
vector<int> component_class::an_edge_area(int e_id, vector<bool> &arrayE_test)
{
    vector<int> es, candiates, candidates_temp;
    candiates.push_back(e_id);
    bool stillhave = true;
    while (stillhave)
    {
        for (int i = 0; i < candiates.size(); i++)
        {
            e_id = candiates[i];
            if (arrayE_test[e_id])
                continue;
            es.push_back(e_id);

            arrayE_test[e_id] = true;

            for (int j = 0; j < FrameI.FEs[e_id].neighbor_Fs.size(); j++)
            {
                int nf_id = FrameI.FEs[e_id].neighbor_Fs[j];
                int which_e_ind;
                for (int k = 0; k < 4; k++)
                {
                    if (FrameI.FFs[nf_id].fe_Ids[k] == e_id)
                    {
                        which_e_ind = k;
                        break;
                    }
                }

                int op_eid = FrameI.FFs[nf_id].fe_Ids[(which_e_ind + 2) % 4];
                if (!arrayE_test[op_eid])
                    candidates_temp.push_back(op_eid);
            }
        }
        if (candidates_temp.size() == 0)
            stillhave = false;
        else
        {
            candiates.clear();
            candiates = candidates_temp;
            candidates_temp.clear();
        }
    }

    return es;
}
void component_class::build_component_information(vector<int> &edges_middle, vector<int> &edges_left, vector<int> &edges_right, vector<int> &faces_middle,
        vector<int> &faces_left, vector<int> &faces_right, vector<int> &components, vector<int> &nodes_left, vector<int> &nodes_right)
{
    for (int i = 0; i < edges_middle.size(); i++)
    {
        int fe = edges_middle[i];
        for (int j = 0; j < FrameI.FEs[fe].neighbor_Hs.size(); j++)
        {
            int h = FrameI.FEs[fe].neighbor_Hs[j];
            if (!insideVectorT(components, h))
                components.push_back(h);
        }
    }

    vector<bool> arrayes;
    initializeVectorT(arrayes, false, FrameI.FEs.size());
    for (int i = 0; i < edges_middle.size(); i++)
        arrayes[edges_middle[i]] = true;

    vector<int> faces;
    for (int i = 0; i < components.size(); i++)
    {
        for (int j = 0; j < 6; j++)
        {
            int f = FrameI.FHs[components[i]].neighbor_FS[j];
            bool notmiddle = true;
            for (int k = 0; k < 4; k++)
            {
                int eid = FrameI.FFs[f].fe_Ids[k];
                if (arrayes[eid])
                    notmiddle = false;
            }
            if (notmiddle)
            {
                bool already = insideVectorT(faces, f);
                if (!already)
                    faces.push_back(f);
            }
        }
    }

    vector<bool> arrayfs;
    initializeVectorT(arrayfs, false, FrameI.FFs.size());
    for (int i = 0; i < faces.size(); i++)
        arrayfs[faces[i]] = true;

    vector<int> candiatfs, candidatfs_temp;
    candiatfs.push_back(faces[0]);
    while (candiatfs.size() > 0)
    {
        for (int i = 0; i < candiatfs.size(); i++)
        {
            int f_id = candiatfs[i];
            if (insideVectorT(faces_left, f_id))
                continue;
            faces_left.push_back(f_id);

            for (int j = 0; j < 4; j++)
            {
                int e_id = FrameI.FFs[f_id].fe_Ids[j];
                for (int k = 0; k < FrameI.FEs[e_id].neighbor_Fs.size(); k++)
                {
                    int nf_id = FrameI.FEs[e_id].neighbor_Fs[k];
                    if (arrayfs[nf_id] && !insideVectorT(faces_left, nf_id))
                        candidatfs_temp.push_back(nf_id);
                }
            }
        }
        candiatfs.clear();
        candiatfs = candidatfs_temp;
        candidatfs_temp.clear();
    }

    set_exclusion(faces, faces_left, faces_right);
    vector<int> tempface2;
    for (int i = 0; i < faces_right.size(); i++)
        tempface2.push_back(-1);
    for (int i = 0; i < components.size(); i++)
    {
        int f1_id = -1, f2_trueid = -1;
        for (int j = 0; j < 6; j++)
        {
            int f = FrameI.FHs[components[i]].neighbor_FS[j];
            int whichf = set_contain(faces_left, f);
            if (whichf != -1)
            {
                f1_id = whichf;
            }
            whichf = set_contain(faces_right, f);
            if (whichf != -1)
            {
                f2_trueid = f;
            }
        }
        tempface2[f1_id] = f2_trueid;
    }
    faces_right = tempface2;

    for (int i = 0; i < faces_left.size(); i++)
    {
        for (int j = 0; j < 4; j++)
        {
            nodes_left.push_back(FrameI.FFs[faces_left[i]].fv_Ids[j]);
            edges_left.push_back(FrameI.FFs[faces_left[i]].fe_Ids[j]);
        }
    }
    set_redundent_clearn(nodes_left);
    set_redundent_clearn(edges_left);

    for (int i = 0; i < faces_right.size(); i++)
    {
        for (int j = 0; j < 4; j++)
        {
            nodes_right.push_back(FrameI.FFs[faces_right[i]].fv_Ids[j]);
            edges_right.push_back(FrameI.FFs[faces_right[i]].fe_Ids[j]);
        }
    }
    set_redundent_clearn(nodes_right);
    set_redundent_clearn(edges_right);

}

component_class::~component_class(void)
{
}
