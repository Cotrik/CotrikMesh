/*
 *  This file is part of birch-clustering-algorithm.
 *
 *  birch-clustering-algorithm is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  birch-clustering-algorithm is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with birch-clustering-algorithm.  If not, see <http://www.gnu.org/licenses/>.
 *
 *	Copyright (C) 2011 Taesik Yoon (otterrrr@gmail.com)
 */


/** Simple test code for birch-clustering algorithm
 *
 * BIRCH has 4 phases: building, compacting, clustering, redistribution.
 * 
 * building - building cftree inserting a new data-point
 * compacting - make cftree smaller enlarging the range of sub-clusters
 * clustering - clustering sub-clusters(summarized clusters) using the existing clustering algorithm
 * redistribution - labeling data-points to the closest center
 */

#include "CFTree.h"
//#include "item_type.h"
#include <string>
#include <sstream>
#include <fstream>
#include <time.h>

#include "MeshFileReader.h"
#include "MeshFileWriter.h"

typedef CFTree<3u> cftree_type;
//struct item_type
//{
//	item_type() : id(0) { std::fill( item, item + sizeof(item)/sizeof(item[0]), 0 ); }
//	item_type( double* in_item ) : id(0) { std::copy(in_item, in_item+sizeof(item)/sizeof(item[0]), item); }
//	double& operator[]( int i ) { return item[i]; }
//	double operator[]( int i ) const { return item[i]; }
//	std::size_t size() const { return sizeof(item)/sizeof(item[0]); }
//
//	int& cid() { return id; }
//	const int cid() const { return id; }
//
//	double item[cftree_type::fdim];
//	int id;
//};

static double randf()
{
	return rand()/(double)RAND_MAX;
}

template<typename T>
static void print_items( const std::string fname, T& items )
{
	struct _compare_item_id
	{
		bool operator()( const item_type<3>& lhs, const item_type<3>& rhs ) const { return lhs.cid() < rhs.cid(); }
	};

	std::ofstream fout(fname.c_str());
	for( std::size_t i = 0 ; i < items.size() ; i++ )
	{
//		for( std::size_t d = 0 ; d < cftree_type::fdim ; d++ )
//			fout << items[i].item[d] << " ";
		fout << items[i].cid() << std::endl;
	}
	fout.close();
}

//template<boost::uint32_t dim>
//typedef std::vector<item_type<dim> > items_type;

template<boost::uint32_t dim>
static void load_items( const char* fname, std::vector<item_type<dim> >& items )
{
    if (fname) {
        std::ifstream fin(fname);
        std::string line;
        std::size_t cnt = 0;
        while (std::getline(fin, line))
            cnt++;
        items.reserve(cnt);

        fin.clear();
        fin.seekg(0);

        while (std::getline(fin, line)) {
            std::stringstream ss(line);

            cftree_type::item_vec_type item(cftree_type::fdim);
            for (std::size_t k = 0; k < item.size(); k++)
                ss >> item[k];

            items.push_back(&item[0]);
        }

        fin.close();
    } else {
        items.reserve(100000);
        for (std::size_t i = 0; i < items.capacity(); i++) {
            cftree_type::item_vec_type item(cftree_type::fdim);
            for (std::size_t k = 0; k < item.size(); k++)
                item[k] = randf();
            items.push_back(item_type<dim>(&item[0]));
        }
    }
}

template<boost::uint32_t dim>
static void load_items(const GeodesicDistance& gd, std::vector<item_type<dim> >& items)
{
    items.reserve(gd.surface.nVertices);
    for (size_t i = 0; i < gd.surface.nVertices; i++) {
        cftree_type::item_vec_type item(cftree_type::fdim);
        for (std::size_t k = 0; k < item.size(); k++)
            item[k] = gd.surface.vertices[3 * i + k];

        items.push_back(&item[0]);
    }
}


int main( int argc, char* argv[] )
{
    //GeodesicDistance& gd = CFEntry<3>::gd;
    GeodesicDistance gd;
    gd.inputMesh = argv[1];//"cap.faces.obj";
    gd.Init();
    std::cout << "GeodesicDistance(0, 100) = " << gd.GetGeodesicDistance(0, 100) << std::endl;

    std::vector<std::vector<double> > geodesic_distance_matrix(gd.surface.nVertices, std::vector<double>(gd.surface.nVertices, 0.0));
//#pragma omp parallel for
    for (size_t i = 0; i < gd.surface.nVertices; i++) {
        for (size_t j = i + 1; j < gd.surface.nVertices; j++) {
            //GeodesicDistance gd1;
            //gd1.inputMesh = argv[1];//"cap.faces.obj";
            //gd1.Init();
            geodesic_distance_matrix[i][j] = gd.GetGeodesicDistance(i, j);
            geodesic_distance_matrix[j][i] = geodesic_distance_matrix[i][j];
            //gd1.Destroy();
        }
        if (i % gd.surface.nVertices == 100)
            std::cout << "Finish Geodesic Distance Computation = " << double(i) * 100/ gd.surface.nVertices << "%\n";
    }
    std::cout << "Finish Geodesic Distance Computation\n";
    std::ofstream ofs("geodesic_distance.mat");
    for (size_t i = 0; i < gd.surface.nVertices; i++) {
        for (size_t j = 0; j < gd.surface.nVertices; j++) {
            ofs << geodesic_distance_matrix[i][j] << "\t";
        }
        ofs << "\n";
    }

	if( argc != 4 )
	{
		std::cout << "usage: birch (input-file) (range-threshold) (output-file)" << std::endl;
		return 0;
	}

	// load or generate items
	std::vector<item_type<3> > items;
	//load_items( argc >=2 ? argv[1] : NULL, items );
	load_items(gd, items);

	std::cout << items.size() << " items loaded" << std::endl;

	cftree_type::float_type birch_threshold = argc >=3 ? atof(argv[2]) : 0.25/(double)cftree_type::fdim;
	cftree_type tree(birch_threshold, 0);
	//cftree_type::gd = gd;
	// phase 1 and 2: building, compacting when overflows memory limit
	for( std::size_t i = 0 ; i < items.size() ; i++ )
		tree.insert( &items[i][0] );

	// phase 2 or 3: compacting? or clustering?
	// merging overlayed sub-clusters by rebuilding true
	tree.rebuild(false);

	// phase 3: clustering sub-clusters using the existing clustering algorithm
	//cftree_type::cfentry_vec_type entries;
	std::vector<CFEntry<3u> > entries;
	std::vector<int> cid_vec;
	tree.cluster( entries );

	// phase 4: redistribution

	// @comment ts - it is also possible to another clustering algorithm hereafter
	//				for example, we have k initial points for k-means clustering algorithm
	//tree.redist_kmeans( items, entries, 0 );

	std::vector<int> item_cids;
    tree.redist(items.begin(), items.end(), entries, item_cids);
    for (std::size_t i = 0; i < item_cids.size(); i++)
        items[i].cid() = item_cids[i];
    print_items(argc >= 4 ? argv[3] : "item_cid.txt", items);

    MeshFileReader reader(argv[1]);
    Mesh& mesh = (Mesh&)reader.GetMesh();
    MeshFileWriter writer(mesh, (std::string(argv[1]) + ".cluster.vtk").c_str());
    writer.WriteFile();
    writer.WritePointData(item_cids, "cluster");
	return 0;
}
