#include <string.h>
#include <time.h>
#include "io.h"
#include "frame_of_mesh.h"
#include "parameterization.h"
#include "parametric_optimization.h"
#include "component_class.h"
#include <assert.h>

char pathT[300] = "tri.off";
char path_IOH[300] = "hex.off";

//char Choices[300]="VOX";
char Choices[300] = "PO";
char Hex_NUM[300] = "8648";

h_io hio;

void parameter_optimize_a_hexmesh(char *para);
void only_re_parameterize_a_hexmesh(int min_Para);
int main(int argc, char* argv[])
{
    if (argc < 5)
    {
        std::cout << "Usage : opt <PO | P1 | P2 | VOX> <triangle_mesh.off> <out_hex.off> <hex_cell_elements>\n";
        return -1;
    }
	char patht[300] = "./";
	char pathioh[300] = "./";
	sprintf(Choices, "%s", argv[1]);
	sprintf(pathT, "%s", argv[2]);
	sprintf(path_IOH, "%s", argv[3]);
	sprintf(Hex_NUM,"%s",argv[4]);

	printf("reading triangle mesh...\n");
	hio.read_triangle_mesh_off(tri_mesh.Vs, tri_mesh.Ts, pathT);

	printf("reading hex mesh...\n");
	char path_temp[300];
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
//	hio.save_base_complex_VTK("base_complex.vtk");
//	hio.save_singularG_VTK("singularity.vtk");

	std::string::size_type sz;   // alias of size_t
	if (strcmp(Hex_NUM, " ") != 0)
		Para_Total_N = std::stoi(Hex_NUM, &sz);
	else
		Para_Total_N = hex_mesh.HHs.size();

	if (strcmp(Choices, "PO") == 0)		    parameter_optimize_a_hexmesh(Choices);
	else if (strcmp(Choices, "POF") == 0)	parameter_optimize_a_hexmesh(Choices);
	else if (strcmp(Choices, "P1") == 0)	only_re_parameterize_a_hexmesh(1);
	else if (strcmp(Choices, "P2") == 0)	only_re_parameterize_a_hexmesh(2);

	return 0;
}
void parameter_optimize_a_hexmesh(char *para)
{
	parametric_optimization po;
	char path_temp[300];
	sprintf(path_temp, "%s%s%s", "hex", Choices, "_Iter");
	if (strcmp(Choices, "PO") == 0)
		po.optimization_pipeline(path_temp);
	else if (strcmp(Choices, "POF") == 0)
		po.optimization_pipelineF(path_temp);
}
void only_re_parameterize_a_hexmesh(int min_Para)
{
	Para_min_N = min_Para;
	parameterization para;
	printf("re-meshing based on base-complex...\n");
	char path_temp[300];
	sprintf(path_temp, "%s%s%s", "hex", Choices, "_para");
	para.parameterization_main(path_temp);
}
