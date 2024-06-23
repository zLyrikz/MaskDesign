#include "./MeshViewer/MeshDefinition.h"
#include "./MaskInterface/InterfaceDesign.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <Windows.h>
#include <filesystem>
using std::cout;
using std::endl;

void ReadParameters(const std::string& _file, std::vector<double>& _parameters)
{
	std::ifstream file(_file);
	if (file.is_open())
	{
		_parameters.clear();
		_parameters.reserve(21);
		std::string line;

		while (std::getline(file, line)) {
			std::istringstream iss(line);
			double value;
			std::string comment;
			if (!(iss >> value)) { continue; } // error or line without number
			_parameters.push_back(value);
		}

		file.close();
	}
	else
	{
		cout << "ERROR cannot open parameters.txt file: " << _file << endl;
	}
}

int main(int argc, char* argv[])
{
	// provide input instruction
	if (argc < 2) {
		std::cout << "Please provide the path to an aligned human face as an argument. Note that the face scale is measured in millimeters.\n";
		std::cout << "Usage: " << argv[0] << " [path_to_human_face_file]\n";
		std::cout << "Example: " << argv[0] << " D:/MaskDesign/data/face1.obj\n";
		return 1; // exit the program
	}

	// load an human face model (scale is measured in millimeters)
	// it needs to be already aligned with the generic mask
	std::string human_face_file(argv[1]);
	Mesh human_face;
	OpenMesh::IO::read_mesh(human_face, human_face_file);
	human_face.request_face_normals();
	human_face.request_vertex_normals();
	human_face.update_face_normals();
	human_face.update_vertex_normals();

	// find our working path
	char CurPath[MAX_PATH];
	GetModuleFileNameA(NULL, CurPath, 100);
	std::string path(CurPath);
	std::filesystem::path p = path;
	for (int i = 0; i < 3; i++) {
		p = p.parent_path();
	}
	std::string base_path = p.string();
	//cout << "base_path=" << base_path << endl;

	// parameters
	std::vector<double> parameters;
	ReadParameters(base_path + "/parameters.txt", parameters);
	int max_cushion_optimization_iteration = parameters[0];
	// (cushion initialization) trajectory curve adjustment parameters
	double objective = parameters[1];
	double width = parameters[2];
	double curvature = parameters[3];
	double symmetry = parameters[4];
	double align = parameters[5];
	double angle = parameters[6];
	//FEM simulation parameters
	double E_module = parameters[7];
	double poisson_ratio = parameters[8];
	int num_layer = parameters[9];
	double thickness = parameters[10];
	//cushion optimization paremeters 
	double weight_average_pressure = parameters[11];
	double force_distribution = parameters[12];
	double area_distribution = parameters[13];
	double weight_convexity = parameters[14];
	double weight_tetrahedral_inversion = parameters[15];
	double weight_height_field = parameters[16];
	// cushion optimization regularization parameters
	double weight_width = parameters[17];
	double weight_length = parameters[18];
	double weight_similarity = parameters[19];
	double weight_size = parameters[20];

	InterfaceDesign interface_design(&human_face, base_path);
	Mesh custom_cushion_surface;
	Mesh custom_mask_interface;
	interface_design.Design(custom_cushion_surface, custom_mask_interface, max_cushion_optimization_iteration,
		objective, width, curvature, symmetry, align, angle,
		E_module, poisson_ratio, num_layer, thickness,
		weight_average_pressure, force_distribution, area_distribution,
		weight_convexity, weight_tetrahedral_inversion, weight_height_field,
		weight_width, weight_length, weight_similarity, weight_size);

	// output designed custom cushion_surface and mask interface to obj file
	OpenMesh::IO::write_mesh(custom_cushion_surface, base_path + "/output/custom_cushion_surface.obj");
	OpenMesh::IO::write_mesh(custom_mask_interface, base_path + "/output/custom_mask_interface.obj");
	cout << "cushion surface and mask interface saved to " << base_path << "/output" << endl;
 
	
	return 0;
}
