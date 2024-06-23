#pragma once
#include "../MeshViewer/MeshDefinition.h"

class InterfaceDesign
{
public:
	InterfaceDesign(const Mesh* _human_face, std::string _base_path);
	~InterfaceDesign();

	void Design(Mesh& _cushion_surface, Mesh& _mask_interface, int _max_cushion_iterations, 
		/*rail curve adjustment parameters*/double _objective, double _width, double _curvature, double _symmetry, double _align, double _angle,
		/*FEM simulation parameters*/ 	double _E_module, double _poisson_ratio, int _num_layer, double _thickness,
		/*cushion optimization paremeters*/double _weight_average_pressure, double _force_distribution, double _area_distribution,
		double _convexity, double _weight_tetrahedral_inversion, double _height_field_weight,
		/*cushion optimization regularization*/ double _weight_width, double _weight_length, double _weight_similarity, double _weight_size) const;

private:
	const Mesh* human_face_; // input human face for mask customization
	std::string base_path_;  // for aquiring stored data
	std::string data_path_;
};

