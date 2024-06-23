#pragma once
#include "CrossSectionOptimizeInfo.h"
#include "../Alglib/ap.h"

class CushionOptimize
{
public:
	CushionOptimize(CrossSectionOptimizeInfo* _cross_section_info);
	~CushionOptimize();

	void BfgsOptimize(CushionSurface& _optimal_cushion, int _max_iteration, double _epsilon_gradient);

	bool GetFEMVisualizationInfo(
		std::vector<double>& _area_distribution_color, std::vector<double>& _force_distribution_color, std::vector<double>& _sigmoid_color,
		TetrahedralMesh& _tetrahedral_cushion, Mesh& _deformed_cushion_surface,
		Mesh& _force_mesh, std::vector<double>& _node_force_magnitude, double _scale) const;

private:
	void ObjectiveFunctionAndGradient(const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr);

private:
	CrossSectionOptimizeInfo* cross_section_info_;

};

