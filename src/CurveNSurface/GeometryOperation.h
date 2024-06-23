#pragma once
#include <Eigen/Core>

class GeometryOperation
{
public:
	GeometryOperation();

	// mask sure the input plane normal is normalized
	void ProjectVector2Plane(const Eigen::Vector3d& _vector, const Eigen::Vector3d& _plane_normal, 
		Eigen::Vector3d& _projection);
};

