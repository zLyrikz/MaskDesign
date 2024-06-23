#pragma once
#include <Eigen/Geometry>

class FindRotation
{
public:
	FindRotation();

	// rotate vector a to have the same direction as vector b
	// choose rotation axis to be the cross product of vector a and b
	// if any of the vector a,b is zero, set rotation = identity
	void VectorA2DirectionB_3D(const Eigen::Vector3d& _vector_a, const Eigen::Vector3d& _direction_b,
		Eigen::AngleAxisd& _rotation);
	// return type as Quaternion
	void VectorA2DirectionB_3D(const Eigen::Vector3d& _vector_a, const Eigen::Vector3d& _direction_b,
		Eigen::Quaterniond& _rotation);
};

