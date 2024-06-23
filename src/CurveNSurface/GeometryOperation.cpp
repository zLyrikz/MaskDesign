#include "GeometryOperation.h"
#include <Eigen/Geometry>

GeometryOperation::GeometryOperation()
{
}

void GeometryOperation::ProjectVector2Plane(const Eigen::Vector3d& _vector, const Eigen::Vector3d& _plane_normal, Eigen::Vector3d& _projection)
{
	_projection = _vector - 
		_plane_normal.dot(_vector) * _plane_normal;
}
