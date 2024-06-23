#pragma once
#include "../MeshViewer/MeshDefinition.h"
#include <Eigen/Geometry>

/* 
* for transforming mesh, given transformation matrix
* and some other mesh operations
*/

class TransformMesh
{
public:
	TransformMesh(Mesh& _mesh);

	// transform mesh in YOZ plane (2D) 
	// new Point = rotation * old Point + translation
	void TransformYOZ(const Eigen::Rotation2D<float>& _rotation, const Eigen::Vector2f& _translate);

	// input 4x4 rigid transformation matrix:
	// (R R R | T)
	// (R R R | T)
	// (R R R | T)
	// ----------
	// (0 0 0 | 1)
	void TransformMatrix4(const Eigen::Matrix4d& _transformation);

	void RotateTranslate(const Eigen::Matrix3d& _rotate, const Eigen::Vector3d& _translate);
private:
	Mesh& mesh_;
};

