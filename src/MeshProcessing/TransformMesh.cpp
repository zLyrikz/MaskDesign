#include "TransformMesh.h"
#include "MeshInfo.h"
#include "SubMesh.h"
#include <iostream>
#include <vector>


TransformMesh::TransformMesh(Mesh& _mesh)
	: mesh_(_mesh) // coding note: to initialize reference member mesh_, we have to do it here
{
}

void TransformMesh::TransformYOZ(const Eigen::Rotation2D<float>& _rotation, const Eigen::Vector2f& _translate)
{
	for (Mesh::VertexHandle iVertex : mesh_.vertices())
	{
		Mesh::Point point = mesh_.point(iVertex);
		Eigen::Vector2f point_yz(point[1], point[2]);

		point_yz = _rotation * point_yz + _translate;

		mesh_.set_point(iVertex, Mesh::Point(point[0], point_yz[0], point_yz[1]));
	}
}

void TransformMesh::TransformMatrix4(const Eigen::Matrix4d& _transformation)
{
	for (auto iterVertex = mesh_.vertices_begin(); iterVertex != mesh_.vertices_end(); ++iterVertex)
	{
		Mesh::Point new_point = mesh_.point(iterVertex);
		Eigen::Vector3d new_point_eigen = Eigen::Map<Eigen::Vector3d>(new_point.data());

		Eigen::Matrix3d rotate = _transformation.block<3, 3>(0, 0);
		Eigen::Vector3d translate = _transformation.block<3, 1>(0, 3);
		new_point_eigen = rotate * new_point_eigen + translate;

		new_point[0] = new_point_eigen[0];
		new_point[1] = new_point_eigen[1];
		new_point[2] = new_point_eigen[2];

		mesh_.set_point(*iterVertex, new_point);
	}
}

void TransformMesh::RotateTranslate(const Eigen::Matrix3d& _rotate, const Eigen::Vector3d& _translate)
{
	Eigen::Matrix4d transformation(Eigen::Matrix4d::Identity());
	transformation.block(0, 0, 3, 3) = _rotate;
	transformation.block(0, 3, 3, 1) = _translate;
	TransformMatrix4(transformation);
}

