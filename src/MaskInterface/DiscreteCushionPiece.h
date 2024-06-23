#pragma once
#include "../MeshViewer/MeshDefinition.h"
#include <Eigen/Geometry>
#include <vector>

// the frame and sampled cross section points at rail u
struct DiscreteCushionPiece
{
	// the standard frame
	Eigen::Vector3d frame_value_;
	Eigen::Vector3d frame_tangent_;
	Eigen::Vector3d frame_normal_;
	Eigen::Vector3d frame_binormal_;

	double u_ = 0.0;
	std::vector<Eigen::Vector3d> discrete_cross_section_;

	// the relaxed frame
	Eigen::Quaterniond rotation_;// the rotation to the standard frame tangent(and all other along)
	Eigen::Vector3d rotated_tangent_;//rotated_tangent_ = rotation_*frame_tangent_
	Eigen::Vector3d new_frame_normal_;
	Eigen::Vector3d new_frame_binormal_;

};