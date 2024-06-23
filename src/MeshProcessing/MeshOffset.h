#pragma once
#include "../MeshViewer/MeshDefinition.h"

class MeshOffset
{
public:
	MeshOffset();
	~MeshOffset();

	// 
	// compute normal, then offset each vertex along the normal
	void OffsetSurface(const Mesh& _original, Mesh& _offset, double _distance) const;
	void OffsetSurface(const std::vector<std::vector<Mesh::Point>>& _vu_points, std::vector<std::vector<Mesh::Point>>& _offset_vu_points, double _distance) const;
	// connect original and the offset, get a thickened "watertight" mesh
	void ThickenSurface(const Mesh& _original, Mesh& _watertight, double _distance) const;
	void ThickenSurface_Libigl(const Mesh& _original, Mesh& _watertight, double _distance) const;
	// we offset the vu points then append the offset u points at the back of the original u points
	// this data structure is for silicone printing
	void ThickenSurface(const std::vector<std::vector<Mesh::Point>>& _vu_points, std::vector<std::vector<Mesh::Point>>& _merged_vu_points, double _distance) const;

	// we depend this on the natrual vertex index correspondence of the original and offset mesh
	void WatertightFromOffset(const Mesh& _original, const Mesh& _offset, Mesh& _watertight) const;
	
public:
	//Utility functions
	void FlipNormal(Mesh& _mesh) const;
	void ComputeMeshNormal(Mesh& _mesh) const;
};

