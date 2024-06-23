#pragma once
#include "../MeshViewer/MeshDefinition.h"


// notice here we always assume cushion surface is on the outside, and offset inward; normal of the cushion is pointing inward
class ConnectFrontNBack
{
public:
	ConnectFrontNBack();
	~ConnectFrontNBack();

	// offset the cushion boundary, to connect with the connector
	// notice,we need to make the boundary polyline sequence in same order
	void CreateIntermediatePart(const std::vector<Mesh::Point>& _back_inside_boundary,
		const std::vector<Mesh::Point>& _back_outside_boundary,
		const Mesh& _cushion_surface, double _distance, 
		Mesh& _intermediate, double _csg_tolerance = 0.2) const;
	// given the four boundarys, we connect them
	// notice,we need to make the boundary polyline sequence in same order
	void CreateIntermediatePart(const std::vector<Mesh::Point>& _back_inside_boundary,
		const std::vector<Mesh::Point>& _back_outside_boundary, 
		const std::vector<Mesh::Point>& _cushion_inside_boundary, const std::vector<Mesh::Point>& _cushion_outside_boundary,
		Mesh& _intermediate, double _csg_tolerance = 0.2) const;

	// reprensent the intermedite part by laryers of polyline, this is suitable for silicone printing
	// we only create polyline layers strictly between the cushion and connector boundarys (that means not including the boundary)
	// polyline layer order is from the cushion to the connector
	// 
	// NOT offset the cushion boundary to get the inside cushion boundary
	void ComputeIntermediatePartPolylineLayers(const std::vector<Mesh::Point>& _back_inside_boundary,
		const std::vector<Mesh::Point>& _back_outside_boundary,
		const Mesh& _cushion_surface, double _distance,
		std::vector<std::vector<Mesh::Point>>& _polyline_layers, int _num_layer = 10, bool _include_boundary = false) const;
	void ComputeIntermediatePartPolylineLayers(const std::vector<Mesh::Point>& _back_inside_boundary,
		const std::vector<Mesh::Point>& _back_outside_boundary,
		const std::vector<Mesh::Point>& _cushion_inside_boundary, const std::vector<Mesh::Point>& _cushion_outside_boundary,
		std::vector<std::vector<Mesh::Point>>& _polyline_layers, int _num_layer = 10, bool _include_boundary = false) const;

	//for each boundary vertex of cushion, shoot a -z direction ray.
	//  no intersection with back iff it's not overlaping with back
	// then translate the vertex along z direction ray to create overlap
	// cases when the trick fails: 
	//  1.vertex too far away from the back, so that z direction line does not have intersection with back
	//  2.vertex on the outside, but in the concave groove, it finds intersection using -z direction ray (can be addressed by computing signed distance with the back)
	void CreateOverlap(const Mesh& _back, Mesh& _cushion_surface) const;

	// offset the cushion surface to make it watertight, then connect with the back part
	void CreateWatertightOverlap(const Mesh& _back, const Mesh& _cushion_surface, double _distance, Mesh& _watertight) const;

	void CsgConnect(const Mesh& _back, const Mesh& _cushion_watertight, Mesh& _csg_connected) const;
	void CsgIntersect(const Mesh& _mesh1, const Mesh& _mesh2, Mesh& _csg) const;

	// all curves projected to XOY plane,
	// from target curve center, shoots a ray through target curve, trying to hit a point on the source curve
	// we approximate the hit point by two closest points on the source curve
	void GetBoundaryCorrespondenceInXOY(
		const std::vector<Mesh::Point>& _source_curve,
		const std::vector<Mesh::Point>& _target_curve, std::vector<Mesh::Point>& _correspond_source) const;

private:
	// utility functions

	void GetBoundaryCorrespondenceAndAddTolerance(
		const std::vector<Mesh::Point>& _back_inside_boundary, const std::vector<Mesh::Point>& _back_outside_boundary,
		std::vector<Mesh::Point>& _cushion_inside_boundary, std::vector<Mesh::Point>& _cushion_outside_boundary, 
		std::vector<Mesh::Point>& _correspond_back_outside, std::vector<Mesh::Point>& _correspond_back_inside, 
		double _csg_tolerance = 0.2) const;
	void GetCushionConnectingBoundary(const Mesh& _cushion_surface,
		std::vector<Mesh::Point>& _cushion_boundary, std::vector<Mesh::Point>& _offset_boundary, double _distance) const;
	void GetPointsPolarAngleInXOY(std::vector<double>& _angles, const std::vector<Mesh::Point>& _points, const OpenMesh::Vec2d& _origin) const;// note: input _angles must be empty; we do not clear it!

public:
	// utility functions

	// find the boundary with a larger z value(so this boundary should be connecting with the back part)
	// output an initial halfedge for iteration
	void FindConnectingBoundary(const Mesh& _cushion_surface, Mesh::HalfedgeHandle& _boundary) const;
	// output the whole boundary points
	void FindConnectingBoundary(const Mesh& _cushion_surface, std::vector<Mesh::VertexHandle>& _boundary) const;

	// this is for rendering
	void FlattenPolylineLayers(const std::vector<std::vector<Mesh::Point>>& _polyline_layers, Mesh& _flattened_polyline) const;

};

