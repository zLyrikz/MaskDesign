#pragma once
#include "../MeshViewer/MeshDefinition.h"

/* compute intersection of a triangle mesh and a plane
	plane: Ax + By + Cz = D

*/
class MeshPlaneIntersection
{
public:
	MeshPlaneIntersection();
	void setData(const Mesh* _mesh, const Mesh::Point& _normal, double _D);
	void setData(const Mesh* _mesh, const Mesh::Point& _normal, const Mesh::Point& _point_on_plane);

	void FindIntersectBoundaryFaces();

	// store the intersected faces in their adjacent order
	// start from a interescted boundary face, end at the boundary (so only find one countinous part of the faces, not all parts)
	// search neighbors and neighbor's neighbors that has intersection with the plane (like Breadth-first search)
	// I think this only works, when the intersection is a curve and the intersected mesh faces are not too sparse to be like a curve
	void FindIntersectFaces();

	// this is an extension of the original FindIntersectFaces() function
	void FindIntersectFaces_AllParts();

public:
	// TODO: use a new data structure, store two vertex on same side, one vertex on the other side, one sign to tell the vertex side.
	// return false iff can't hit end
	bool FindLocalIntersect(const Mesh::FaceHandle& _start_face, const Mesh::FaceHandle& _end_face,
		const Mesh::Point& _direction, std::vector<Mesh::Point>& _intersect_polyline) const;
	// end until the intersect polyline reach a given length
	// start point given must be on the start face, it's for computing the polyline length
	bool FindLocalIntersect(const Mesh::FaceHandle& _start_face, const Mesh::Point& _start_point, double _end_length,
		const Mesh::Point& _direction, std::vector<Mesh::Point>& _intersect_polyline, double& _actual_length) const;
public:
	// get the intersection mesh faces that's also on the mesh boundary
	void getIntersectBoundaryFaces(std::vector<Mesh::FaceHandle>& _intersect_boundary_face) const;
	void getIntersectFaces(std::vector<Mesh::FaceHandle>& _intersect_face) const;

	// one point for each intersect face
	void getIntersectPoints(std::vector<Mesh::Point>& _points) const;

	// from all the parts of intersection, get the one closet to a given point
	void getIntersectPoints_ClosetToAPoint(const Mesh::Point& _point, std::vector<Mesh::Point>& _points) const;
	void getIntersectPoints_AllParts(std::vector<std::vector<Mesh::Point>>& _points) const;


public:
	// previous interface
	/*
	* notice! in the new version, FindIntersectFaces has done the work of FindIntersectPolyline
	* so FindIntersectPolyline works exactly as function FindIntersectFaces now
	*/
	// according to the intersected faces found
	// take one step further to find the intersect points on edges of each face, so they should form a polyline
	void FindIntersectPolyline();
	// this is an extension of the original FindIntersectPolyline() function
	void FindIntersectPolyline_AllParts();

private:
	// record intersection information
	// note that for each intersected face, there are two intersected edges
	struct IntersectFace
	{
		Mesh::FaceHandle face_;

		// the halfedges of this face that intersects with the plane
		// and the corresponding intersect points on halfedges
		std::vector<std::pair<OpenMesh::SmartHalfedgeHandle, Mesh::Point>> halfedges_points_;
	};

private:
	void Init();

	// test if a face has intersection with the plane
	bool OneFaceIntersection(Mesh::FaceHandle _face) const; 

	void FindIntersectFaces(const std::vector<Mesh::FaceHandle>& _intersect_boundary_face, std::vector<IntersectFace>& _intersect_elements) const;
	void FindIntersectPolyline(std::vector<IntersectFace>& _intersect_elements) const;

	// one point for each intersect face
	void getIntersectPoints(const std::vector<IntersectFace>& _intersect_elements, std::vector<Mesh::Point>& _points) const;

	bool FindOneIntersectFaceExceptSome(const std::vector<IntersectFace>& _intersect_elements, Mesh::FaceHandle& _face) const;

	bool HalfedgeIntersectPlane(const OpenMesh::SmartHalfedgeHandle& _halfedge, const Mesh::Point& _normal, double _D) const;
	void LineIntersectPlane(const Mesh::Point& _start, const Mesh::Point& _end, const Mesh::Point& _normal, double _D, Mesh::Point& _intersect) const;
private:
	// used in the function FindLcoalIntersect
	// given the intersect element of the start face, the halfedge in the vector has same order with the search direction,
	// search on until hit the end face,
	// iff cannot hit the end face, return false
	bool FindLocalIntersect_WithStartInfo(const Mesh::FaceHandle& _start_face, const Mesh::FaceHandle& _end_face,
		std::vector<IntersectFace>& _intersect_element) const;
	bool FindLocalIntersect_WithStartInfo(const Mesh::FaceHandle& _start_face, const Mesh::Point& _start_point, double _end_length,
		std::vector<IntersectFace>& _intersect_element, double& _actual_length) const;
	// used in the function FindLcoalIntersect
	void PrepareStartFace(const Mesh::FaceHandle& _start_face, const Mesh::Point& _direction, std::vector<IntersectFace>& _intersect_element) const;
	// return false if next face is invalid
	// if is valid, store the intersect information in the _intersect_element and update the current face and halfedge to the new face
	bool SearchNextFace(Mesh::FaceHandle& _current_face, OpenMesh::SmartHalfedgeHandle& _current_halfedge, std::vector<IntersectFace>& _intersect_element) const;
	double ComputeIntersectLineSegmentLength(const IntersectFace& _face) const;
	void getIntersectPolyline(const std::vector<IntersectFace>& _intersect_elements, std::vector<Mesh::Point>& _polyline) const;
private:
	const Mesh* mesh_;

	// plane: Ax + By + Cz = D
	Mesh::Point normal_; // (A, B, C)
	double D_; // plane parameter D

	std::vector<IntersectFace> intersect_elements_;
	std::vector<Mesh::FaceHandle> intersect_boundary_face_;// the intersection mesh faces that's also on the mesh boundary
	//std::vector<std::array<Mesh::FaceHandle, 2>> boundary_pair_;// TODO: find all the intersect parts!

	std::vector<std::vector<IntersectFace>> intersect_elements_all_parts_;// find all discontinuous parts of the intersection elements
};

