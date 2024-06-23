#pragma once
#include "../MeshViewer/MeshDefinition.h"

// notice we don't handle the situation where a part of the polygon is on the plane!
// if a continuous sequence of vertex is on the plane, then this function skips all the vertex in the middle
class PolylinePlaneIntersection
{
public:
	PolylinePlaneIntersection();
	~PolylinePlaneIntersection();

	// get the point closest to the given point on plane
	// this is for polygon, i.e. first and last point of the input polyline connects
	void ComputeInterectionPoints_Polygon_ClosestOne(const std::vector<Mesh::Point>& _polygon, const Mesh::Point& _normal, const Mesh::Point& _point_on_plane,
		Mesh::Point& _intersect) const;

	void ComputeInterectionPoints_Polygon(const std::vector<Mesh::Point>& _polygon, const Mesh::Point& _normal, const Mesh::Point& _point_on_plane,
		std::vector<Mesh::Point>& _intersects) const;

public:
	// utility functions

	// -1 = negative side
	//  0 = on the plane
	//  1 = positive side
	// require input normal already normalized 
	int ComputePointSide(const Mesh::Point& _point, const Mesh::Point& _normal, const Mesh::Point& _point_on_plane) const;

	// return true iff points are on the same side
	// if false, the point_begin iterator points to the point that's firstly found side different with the first point;
	// if false, the output iterator will NOT be the same with the input begin or end
	bool IsPointsStrictlyOnSameSide(std::vector<Mesh::Point>::const_iterator& _points_begin,
		const std::vector<Mesh::Point>::const_iterator& _points_end,
		const Mesh::Point& _normal, const Mesh::Point& _point_on_plane) const;

	// input end points a line segment, output the intersection point of the line segment and the plane
	// require normal already normalized
	// return true if
	void LineSegmentPlaneIntersection(const Mesh::Point& _end1, const Mesh::Point& _end2,
		const Mesh::Point& _normal, const Mesh::Point& _point_on_plane, Mesh::Point& _intersect) const;
};

