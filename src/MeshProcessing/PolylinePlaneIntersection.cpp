#include "PolylinePlaneIntersection.h"
#include <iostream>
using std::cout;
using std::endl;
PolylinePlaneIntersection::PolylinePlaneIntersection()
{
}

PolylinePlaneIntersection::~PolylinePlaneIntersection()
{
}

void PolylinePlaneIntersection::ComputeInterectionPoints_Polygon_ClosestOne(const std::vector<Mesh::Point>& _polygon,
	const Mesh::Point& _normal, const Mesh::Point& _point_on_plane, Mesh::Point& _intersect) const
{
	std::vector<Mesh::Point> intersects;
	ComputeInterectionPoints_Polygon(_polygon, _normal, _point_on_plane, intersects);

	//cout << "final intersects.size()" << intersects.size() << endl;

	double min_distance = FLT_MAX;
	for (auto& iIntersect : intersects)
	{
		double distance = (_point_on_plane - iIntersect).sqrnorm();
		if (min_distance > distance)//a better result is found
		{
			min_distance = distance;
			_intersect = iIntersect;

			//cout << "_intersect=\n" << _intersect << endl;
		}
	}

}

void PolylinePlaneIntersection::ComputeInterectionPoints_Polygon(const std::vector<Mesh::Point>& _polygon,
	const Mesh::Point& _normal, const Mesh::Point& _point_on_plane, std::vector<Mesh::Point>& _intersects) const
{
	if (_polygon.size() > 1)
	{
		// 1.
		// find jumping points, which is on the different side with the last point
		//

		std::vector<Mesh::Point>::const_iterator iter_begin = _polygon.begin();
		std::vector<Mesh::Point>::const_iterator iter_end = _polygon.end();
		std::vector<std::vector<Mesh::Point>::const_iterator> jumping_point;
		jumping_point.reserve(2);//generally there are 2 intersections
		do
		{
			// note iter_begin is updated here
			bool is_same_side = IsPointsStrictlyOnSameSide(iter_begin, iter_end, _normal, _point_on_plane);
			if (is_same_side == false)
			{
				jumping_point.push_back(iter_begin);
				assert(iter_begin != iter_end);
			}
		} while (iter_begin != iter_end);
		//cout << "jumping_point.size() " << jumping_point.size() << endl;

		//2.
		// get the intersection points
		//
		_intersects.reserve(2);
		for (auto iterJump = jumping_point.begin(); iterJump != jumping_point.end(); ++iterJump)
		{
			std::vector<Mesh::Point>::const_iterator iter_point = *iterJump;
			Mesh::Point end1 = *(iter_point - 1);
			Mesh::Point end2 = *iter_point;
			Mesh::Point intersect;
			LineSegmentPlaneIntersection(end1, end2, _normal, _point_on_plane, intersect);
			//cout << "LineSegmentPlaneIntersection intersect\n" << intersect << endl;
			_intersects.push_back(intersect);
		}

		// note for polygon, we need check if first and last point is on different side

		int side_first = ComputePointSide(_polygon[0], _normal, _point_on_plane);
		int side_last = ComputePointSide(_polygon.back(), _normal, _point_on_plane);
		if (side_first != side_last && side_first != 0)
		{
			Mesh::Point intersect;
			LineSegmentPlaneIntersection(_polygon.back(), _polygon[0], _normal, _point_on_plane, intersect);
			_intersects.push_back(intersect);
		}
		//cout << "_intersects.size() " << _intersects.size() << endl;


		//3.
		// clean the intersect points:
		// i.e. if two points are too close to each other, we merge them as one; 
		// this can happen if a vertex of the polygon is on the plane, then this point will be appended twice
		//
		for (auto iter = _intersects.begin(); iter != _intersects.end() - 1;)
		{
			if ((*iter - *(iter + 1)).norm() < 1e-5)//TODO: consider what is the best epsilon value to choose?
			{
				iter = _intersects.erase(iter);
			}
			else
			{
				++iter;
			}
		}

	}
	else
	{
		cout << "[WARNING FROM PolylinePlaneIntersection::ComputeInterectionPoints_Polygon] input polygon less than 2 points" << endl;
	}

}

int PolylinePlaneIntersection::ComputePointSide(const Mesh::Point& _point, 
	const Mesh::Point& _normal, const Mesh::Point& _point_on_plane) const
{
	Mesh::Point direction = _point - _point_on_plane;
	double inner_product = _normal.dot(direction);

	int side = 1;
	if (abs(inner_product) < 1e-10)
	{
		side = 0;
	}
	else if (inner_product < 0.0)
	{
		side = -1;
	}

	return side;
}

bool PolylinePlaneIntersection::IsPointsStrictlyOnSameSide(std::vector<Mesh::Point>::const_iterator& _points_begin,
	const std::vector<Mesh::Point>::const_iterator& _points_end,
	const Mesh::Point& _normal, const Mesh::Point& _point_on_plane) const
{
	// test if other points are on the same side with the first one
	int side1 = ComputePointSide(*_points_begin, _normal, _point_on_plane);
	++_points_begin;// traverse start from the second point
	while (_points_begin != _points_end)
	{ 
		int side = ComputePointSide(*_points_begin, _normal, _point_on_plane);
		if (side != side1)// find current point is not on the same side with the first one
		{
			return false;// now the points_begin points to the point that's firstly found side different with the first point
		}
		++_points_begin;
	}

	return true;
}

void PolylinePlaneIntersection::LineSegmentPlaneIntersection(const Mesh::Point& _end1, const Mesh::Point& _end2, 
	const Mesh::Point& _normal, const Mesh::Point& _point_on_plane, Mesh::Point& _intersect) const
{
	// get t s.t.
	// n * [t(P2-P1)+P1 - P] = 0
	// so t(P2-P1)+P1 is the intersection point
	double denominator = _normal.dot(_end2 - _end1);
	if (abs(denominator) < 1e-10)// the line is on the Plane, return one point
	{
		_intersect = _end1;
	}
	else
	{
		//double numerator = _normal.dot(_end1 - _point_on_plane);
		double numerator = _normal.dot(_point_on_plane - _end1);
		double t = numerator / denominator;
		_intersect = t * (_end2 - _end1) + _end1;
	}
}
