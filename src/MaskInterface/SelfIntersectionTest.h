#pragma once
#include "DiscreteCushionPiece.h"
#include <Eigen/Core>
#include <vector>

class SelfIntersectionTest
{
public:
	SelfIntersectionTest();

	// output (begin,end) interval groups indiacting the self intersection parts
	// notice since the rail is periodic, the interval could also be periodic(begin index > end index)
	void TestIntersection(const std::vector<DiscreteCushionPiece>& _cross_sections, std::vector<std::pair<int, int>>& _intervals) const;
	bool IsIntersect(const DiscreteCushionPiece& _cross_section1, const DiscreteCushionPiece& _cross_section2) const;

	// take cross sections as circles with given a radius, check out the intersection state
	void TestCircleIntersection(float _radius, const std::vector<DiscreteCushionPiece>& _circles, std::vector<std::pair<int, int>>& _intervals) const;
	bool IsCircleIntersect(float _radius, const DiscreteCushionPiece& _circle1, const DiscreteCushionPiece& _circle2) const;

public:
	//utility functions

	// false iff cross section 2 completely in the positive side of 1
	bool CrossSection2Intersect1(const DiscreteCushionPiece& _cross_section1, const DiscreteCushionPiece& _cross_section2) const;
	// the special case when the cross sections are circles
	bool Circle2Intersect1(float _radius, const DiscreteCushionPiece& _circle1, const DiscreteCushionPiece& _circle2) const;

	// check if a point is in the positive side of a plane
	// true iff (point - plane_point)*plane_tangent >= 0
	bool IsInPositiveSide(const Eigen::Vector3d& _point, const Eigen::Vector3d& _plane_point, const Eigen::Vector3d& _plane_tangent) const;
};