#include "SelfIntersectionTest.h"
#include <Eigen/Geometry>

SelfIntersectionTest::SelfIntersectionTest()
{

};

void SelfIntersectionTest::TestIntersection(const std::vector<DiscreteCushionPiece>& _cross_sections, std::vector<std::pair<int, int>>& _intervals) const
{
	_intervals.clear();
	_intervals.reserve(4);

	int num_cross_section = _cross_sections.size();
	{
		int begin_idx = 0;
		while ((begin_idx + 1) < num_cross_section)// find all parts of intersection
		{
			// find the begin index of the current part
			while ((begin_idx + 1) < num_cross_section
				&&
				!IsIntersect(_cross_sections[begin_idx], _cross_sections[begin_idx + 1]))
			{
				++begin_idx;
			}

			// find the end index of the current part
			if ((begin_idx + 1) < num_cross_section)
			{
				int counter = 1;
				while ((begin_idx + counter + 1) < num_cross_section
					&&
					IsIntersect(_cross_sections[begin_idx + counter], _cross_sections[begin_idx + counter + 1]))
				{
					++counter;
				}
				int end_idx = begin_idx + counter;// max possible value would be [num_cross_section - 1]

				_intervals.push_back(std::make_pair(begin_idx, end_idx));

				//prepare for the next part
				begin_idx = end_idx + 1;
			}
		}
	}

	// deal with two special cases (arise from the periodic problem)
	if (_intervals.size() >= 2)
	{
		if (_intervals[0].first == 0 && _intervals.back().second == num_cross_section - 1)
		{
			// check if first and last intersect
			if (IsIntersect(_cross_sections.back(), _cross_sections[0]))
			{
				// join the first and last parts
				_intervals[0].first = _intervals.back().first;
				_intervals.pop_back();
			}
		}
	}
	else if (_intervals.size() == 0)
	{
		// check if first and last intersect
		if (IsIntersect(_cross_sections.back(), _cross_sections[0]))
		{
			_intervals.push_back(std::make_pair(num_cross_section - 1, 0));
		}
	}
}

void SelfIntersectionTest::TestCircleIntersection(float _radius, const std::vector<DiscreteCushionPiece>& _circles, std::vector<std::pair<int, int>>& _intervals) const
{
	//same implementation with TestIntersection
	_intervals.clear();
	_intervals.reserve(4);

	int num_cross_section = _circles.size();
	{
		int begin_idx = 0;
		// find all parts of intersection
		while ((begin_idx + 1) < num_cross_section)
		{
			// find the begin index of the current part
			while ((begin_idx + 1) < num_cross_section
				&&
				!IsCircleIntersect(_radius, _circles[begin_idx], _circles[begin_idx + 1]))
			{
				++begin_idx;
			}

			// find the end index of the current part
			if ((begin_idx + 1) < num_cross_section)
			{
				int counter = 1;
				while ((begin_idx + counter + 1) < num_cross_section
					&&
					IsCircleIntersect(_radius, _circles[begin_idx + counter], _circles[begin_idx + counter + 1]))
				{
					++counter;
				}
				int end_idx = begin_idx + counter;// max possible value would be [num_cross_section - 1]

				_intervals.push_back(std::make_pair(begin_idx, end_idx));

				//prepare for the next part
				begin_idx = end_idx + 1;
			}
		}
	}

	// deal with two special cases (arise from the periodic problem)
	if (_intervals.size() >= 2)
	{
		if (_intervals[0].first == 0 && _intervals.back().second == num_cross_section - 1)
		{
			// check if first and last intersect
			if (IsCircleIntersect(_radius, _circles.back(), _circles[0]))
			{
				// join the first and last parts
				_intervals[0].first = _intervals.back().first;
				_intervals.pop_back();
			}
		}
	}
	else if (_intervals.size() == 0)
	{
		// check if first and last intersect
		if (IsCircleIntersect(_radius, _circles.back(), _circles[0]))
		{
			_intervals.push_back(std::make_pair(num_cross_section - 1, 0));
		}
	}
}


bool SelfIntersectionTest::IsIntersect(const DiscreteCushionPiece& _cross_section1, const DiscreteCushionPiece& _cross_section2) const
{
	// no intersect for both sides
	DiscreteCushionPiece cross_section2 = _cross_section2;
	cross_section2.rotated_tangent_ = -cross_section2.rotated_tangent_;
	return CrossSection2Intersect1(_cross_section1, _cross_section2) || CrossSection2Intersect1(cross_section2, _cross_section1);

	//return CrossSection2Intersect1(_cross_section1, _cross_section2);
}

bool SelfIntersectionTest::IsCircleIntersect(float _radius, const DiscreteCushionPiece& _circle1, const DiscreteCushionPiece& _circle2) const
{
	// no intersect for both sides
	DiscreteCushionPiece circle2 = _circle2;
	circle2.rotated_tangent_ = -circle2.rotated_tangent_;
	return Circle2Intersect1(_radius, _circle1, _circle2) || Circle2Intersect1(_radius, circle2, _circle1);
}

bool SelfIntersectionTest::CrossSection2Intersect1(const DiscreteCushionPiece& _cross_section1, const DiscreteCushionPiece& _cross_section2) const
{
	bool is_intersect = false;// will be true if find a point in cross section2 not on the positive side of the cross section1
	int num_points_cross_section2 = _cross_section2.discrete_cross_section_.size();

	// notice the trick here is to traverse cross section2 from both sides to the center (also not paying extra attention on the center point if num_points is odd)
	// because I guess there is a larger chance that they intersects closer to the sides
	for (int iPoint = 0; iPoint < num_points_cross_section2 / 2.0; ++iPoint)
	{
		bool not_intersect1 = IsInPositiveSide(_cross_section2.discrete_cross_section_[iPoint],
			_cross_section1.frame_value_, _cross_section1.rotated_tangent_);
		bool not_intersect2 = IsInPositiveSide(_cross_section2.discrete_cross_section_[num_points_cross_section2 - 1 - iPoint],
			_cross_section1.frame_value_, _cross_section1.rotated_tangent_);
		if (not_intersect1 == false || not_intersect2 == false)// current point of cross section2 in the negative side of cross section1 
		{
			is_intersect = true;
			break;
		}
	}
	return is_intersect;
}

bool SelfIntersectionTest::Circle2Intersect1(float _radius, const DiscreteCushionPiece& _circle1, const DiscreteCushionPiece& _circle2) const
{
	bool is_intersect = true;// will be false if all points in circle2 on the positive side of the circle1
	
	// only need to check if the furtherst points on circle2 is in the negative side of the circle1
	Eigen::Vector3d binormal = _circle1.rotated_tangent_.cross(_circle2.rotated_tangent_).normalized();
	Eigen::Vector3d normal2 = binormal.cross(_circle2.rotated_tangent_);
	Eigen::Vector3d furthest_point1 = _circle2.frame_value_ + _radius * normal2;
	Eigen::Vector3d furthest_point2 = _circle2.frame_value_ - _radius * normal2;

	if (IsInPositiveSide(furthest_point1, _circle1.frame_value_, _circle1.rotated_tangent_) &&
		IsInPositiveSide(furthest_point2, _circle1.frame_value_, _circle1.rotated_tangent_))
	{
		is_intersect = false;
	}

	return is_intersect;
}


bool SelfIntersectionTest::IsInPositiveSide(const Eigen::Vector3d& _point, const Eigen::Vector3d& _plane_point, const Eigen::Vector3d& _plane_tangent) const
{
	float inner_prodcut = (_point - _plane_point).dot(_plane_tangent);
	if (inner_prodcut < 0)
	{
		return false;
	}
	else
	{
		return true;
	}
}