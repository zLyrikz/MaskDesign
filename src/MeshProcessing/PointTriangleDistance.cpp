#include "PointTriangleDistance.h"
#include "../MeshProcessing/FindRotation.h"

#include <iostream>
using std::cout;
using std::endl;
using Eigen::Vector3d;
using Eigen::Vector3f;
using Eigen::Vector2d;
using Eigen::AngleAxisd;

PointTriangleDistance::PointTriangleDistance()
{
}

void PointTriangleDistance::ClosestPointOnTriangle(const Eigen::Vector3f& _point, 
	const Eigen::Vector3f& _a, const Eigen::Vector3f& _b, const Eigen::Vector3f& _c,
	Eigen::Vector3f& _closest_point)
{
	// Check if P in vertex region outside A
	Vector3f ab = _b - _a;
	Vector3f ac = _c - _a;
	Vector3f ap = _point - _a;
	float d1 = ab.dot(ap);
	float d2 = ac.dot(ap);
	if (d1 <= 0.0f && d2 <= 0.0f)
	{
		_closest_point = _a; // barycentric coordinates (1,0,0)
		return; 
	}

	// Check if P in vertex region outside B
	Vector3f bp = _point - _b;
	float d3 = ab.dot(bp);
	float d4 = ac.dot(bp);
	if (d3 >= 0.0f && d4 <= d3) // notice d4 - d3 = bc.dot(bp)
	{
		_closest_point = _b; // barycentric coordinates (0,1,0)
		return;
	}

	// Check if P in edge region of AB, if so return projection of P onto AB
	float vc = d1 * d4 - d3 * d2; // signed area of RAB (cross converted to dot)
	if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f)
	{
		float v = d1 / (d1 - d3);
		_closest_point = _a + v * ab;// barycentric coordinates (1-v,v,0)
		return; 
	}

	// Check if P in vertex region outside C
	Vector3f cp = _point - _c;
	float d5 = ab.dot(cp);
	float d6 = ac.dot(cp);
	if (d6 >= 0.0f && d5 <= d6)
	{
		_closest_point = _c; // barycentric coordinates (0,0,1)
		return; 
	}

	// Check if P in edge region of AC, if so return projection of P onto AC
	float vb = d5 * d2 - d1 * d6;
	if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) 
	{
		float w = d2 / (d2 - d6);
		_closest_point = _a + w * ac; // barycentric coordinates (1-w,0,w)
		return; 
	}

	// Check if P in edge region of BC, if so return projection of P onto BC
	float va = d3 * d6 - d5 * d4;
	if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f)
	{
		float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		_closest_point = _b + w * (_c - _b); // barycentric coordinates (0,1-w,w)
		return; 
	}

	// P inside face region. Compute Q through its barycentric coordinates (u,v,w)
	float denom = 1.0f / (va + vb + vc);
	float v = vb * denom;
	float w = vc * denom;
	_closest_point = _a + ab * v + ac * w; // = u*a + v*b + w*c, u = va * denom = 1.0f-v-w
	return; 
}

float PointTriangleDistance::Point2Triangle(const Eigen::Vector3f& _point, 
	const Eigen::Vector3f& _triangle_a, const Eigen::Vector3f& _triangle_b, const Eigen::Vector3f& _triangle_c, 
	const Eigen::Vector3f& _normal)
{
	Eigen::Vector3f closest_point;
	ClosestPointOnTriangle(_point, _triangle_a, _triangle_b, _triangle_c, closest_point);
	Vector3f vector = _point - closest_point;
	float distance = vector.norm();

	// decide sign
	bool negative = (vector.dot(_normal) < 0);
	if (negative)
	{
		distance = -distance;
	}

	return distance;
}

double PointTriangleDistance::Point2Triangle(const Eigen::Vector3d& _point,
	const Eigen::Vector3d& _triangle_a, const Eigen::Vector3d& _triangle_b, const Eigen::Vector3d& _triangle_c,
	bool _is_signed)
{
	/* distance1 = point project to the triangle plane
	* distance2 = projected point to triangle
	* point to triangle distance = sqrt(distance1 ^ 2 + distance2 ^ 2)
	*/ 


	// transform triangle to XOY plane
	AngleAxisd rotate;
	Vector3d translate;
	transformTriangle(_triangle_a, _triangle_b, _triangle_c, rotate, translate);

	// get the point and triangle position after transformation
	// notice the query point's z coordinate may not equal to 0
	Vector3d point = rotate * _point + translate;
	Vector3d triangle_a = rotate * _triangle_a + translate;
	Vector3d triangle_b = rotate * _triangle_b + translate;
	Vector3d triangle_c = rotate * _triangle_c + translate;

	//cout << "point\n " << point << endl;
	//cout << "triangle_a\n " << triangle_a << endl;
	//cout << "triangle_b\n " << triangle_b << endl;
	//cout << "triangle_c\n " << triangle_c << endl;

	// reduce the problem to 2D
	Vector2d point_xy = Vector2d(point.x(), point.y());
	Vector2d triangle_a_2d = Vector2d(triangle_a.x(), triangle_a.y());
	Vector2d triangle_b_2d = Vector2d(triangle_b.x(), triangle_b.y());
	Vector2d triangle_c_2d = Vector2d(triangle_c.x(), triangle_c.y());

	double distance1 = abs(point.z());
	double distance2 = Point2Triangle_2D(point_xy, triangle_a_2d, triangle_b_2d, triangle_c_2d);
	double distance = sqrt(distance1 * distance1 + distance2 * distance2);

	if (_is_signed)
	{
		if (point.z() < 0)
		{
			distance = -distance;
		}
	}

	return distance;
}

void PointTriangleDistance::transformTriangle(
	const Eigen::Vector3d& _triangle_a, const Eigen::Vector3d& _triangle_b, const Eigen::Vector3d& _triangle_c,
	Eigen::AngleAxisd& _rotation, Eigen::Vector3d& _translation)
{
	// 1 translate point A to origin
	Vector3d translate1 = -_triangle_a;

	// 2 rotate edge AB to have the same direction as the axis x
	AngleAxisd rotate2;

	FindRotation find_rotation;
	find_rotation.VectorA2DirectionB_3D(_triangle_b - _triangle_a, Eigen::Vector3d::UnitX(), rotate2);

	// 3 rotate along axis x, make point C lie in XOY plane
	AngleAxisd rotate3;

	Vector3d triangle_c_12 = rotate2 * (_triangle_c + translate1);// C point after transformation 1, 2
	Vector3d triangle_c_12_yz(0.0, triangle_c_12.y(), triangle_c_12.z());// y,z part of the point
	double c_12_yz_length = triangle_c_12_yz.norm();
	Vector3d triangle_c_123_yz = Vector3d(0.0, c_12_yz_length, 0.0);// y,z part of C point after transformation 1, 2, 3
	find_rotation.VectorA2DirectionB_3D(triangle_c_12_yz, triangle_c_123_yz, rotate3);

	// combine the transformations above
	// X -> R3*R2*(X + T1)
	// final rotation and translation: R = R3*R2; T = R3*R2*T1
	_rotation = rotate3 * rotate2;
	_translation = _rotation * translate1;
}

double PointTriangleDistance::Point2Triangle_2D(const Eigen::Vector2d& _point, 
	const Eigen::Vector2d& _triangle_a, const Eigen::Vector2d& _triangle_b, const Eigen::Vector2d& _triangle_c)
{
	double distance_ab1 = 0.0;
	double distance_bc1 = 0.0;
	double distance_ca1 = 0.0;
	double distance_ab2 = 0.0;
	double distance_bc2 = 0.0;
	double distance_ca2 = 0.0;
	Point2LineSegment_2D(_point, _triangle_a, _triangle_b, distance_ab1, distance_ab2);
	Point2LineSegment_2D(_point, _triangle_b, _triangle_c, distance_bc1, distance_bc2);
	Point2LineSegment_2D(_point, _triangle_c, _triangle_a, distance_ca1, distance_ca2);

	double final_distance = 0.0;

	if (!(
		((distance_ab1 > 0 && distance_bc1 > 0 && distance_ca1 > 0) || 
		(distance_ab1 < 0 && distance_bc1 < 0 && distance_ca1 < 0))
		)) // in this case, point is outside of the triangle
	{
		double distance_ab = sqrt(distance_ab1 * distance_ab1 + distance_ab2 * distance_ab2);
		double distance_bc = sqrt(distance_bc1 * distance_bc1 + distance_bc2 * distance_bc2);
		double distance_ca = sqrt(distance_ca1 * distance_ca1 + distance_ca2 * distance_ca2);

		// to find the minimum one
		final_distance = distance_ab;
		if (distance_bc < final_distance)
		{
			final_distance = distance_bc;
		}
		if (distance_ca < final_distance)
		{
			final_distance = distance_ca;
		}
	}

	return final_distance;
}

void PointTriangleDistance::Point2LineSegment_2D(
	const Eigen::Vector2d& _point, const Eigen::Vector2d& _line_start, const Eigen::Vector2d& _line_end, 
	double& _distance_point2line, double& _distance_project2segment)
{
	Vector2d line_direction = _line_end - _line_start;
	Vector2d left_normal(-line_direction.y(), line_direction.x());
	left_normal.normalize();

	Vector2d point_relative = _point - _line_start; // relative position of the point
	_distance_point2line = point_relative.dot(left_normal);

	// next, find projected point to line segment distance

	// relative projected point parallel to line_direction
	Vector2d projected_relative = point_relative - _distance_point2line * left_normal;// projected_relative + _line_start = projected point
	
	double ratio = 0.0; // projected_relative = ratio * line_direction
	int non_zero_element = 0; // indicate a non zero element of line_direction
	if (abs(line_direction[0]) > 1e-10)
	{
		non_zero_element = 0;
	}
	else if(abs(line_direction[1]) > 1e-10)
	{
		non_zero_element = 1;
	}
	else // in the input, two points of line segment are equal
	{
		_distance_project2segment = 0.0;
		return;
	}
	ratio = projected_relative[non_zero_element] / line_direction[non_zero_element];

	if (ratio < 0)
	{
		_distance_project2segment = projected_relative.norm();
	}
	else if (ratio > 1)
	{
		_distance_project2segment = projected_relative.norm() - line_direction.norm();
	}
	else
	{
		_distance_project2segment = 0.0;
	}
}




