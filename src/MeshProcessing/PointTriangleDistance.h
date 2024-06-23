#pragma once
#include <Eigen/Geometry> // Code note: how to avoid include it in header file


class PointTriangleDistance
{
public:
	PointTriangleDistance();

	// find closest point on triangle to the query point (very efficient way)
	// input triangle in counter clockwise
	// see P141 of http://staff.ustc.edu.cn/~zhuang/hpg/books/Morgan_Kaufmann.Real-Time%20Collision%20Detection(2005).pdf
	void ClosestPointOnTriangle(const Eigen::Vector3f& _point,
		const Eigen::Vector3f& _triangle_a, const Eigen::Vector3f& _triangle_b, const Eigen::Vector3f& _triangle_c,
		Eigen::Vector3f& _closest_point);

	// point to triangle distance (3D)
	// signed distance
	// input counter clockwise triangle vertices; distance < 0 if behind the triangle
	float Point2Triangle(const Eigen::Vector3f& _point,
		const Eigen::Vector3f& _triangle_a, const Eigen::Vector3f& _triangle_b, const Eigen::Vector3f& _triangle_c,
		const Eigen::Vector3f& _normal);

	// my previously implemented version
	double Point2Triangle(const Eigen::Vector3d& _point,
		const Eigen::Vector3d& _triangle_a, const Eigen::Vector3d& _triangle_b, const Eigen::Vector3d& _triangle_c,
		bool _is_signed = false);

	// transform triangle to XOY plane (z coordinate will be 0)
	// after transformation; A = origin, AB = positive x direction.
	// input triangle, output Roation and Translation: X -> Roation * X + Translation 
	void transformTriangle(const Eigen::Vector3d& _triangle_a, const Eigen::Vector3d& _triangle_b, const Eigen::Vector3d& _triangle_c,
		Eigen::AngleAxisd& _rotation, Eigen::Vector3d& _translation);

	// point to triangle distance (2D)
	double Point2Triangle_2D(const Eigen::Vector2d& _point,
		const Eigen::Vector2d& _triangle_a, const Eigen::Vector2d& _triangle_b, const Eigen::Vector2d& _triangle_c);

	// project the point to the line
	// output: 
	// 1 point to line distance (distance > 0 if point on the left of the line)
	// 2 projected point to line segment distance
	// note: the point to line segment distance = sqrt(distance1 ^ 2 + distance2 ^ 2)
	void Point2LineSegment_2D(const Eigen::Vector2d& _point,
		const Eigen::Vector2d& _line_start, const Eigen::Vector2d& _line_end,
		double& _distance_point2line, double& _distance_project2segment);

	//// although given points in 2D space, assume here point lies on the line 
	//double Point2LineSegment_1D(const Eigen::Vector2d& _point,
	//	const Eigen::Vector2d& _line_a, const Eigen::Vector2d& _line_b);

	// point to plane distance 
	// input plane by giving three noncollinear points on the plane
	void Point2Plane(const Eigen::Vector3d& _point,
		const Eigen::Vector3d& _plane_a, const Eigen::Vector3d& _plane_b, const Eigen::Vector3d& _plane_c,
		double& _distance, Eigen::Vector3d& _closest_point);
};

