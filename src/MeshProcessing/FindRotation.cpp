#include "FindRotation.h"
#include <iostream>
using std::cout;
using std::endl;
using Eigen::Vector3d;
using Eigen::AngleAxisd;

FindRotation::FindRotation()
{
}

void FindRotation::VectorA2DirectionB_3D(const Eigen::Vector3d& _vector_a, const Eigen::Vector3d& _direction_b, Eigen::AngleAxisd& _rotation)
{
	_rotation = AngleAxisd(Eigen::Matrix3d::Identity());

	if (!_vector_a.isZero(1e-10) && !_direction_b.isZero(1e-10))
	{
		Vector3d normlized_a = _vector_a.normalized();
		Vector3d normlized_b = _direction_b.normalized();
		double inner_product = normlized_a.dot(normlized_b);

		double angle = acos(inner_product);
		Vector3d rotation_axis = normlized_a.cross(normlized_b);

		if (inner_product < 1.0 && inner_product > -1.0 && !rotation_axis.isZero(1e-10)) // else _vector_a already parallel to _direction_b
		{
			rotation_axis.normalize();
			_rotation = AngleAxisd(angle, rotation_axis);
		}
		else if (inner_product < 0)// parallel but in reverse direction
		{
			if (abs(_vector_a[0]) > 1e-10 || abs(_vector_a[1]) > 1e-10)
			{
				// a none zero vector that's perpendicular to vector a
				rotation_axis = Eigen::Vector3d(-_vector_a[1], _vector_a[0], 0);
				rotation_axis.normalize();

				_rotation = AngleAxisd(M_PI, rotation_axis);
			}
			else // only _vector_a[2] != 0 (parallel to axis z)
			{
				_rotation = AngleAxisd(M_PI, Eigen::Vector3d::UnitX());
			}
		}
	}
	else
	{
		cout << "[NOTE from FindRotation::VectorA2DirectionB_3D] _vector_a.isZero(1e-10) || _direction_b.isZero(1e-10)" << endl;
	}
}

void FindRotation::VectorA2DirectionB_3D(const Eigen::Vector3d& _vector_a, const Eigen::Vector3d& _direction_b, Eigen::Quaterniond& _rotation)
{
	_rotation = Eigen::Quaterniond(1.0, 0.0, 0.0, 0.0);

	if (!_vector_a.isZero(1e-10) && !_direction_b.isZero(1e-10))
	{
		Vector3d normlized_a = _vector_a.normalized();
		Vector3d normlized_b = _direction_b.normalized();
		double inner_product = normlized_a.dot(normlized_b);

		double angle = acos(inner_product);
		Vector3d rotation_axis = normlized_a.cross(normlized_b);

		if (inner_product < 1.0 && inner_product > -1.0 && !rotation_axis.isZero(1e-10)) // else _vector_a already parallel to _direction_b
		{
			rotation_axis.normalize();
			_rotation = AngleAxisd(angle, rotation_axis);
		}
		else if (inner_product < 0)// parallel but in reverse direction
		{
			if (abs(_vector_a[0]) > 1e-10 || abs(_vector_a[1]) > 1e-10)
			{
				// a none zero vector that's perpendicular to vector a
				rotation_axis = Eigen::Vector3d(-_vector_a[1], _vector_a[0], 0);
				rotation_axis.normalize();

				_rotation = AngleAxisd(M_PI, rotation_axis);
			}
			else // only _vector_a[2] != 0 (parallel to axis z)
			{
				_rotation = AngleAxisd(M_PI, Eigen::Vector3d::UnitX());
			}
		}
	}
	else
	{
		cout << "[NOTE from FindRotation::VectorA2DirectionB_3D] _vector_a.isZero(1e-10) || _direction_b.isZero(1e-10)" << endl;
	}
}
