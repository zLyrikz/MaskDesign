#include "CurveInfo.h"
#include <Eigen/Geometry>
#include <iostream>
using std::cout;
using std::endl;

CurveInfo::CurveInfo()
{
}

void CurveInfo::getFrenetFrame_3DCurve(const double& _t, const alglib::spline1dinterpolant* _curve,
	Eigen::Matrix4d& _frame)
{
	// the frame
	Eigen::Vector3d origin;
	Eigen::Vector3d tangent;
	Eigen::Vector3d normal;
	Eigen::Vector3d binormal;

	// get some derivatives of the spline
	Eigen::Vector3d second_derivative;
	for (int iDimension = 0; iDimension < 3; ++iDimension)
	{
		alglib::spline1ddiff(_curve[iDimension], _t,
			origin[iDimension], tangent[iDimension], second_derivative[iDimension]);
	}
	tangent.normalize();
	binormal = tangent.cross(second_derivative).normalized();
	normal = binormal.cross(tangent);

	//if (abs(normal.norm() - 1) > 1e-5)
	//{
	//	std::cout << "normal.norm()-1=" << normal.norm() - 1 << endl;
	//}

	_frame = Eigen::Matrix4d::Identity();
	_frame.block(0, 0, 3, 1) = tangent;
	_frame.block(0, 1, 3, 1) = normal;
	_frame.block(0, 2, 3, 1) = binormal;
	_frame.block(0, 3, 3, 1) = origin;
}

void CurveInfo::getFrameGivenNormal_3DCurve(const double& _t, const alglib::spline1dinterpolant* _curve,
	const Eigen::Vector3d& _normal, Eigen::Matrix4d& _frame)
{
	// the frame
	Eigen::Vector3d origin;
	Eigen::Vector3d tangent;
	Eigen::Vector3d normal = _normal.normalized();
	Eigen::Vector3d binormal;

	// get some derivatives of the spline
	Eigen::Vector3d second_derivative;
	for (int iDimension = 0; iDimension < 3; ++iDimension)
	{
		alglib::spline1ddiff(_curve[iDimension], _t,
			origin[iDimension], tangent[iDimension], second_derivative[iDimension]);
	}
	tangent.normalize();
	binormal = tangent.cross(normal).normalized();

	if (abs(tangent.dot(normal)) > 1e-5)
	{
		std::cout << "WARNING CurveInfo::getFrameGivenNormal_3DCurve -- tangent not perpendicular to the input normal" << endl;
	}

	_frame = Eigen::Matrix4d::Identity();
	_frame.block(0, 0, 3, 1) = tangent;
	_frame.block(0, 1, 3, 1) = normal;
	_frame.block(0, 2, 3, 1) = binormal;
	_frame.block(0, 3, 3, 1) = origin;
}

void CurveInfo::getFrameGivenNormal_3DCurve(const double& _t, const BSpline& _curve, const Eigen::Vector3d& _normal, Eigen::Matrix4d& _frame)
{
	// the frame
	Eigen::VectorXd origin;
	Eigen::VectorXd tangent;
	Eigen::Vector3d normal = _normal.normalized();
	Eigen::Vector3d binormal;

	_curve.getValue(_t, origin);
	_curve.getDerivative(_t, tangent);

	Eigen::Vector3d tangent_3d{ tangent[0], tangent [1], tangent[2] };
	tangent_3d.normalize();
	binormal = tangent_3d.cross(normal).normalized();

	if (abs(tangent_3d.dot(normal)) > 1e-5)
	{
		std::cout << "WARNING CurveInfo::getFrameGivenNormal_3DCurve -- tangent not perpendicular to the input normal" << endl;
	}

	_frame = Eigen::Matrix4d::Identity();
	_frame.block(0, 0, 3, 1) = tangent_3d;
	_frame.block(0, 1, 3, 1) = normal;
	_frame.block(0, 2, 3, 1) = binormal;
	_frame.block(0, 3, 3, 1) = origin;
}
