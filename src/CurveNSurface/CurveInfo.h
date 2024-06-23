#pragma once
#include "../Alglib/interpolation.h"
#include "../CurveNSurface/BSpline.h"
#include <Eigen/Core>

/*
* for calculate some curve information:
* Frenet Frame
*/

class CurveInfo
{
public:
	CurveInfo();

	// get the frenet frame of the curve at parameter t
	// (Tx Nx Bx | Ox)
	// (Ty Ny Bz | Oy)
	// (Tz Nz Bz | Oz)
	// ---------------
	// (0  0  0  | 1 )
	void getFrenetFrame_3DCurve(const double& _t, const alglib::spline1dinterpolant* _curve, Eigen::Matrix4d& _frame);

	void getFrameGivenNormal_3DCurve(const double& _t, const alglib::spline1dinterpolant* _curve, 
		const Eigen::Vector3d& _normal, Eigen::Matrix4d& _frame);// will normalize the normal in the function

	void getFrameGivenNormal_3DCurve(const double& _t, const BSpline& _curve,
		const Eigen::Vector3d& _normal, Eigen::Matrix4d& _frame);// will normalize the normal in the function
};

