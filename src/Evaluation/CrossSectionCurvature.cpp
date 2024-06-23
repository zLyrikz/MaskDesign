#include "CrossSectionCurvature.h"
#include <iostream>
#include <functional>
using std::cout;
using std::endl;
CrossSectionCurvature::CrossSectionCurvature()
{
	cross_section_ = nullptr;
	half1_triangle_area_ = 0.0;
	half2_triangle_area_ = 0.0;
}

CrossSectionCurvature::~CrossSectionCurvature()
{
}

void CrossSectionCurvature::ComputeMaxCurvatureCandidates(const CrossSection* _cross_section
	, double& _half1_start_max_curvature, double& _half1_mid_max_curvature, double& _half1_end_max_curvature
	, double& _half2_start_max_curvature, double& _half2_mid_max_curvature, double& _half2_end_max_curvature)
{
	SetCrossSection(_cross_section);

	// half1
	half1_triangle_area_ = ComputeControlPolygonArea(-(*_cross_section->GetHalf1Vector01()), *_cross_section->GetHalf1Vector12());
	Eigen::Vector2d half1_midpiont = 0.5 * (*_cross_section->GetHalf1Point2() + *_cross_section->GetPoint0());// the midpoint between point0 and point2
	double half1_length_point1_to_mid = (half1_midpiont - *_cross_section->GetHalf1Point1()).norm();
	_half1_mid_max_curvature = (half1_length_point1_to_mid * half1_length_point1_to_mid * half1_length_point1_to_mid) / (half1_triangle_area_ * half1_triangle_area_);
	//cout << "half1_length_point1_to_mid=" << half1_length_point1_to_mid << endl;
	//cout << "_half1_mid_max_curvature=" << _half1_mid_max_curvature << endl;
	//cout << "_half1_triangle_area_=" << half1_triangle_area_ << endl;
	_half1_start_max_curvature = abs(half1_triangle_area_) / 
		(_cross_section->GetHalf1Vector01Length() * _cross_section->GetHalf1Vector01Length() * _cross_section->GetHalf1Vector01Length());
	_half1_end_max_curvature = abs(half1_triangle_area_) /
		(_cross_section->GetHalf1Vector12Length() * _cross_section->GetHalf1Vector12Length() * _cross_section->GetHalf1Vector12Length());

	//half2
	half2_triangle_area_ = ComputeControlPolygonArea(*_cross_section->GetHalf2Vector12(), -(*_cross_section->GetHalf2Vector01()));
	Eigen::Vector2d half2_midpiont = 0.5 * (*_cross_section->GetHalf2Point2() + *_cross_section->GetPoint0());// the midpoint between point0 and point2
	double half2_length_point1_to_mid = (half2_midpiont - *_cross_section->GetHalf2Point1()).norm();
	_half2_mid_max_curvature = (half2_length_point1_to_mid * half2_length_point1_to_mid * half2_length_point1_to_mid) / (half2_triangle_area_ * half2_triangle_area_);
	_half2_start_max_curvature = abs(half2_triangle_area_) /
		(_cross_section->GetHalf2Vector01Length() * _cross_section->GetHalf2Vector01Length() * _cross_section->GetHalf2Vector01Length());
	_half2_end_max_curvature = abs(half2_triangle_area_) /
		(_cross_section->GetHalf2Vector12Length() * _cross_section->GetHalf2Vector12Length() * _cross_section->GetHalf2Vector12Length());

}

void CrossSectionCurvature::ComputeMaxCurvature(const CrossSection* _cross_section, double& _half1, double& _half2)
{
	SetCrossSection(_cross_section);

	// half1
	half1_triangle_area_ = ComputeControlPolygonArea(-(*_cross_section->GetHalf1Vector01()), *_cross_section->GetHalf1Vector12());
	Eigen::Vector2d half1_midpiont = 0.5 * (*_cross_section->GetHalf1Point2() + *_cross_section->GetPoint0());// the midpoint between point0 and point2
	Eigen::Vector2d half1_mid_midpiont_point0 = 0.5 * (half1_midpiont + *_cross_section->GetPoint0());// the midpoint between point0 and midpoint
	Eigen::Vector2d half1_mid_midpiont_point2 = 0.5 * (*_cross_section->GetHalf1Point2() + half1_midpiont);// the midpoint between midpoint and point2
	double radius = (half1_midpiont - *_cross_section->GetPoint0()).norm() / 2.0;
	double distance_point1_midmid0 = (*_cross_section->GetHalf1Point1() - half1_mid_midpiont_point0).norm();
	double distance_point1_midmid2 = (*_cross_section->GetHalf1Point1() - half1_mid_midpiont_point2).norm();
	if (distance_point1_midmid0 > radius && distance_point1_midmid2 > radius)
	{
		double half1_length_point1_to_mid = (half1_midpiont - *_cross_section->GetHalf1Point1()).norm();
		_half1 = (half1_length_point1_to_mid * half1_length_point1_to_mid * half1_length_point1_to_mid) / (half1_triangle_area_ * half1_triangle_area_);
	}
	else
	{
		double half1_start_max_curvature = abs(half1_triangle_area_) /
			(_cross_section->GetHalf1Vector01Length() * _cross_section->GetHalf1Vector01Length() * _cross_section->GetHalf1Vector01Length());
		double half1_end_max_curvature = abs(half1_triangle_area_) /
			(_cross_section->GetHalf1Vector12Length() * _cross_section->GetHalf1Vector12Length() * _cross_section->GetHalf1Vector12Length());
		if (half1_start_max_curvature > half1_end_max_curvature)
		{
			_half1 = half1_start_max_curvature;
		}
		else
		{
			_half1 = half1_end_max_curvature;
		}
	}


	//half2
	half2_triangle_area_ = ComputeControlPolygonArea(*_cross_section->GetHalf2Vector12(), -(*_cross_section->GetHalf2Vector01()));
	Eigen::Vector2d half2_midpiont = 0.5 * (*_cross_section->GetHalf2Point2() + *_cross_section->GetPoint0());// the midpoint between point0 and point2
	Eigen::Vector2d half2_mid_midpiont_point0 = 0.5 * (half2_midpiont + *_cross_section->GetPoint0());// the midpoint between point0 and midpoint
	Eigen::Vector2d half2_mid_midpiont_point2 = 0.5 * (*_cross_section->GetHalf2Point2() + half2_midpiont);// the midpoint between midpoint and point2
	radius = (half2_midpiont - *_cross_section->GetPoint0()).norm() / 2.0;
	distance_point1_midmid0 = (*_cross_section->GetHalf2Point1() - half2_mid_midpiont_point0).norm();
	distance_point1_midmid2 = (*_cross_section->GetHalf2Point1() - half2_mid_midpiont_point2).norm();
	if (distance_point1_midmid0 > radius && distance_point1_midmid2 > radius)
	{
		double half2_length_point1_to_mid = (half2_midpiont - *_cross_section->GetHalf2Point1()).norm();
		_half2 = (half2_length_point1_to_mid * half2_length_point1_to_mid * half2_length_point1_to_mid) / (half2_triangle_area_ * half2_triangle_area_);

	}
	else
	{
		double half2_start_max_curvature = abs(half2_triangle_area_) /
			(_cross_section->GetHalf2Vector01Length() * _cross_section->GetHalf2Vector01Length() * _cross_section->GetHalf2Vector01Length());
		double half2_end_max_curvature = abs(half2_triangle_area_) /
			(_cross_section->GetHalf2Vector12Length() * _cross_section->GetHalf2Vector12Length() * _cross_section->GetHalf2Vector12Length());
		if (half2_start_max_curvature > half2_end_max_curvature)
		{
			_half2 = half2_start_max_curvature;
		}
		else
		{
			_half2 = half2_end_max_curvature;
		}
	}
}

void CrossSectionCurvature::ComputeMaxCurvatureCandidates_WithGradient(const CrossSection* _cross_section, 
	double& _half1_start_max_curvature, double& _half1_mid_max_curvature, double& _half1_end_max_curvature, 
	double& _half2_start_max_curvature, double& _half2_mid_max_curvature, double& _half2_end_max_curvature, 
	Eigen::VectorXd& _gradient_half1_start_max_curvature, Eigen::VectorXd& _gradient_half1_mid_max_curvature, Eigen::VectorXd& _gradient_half1_end_max_curvature, 
	Eigen::VectorXd& _gradient_half2_start_max_curvature, Eigen::VectorXd& _gradient_half2_mid_max_curvature, Eigen::VectorXd& _gradient_half2_end_max_curvature)
{
	SetCrossSection(_cross_section);

	// half1
	half1_triangle_area_ = ComputeControlPolygonArea(-(*_cross_section->GetHalf1Vector01()), *_cross_section->GetHalf1Vector12());
	Eigen::Vector2d half1_midpiont = 0.5 * (*_cross_section->GetHalf1Point2() + *_cross_section->GetPoint0());// the midpoint between point0 and point2
	Eigen::Vector2d half1_point1_to_mid = half1_midpiont - *_cross_section->GetHalf1Point1();
	double half1_length_point1_to_mid = half1_point1_to_mid.norm();
	double pow3_half1_length_point1_to_mid = pow(half1_length_point1_to_mid, 3.0);
	double pow2_half1_triangle_area = half1_triangle_area_ * half1_triangle_area_;
	_half1_mid_max_curvature = pow3_half1_length_point1_to_mid / pow2_half1_triangle_area;
	double pow3_half1_vector01_length = _cross_section->GetHalf1Vector01Length() * _cross_section->GetHalf1Vector01Length() * _cross_section->GetHalf1Vector01Length();
	_half1_start_max_curvature = abs(half1_triangle_area_) / pow3_half1_vector01_length;
	double pow3_half1_vector12_length = _cross_section->GetHalf1Vector12Length() * _cross_section->GetHalf1Vector12Length() * _cross_section->GetHalf1Vector12Length();
	_half1_end_max_curvature = abs(half1_triangle_area_) / pow3_half1_vector12_length;

	//half2
	half2_triangle_area_ = ComputeControlPolygonArea(*_cross_section->GetHalf2Vector12(), -(*_cross_section->GetHalf2Vector01()));
	Eigen::Vector2d half2_midpiont = 0.5 * (*_cross_section->GetHalf2Point2() + *_cross_section->GetPoint0());// the midpoint between point0 and point2
	Eigen::Vector2d half2_point1_to_mid = half2_midpiont - *_cross_section->GetHalf2Point1();
	double half2_length_point1_to_mid = half2_point1_to_mid.norm();
	double pow3_half2_length_point1_to_mid = pow(half2_length_point1_to_mid, 3.0);
	double pow2_half2_triangle_area = half2_triangle_area_ * half2_triangle_area_;
	_half2_mid_max_curvature = pow3_half2_length_point1_to_mid / pow2_half2_triangle_area;
	double pow3_half2_vector01_length = _cross_section->GetHalf2Vector01Length() * _cross_section->GetHalf2Vector01Length() * _cross_section->GetHalf2Vector01Length();
	_half2_start_max_curvature = abs(half2_triangle_area_) / pow3_half2_vector01_length;
	double pow3_half2_vector12_length = _cross_section->GetHalf2Vector12Length() * _cross_section->GetHalf2Vector12Length() * _cross_section->GetHalf2Vector12Length();
	_half2_end_max_curvature = abs(half2_triangle_area_) / pow3_half2_vector12_length;

	// get gradient
	
	// half1 area gradient
	// half1 area = 0.5* |Half1Vector10 x Half1Vector12| = 0.5* |v11 x v21|
	// first get v11 and v21 giradient, then combination
	const Eigen::Vector2d& tangent = *(_cross_section->GetTangent());
	Eigen::Matrix2Xd v11_gradient;//gradient_half1_vector10 w.r.t. 7 parameters; col size =7
	v11_gradient.setZero(Eigen::NoChange, 7);
	// l1
	v11_gradient.col(0) = tangent;
	// theta
	v11_gradient(0, 1) = _cross_section->GetHalf1Vector01Length() * (-tangent(1));
	v11_gradient(1, 1) = _cross_section->GetHalf1Vector01Length() * tangent(0);
	Eigen::Matrix2Xd v21_gradient;//gradient_half1_vector12 w.r.t. 7 parameters; col size =7
	v21_gradient.setZero(Eigen::NoChange, 7);
	// l1 (same with v11)
	v21_gradient.col(0) = tangent;
	// theta (same with v11)
	v21_gradient.col(1) = v11_gradient.col(1);
	// half1 rest (col 3,4)
	v21_gradient.block(0, 3, 2, 2) = Eigen::Matrix2d::Identity();
	// point0 (col 5,6)
	v21_gradient.block(0, 5, 2, 2) = -Eigen::Matrix2d::Identity();
	Eigen::Vector2d half1_vector10 = -(*_cross_section->GetHalf1Vector01());
	Eigen::VectorXd half1_area_gradient;
	AreaGradient(half1_area_gradient, half1_triangle_area_, v11_gradient, v21_gradient, half1_vector10, *_cross_section->GetHalf1Vector12());

	// half2 area gradient
	// half2 area = 0.5* |Half2Vector12 x Half2Vector10| = 0.5* |v12 x v22|
	// first get v12 and v22 giradient, then combination
	Eigen::Matrix2Xd v12_gradient;//gradient_half2_vector12 w.r.t. 7 parameters; col size =7
	v12_gradient.setZero(Eigen::NoChange, 7);
	// theta
	v12_gradient(0, 1) = _cross_section->GetHalf2Vector01Length() * tangent(1);
	v12_gradient(1, 1) = _cross_section->GetHalf2Vector01Length() * (-tangent(0));
	// l2 
	v12_gradient.col(2) = -tangent;
	// point0 (col 5,6)
	v12_gradient.block(0, 5, 2, 2) = -Eigen::Matrix2d::Identity();
	Eigen::Matrix2Xd v22_gradient;//gradient_half2_vector10 w.r.t. 7 parameters; col size =7
	v22_gradient.setZero(Eigen::NoChange, 7);
	// theta (same with v12)
	v22_gradient.col(1) = v12_gradient.col(1);
	// l2 (same with v12)
	v22_gradient.col(2) = v12_gradient.col(2);
	Eigen::Vector2d half2_vector10 = -(*_cross_section->GetHalf2Vector01());
	Eigen::VectorXd half2_area_gradient;
	AreaGradient(half2_area_gradient, half2_triangle_area_, v12_gradient, v22_gradient, *_cross_section->GetHalf2Vector12(), half2_vector10);

	// pow3_length_point1_to_mid_gradient
	Eigen::Matrix2Xd half1_point1_to_mid_gradient;
	half1_point1_to_mid_gradient.setZero(Eigen::NoChange, 7);
	//l1(same with v21)
	half1_point1_to_mid_gradient.col(0) = tangent;
	// theta(same with v21)
	half1_point1_to_mid_gradient.col(1) = v11_gradient.col(1);
	//half1 rest
	half1_point1_to_mid_gradient.block(0, 3, 2, 2) = 0.5 * Eigen::Matrix2d::Identity();
	//point0
	half1_point1_to_mid_gradient.block(0, 5, 2, 2) = -0.5 * Eigen::Matrix2d::Identity();
	Eigen::Vector2d norm_x_pow3_gradient_half1_point1_to_mid;
	NormXPow3Gradient(norm_x_pow3_gradient_half1_point1_to_mid, half1_point1_to_mid, half1_length_point1_to_mid);
	Eigen::VectorXd pow3_half1_length_point1_to_mid_gradient = (norm_x_pow3_gradient_half1_point1_to_mid.transpose() * half1_point1_to_mid_gradient).transpose();

	Eigen::Matrix2Xd half2_point1_to_mid_gradient;
	half2_point1_to_mid_gradient.setZero(Eigen::NoChange, 7);
	// theta(same with v12)
	half2_point1_to_mid_gradient.col(1) = v12_gradient.col(1);
	//l2
	half2_point1_to_mid_gradient.col(2) = v12_gradient.col(2);
	//point0
	half2_point1_to_mid_gradient.block(0, 5, 2, 2) = -0.5 * Eigen::Matrix2d::Identity();
	Eigen::Vector2d norm_x_pow3_gradient_half2_point1_to_mid;
	NormXPow3Gradient(norm_x_pow3_gradient_half2_point1_to_mid, half2_point1_to_mid, half2_length_point1_to_mid);
	Eigen::VectorXd pow3_half2_length_point1_to_mid_gradient = (norm_x_pow3_gradient_half2_point1_to_mid.transpose() * half2_point1_to_mid_gradient).transpose();

	// finally 
	// mid
	MidMaxCurvatureDerivative(_gradient_half1_mid_max_curvature,
		abs(half1_triangle_area_), pow2_half1_triangle_area, half1_area_gradient, pow3_half1_length_point1_to_mid_gradient, pow3_half1_length_point1_to_mid);
	MidMaxCurvatureDerivative(_gradient_half2_mid_max_curvature,
		abs(half2_triangle_area_), pow2_half2_triangle_area, half2_area_gradient, pow3_half2_length_point1_to_mid_gradient, pow3_half2_length_point1_to_mid);
	
	// start
	Eigen::Vector2d norm_x_pow3_gradient_half1_vector10;
	NormXPow3Gradient(norm_x_pow3_gradient_half1_vector10, half1_vector10, _cross_section->GetHalf1Vector01Length());
	Eigen::VectorXd pow3_half1_vector10_gradient = (norm_x_pow3_gradient_half1_vector10.transpose() * v11_gradient).transpose();
	TwoEndsCurvatuveDerivative(_gradient_half1_start_max_curvature, 
		_half1_start_max_curvature, half1_area_gradient, pow3_half1_vector10_gradient, pow3_half1_vector01_length);

	Eigen::Vector2d norm_x_pow3_gradient_half2_vector10;
	NormXPow3Gradient(norm_x_pow3_gradient_half2_vector10, half2_vector10, _cross_section->GetHalf2Vector01Length());
	Eigen::VectorXd pow3_half2_vector10_gradient = (norm_x_pow3_gradient_half2_vector10.transpose() * v22_gradient).transpose();
	TwoEndsCurvatuveDerivative(_gradient_half2_start_max_curvature,
		_half2_start_max_curvature, half2_area_gradient, pow3_half2_vector10_gradient, pow3_half2_vector01_length);

	// end
	Eigen::Vector2d norm_x_pow3_gradient_half1_vector12;
	NormXPow3Gradient(norm_x_pow3_gradient_half1_vector12, *_cross_section->GetHalf1Vector12(), _cross_section->GetHalf1Vector12Length());
	Eigen::VectorXd pow3_half1_vector12_gradient = (norm_x_pow3_gradient_half1_vector12.transpose() * v21_gradient).transpose();
	TwoEndsCurvatuveDerivative(_gradient_half1_end_max_curvature,
		_half1_end_max_curvature, half1_area_gradient, pow3_half1_vector12_gradient, pow3_half1_vector12_length);

	Eigen::Vector2d norm_x_pow3_gradient_half2_vector12;
	NormXPow3Gradient(norm_x_pow3_gradient_half2_vector12, *_cross_section->GetHalf2Vector12(), _cross_section->GetHalf2Vector12Length());
	Eigen::VectorXd pow3_half2_vector12_gradient = (norm_x_pow3_gradient_half2_vector12.transpose() * v12_gradient).transpose();
	TwoEndsCurvatuveDerivative(_gradient_half2_end_max_curvature,
		_half2_end_max_curvature, half2_area_gradient, pow3_half2_vector12_gradient, pow3_half2_vector12_length);
}

void CrossSectionCurvature::SetCrossSection(const CrossSection* _cross_section)
{
	cross_section_ = _cross_section;
}

double CrossSectionCurvature::ComputeControlPolygonArea(const Eigen::Vector2d& _point1, const Eigen::Vector2d& _point2) const
{
	return 0.5 * (_point1[0] * _point2[1] - _point1[1] * _point2[0]);
}

void CrossSectionCurvature::MidMaxCurvatureDerivative(Eigen::VectorXd& _gradient, 
	double _area, double _pow2_area, const Eigen::VectorXd& _derivative_area,
	const Eigen::VectorXd& _derivative_pow3_length_point1_to_mid, double _pow3_length_point1_to_mid) const
{
	_gradient = (_derivative_pow3_length_point1_to_mid * _pow2_area - _pow3_length_point1_to_mid * 2.0 * _area * _derivative_area) /
		(_pow2_area * _pow2_area);
}

void CrossSectionCurvature::TwoEndsCurvatuveDerivative(Eigen::VectorXd& _gradient, 
	double _curvatuve, const Eigen::VectorXd& _derivative_area, const Eigen::VectorXd& _derivative_pow3_vector_length, double _pow3_vector_length) const
{
	_gradient = (_derivative_area - _curvatuve * _derivative_pow3_vector_length)
		/ _pow3_vector_length;
}

void CrossSectionCurvature::AreaGradient(Eigen::VectorXd& _gradient, double _signed_area,
	const Eigen::Matrix2Xd& _v1_graident, const Eigen::Matrix2Xd& _v2_graident,
	const Eigen::Vector2d& _v1, const Eigen::Vector2d& _v2) const
{
	// not include common coefficient yet
	Eigen::Vector2d area_gradient_v1(_v2.y(), -_v2.x());
	Eigen::Vector2d area_gradient_v2(-_v1.y(), _v1.x());

	_gradient = 0.5 * (area_gradient_v1.transpose() * _v1_graident + area_gradient_v2.transpose() * _v2_graident).transpose();
	
	if (_signed_area < 0.0)
	{
		_gradient *= (-1.0);
	}
}

void CrossSectionCurvature::NormXPow3Gradient(Eigen::Vector2d& _gradient, const Eigen::Vector2d& _x, double _norm_x) const
{
	_gradient = 3.0 * _norm_x * _x;
}

double CrossSectionCurvature::Half1StartMaxCurvatureFunction(const Eigen::VectorXd& _x) const
{
	double l1 = _x(0);
	double theta = _x(1);
	double l2 = _x(2);
	std::vector<Eigen::Vector2d> half1_rest(1, Eigen::Vector2d(_x(3), _x(4)));
	Eigen::Vector2d point0(_x(5), _x(6));

	CrossSection current;
	current = *cross_section_;
	current.setControlPoints_SameDegree(half1_rest, l1, theta, l2, std::vector<Eigen::Vector2d>(), point0);
	current.UpdateMoreInfo();

	// half1
	double half1_triangle_area = ComputeControlPolygonArea(-(*current.GetHalf1Vector01()), *current.GetHalf1Vector12());
	Eigen::Vector2d half1_midpiont = 0.5 * (*current.GetHalf1Point2() + *current.GetPoint0());// the midpoint between point0 and point2
	double half1_length_point1_to_mid = (half1_midpiont - *current.GetHalf1Point1()).norm();
	//double half1_mid_max_curvature = (half1_length_point1_to_mid * half1_length_point1_to_mid * half1_length_point1_to_mid) / (half1_triangle_area * half1_triangle_area);
	double half1_start_max_curvature = abs(half1_triangle_area) /
		(current.GetHalf1Vector01Length() * current.GetHalf1Vector01Length() * current.GetHalf1Vector01Length());
	//double half1_end_max_curvature = abs(half1_triangle_area) /

	return half1_start_max_curvature;
}

double CrossSectionCurvature::Half1MidMaxCurvatureFunction(const Eigen::VectorXd& _x) const
{
	double l1 = _x(0);
	double theta = _x(1);
	double l2 = _x(2);
	std::vector<Eigen::Vector2d> half1_rest(1, Eigen::Vector2d(_x(3), _x(4)));
	Eigen::Vector2d point0(_x(5), _x(6));

	CrossSection current;
	current = *cross_section_;
	current.setControlPoints_SameDegree(half1_rest, l1, theta, l2, std::vector<Eigen::Vector2d>(), point0);
	current.UpdateMoreInfo();

	// half1
	double half1_triangle_area = ComputeControlPolygonArea(-(*current.GetHalf1Vector01()), *current.GetHalf1Vector12());
	Eigen::Vector2d half1_midpiont = 0.5 * (*current.GetHalf1Point2() + *current.GetPoint0());// the midpoint between point0 and point2
	double half1_length_point1_to_mid = (half1_midpiont - *current.GetHalf1Point1()).norm();
	double half1_mid_max_curvature = (half1_length_point1_to_mid * half1_length_point1_to_mid * half1_length_point1_to_mid) / (half1_triangle_area * half1_triangle_area);
	//double half1_start_max_curvature = abs(half1_triangle_area) /

	return half1_mid_max_curvature;
}

double CrossSectionCurvature::Half1EndMaxCurvatureFunction(const Eigen::VectorXd& _x) const
{
	double l1 = _x(0);
	double theta = _x(1);
	double l2 = _x(2);
	std::vector<Eigen::Vector2d> half1_rest(1, Eigen::Vector2d(_x(3), _x(4)));
	Eigen::Vector2d point0(_x(5), _x(6));

	CrossSection current;
	current = *cross_section_;
	current.setControlPoints_SameDegree(half1_rest, l1, theta, l2, std::vector<Eigen::Vector2d>(), point0);
	current.UpdateMoreInfo();

	// half1
	double half1_triangle_area = ComputeControlPolygonArea(-(*current.GetHalf1Vector01()), *current.GetHalf1Vector12());
	Eigen::Vector2d half1_midpiont = 0.5 * (*current.GetHalf1Point2() + *current.GetPoint0());// the midpoint between point0 and point2
	double half1_length_point1_to_mid = (half1_midpiont - *current.GetHalf1Point1()).norm();
	//double half1_mid_max_curvature = (half1_length_point1_to_mid * half1_length_point1_to_mid * half1_length_point1_to_mid) / (half1_triangle_area * half1_triangle_area);
	////double half1_start_max_curvature = abs(half1_triangle_area) /
	////	(current.GetHalf1Vector01Length() * current.GetHalf1Vector01Length() * current.GetHalf1Vector01Length());
	double half1_end_max_curvature = abs(half1_triangle_area) /
		(current.GetHalf1Vector12Length() * current.GetHalf1Vector12Length() * current.GetHalf1Vector12Length());

	return half1_end_max_curvature;
}

double CrossSectionCurvature::Half2StartMaxCurvatureFunction(const Eigen::VectorXd& _x) const
{
	double l1 = _x(0);
	double theta = _x(1);
	double l2 = _x(2);
	std::vector<Eigen::Vector2d> half1_rest(1, Eigen::Vector2d(_x(3), _x(4)));
	Eigen::Vector2d point0(_x(5), _x(6));

	CrossSection current;
	current = *cross_section_;
	current.setControlPoints_SameDegree(half1_rest, l1, theta, l2, std::vector<Eigen::Vector2d>(), point0);
	current.UpdateMoreInfo();

	//half2
	double half2_triangle_area = ComputeControlPolygonArea(*current.GetHalf2Vector12(), -(*current.GetHalf2Vector01()));
	Eigen::Vector2d half2_midpiont = 0.5 * (*current.GetHalf2Point2() + *current.GetPoint0());// the midpoint between point0 and point2
	double half2_length_point1_to_mid = (half2_midpiont - *current.GetHalf2Point1()).norm();
	double half2_mid_max_curvature = (half2_length_point1_to_mid * half2_length_point1_to_mid * half2_length_point1_to_mid) / (half2_triangle_area * half2_triangle_area);
	double half2_start_max_curvature = abs(half2_triangle_area) /
		(current.GetHalf2Vector01Length() * current.GetHalf2Vector01Length() * current.GetHalf2Vector01Length());
	double half2_end_max_curvature = abs(half2_triangle_area) /
		(current.GetHalf2Vector12Length() * current.GetHalf2Vector12Length() * current.GetHalf2Vector12Length());

	return half2_start_max_curvature;
}

double CrossSectionCurvature::Half2MidMaxCurvatureFunction(const Eigen::VectorXd& _x) const
{
	double l1 = _x(0);
	double theta = _x(1);
	double l2 = _x(2);
	std::vector<Eigen::Vector2d> half1_rest(1, Eigen::Vector2d(_x(3), _x(4)));
	Eigen::Vector2d point0(_x(5), _x(6));

	CrossSection current;
	current = *cross_section_;
	current.setControlPoints_SameDegree(half1_rest, l1, theta, l2, std::vector<Eigen::Vector2d>(), point0);
	current.UpdateMoreInfo();

	//half2
	double half2_triangle_area = ComputeControlPolygonArea(*current.GetHalf2Vector12(), -(*current.GetHalf2Vector01()));
	Eigen::Vector2d half2_midpiont = 0.5 * (*current.GetHalf2Point2() + *current.GetPoint0());// the midpoint between point0 and point2
	double half2_length_point1_to_mid = (half2_midpiont - *current.GetHalf2Point1()).norm();
	double half2_mid_max_curvature = (half2_length_point1_to_mid * half2_length_point1_to_mid * half2_length_point1_to_mid) / (half2_triangle_area * half2_triangle_area);
	double half2_start_max_curvature = abs(half2_triangle_area) /
		(current.GetHalf2Vector01Length() * current.GetHalf2Vector01Length() * current.GetHalf2Vector01Length());
	double half2_end_max_curvature = abs(half2_triangle_area) /
		(current.GetHalf2Vector12Length() * current.GetHalf2Vector12Length() * current.GetHalf2Vector12Length());

	return half2_mid_max_curvature;
}

double CrossSectionCurvature::Half2EndMaxCurvatureFunction(const Eigen::VectorXd& _x) const
{
	double l1 = _x(0);
	double theta = _x(1);
	double l2 = _x(2);
	std::vector<Eigen::Vector2d> half1_rest(1, Eigen::Vector2d(_x(3), _x(4)));
	Eigen::Vector2d point0(_x(5), _x(6));

	CrossSection current;
	current = *cross_section_;
	current.setControlPoints_SameDegree(half1_rest, l1, theta, l2, std::vector<Eigen::Vector2d>(), point0);
	current.UpdateMoreInfo();
	//half2
	double half2_triangle_area = ComputeControlPolygonArea(*current.GetHalf2Vector12(), -(*current.GetHalf2Vector01()));
	Eigen::Vector2d half2_midpiont = 0.5 * (*current.GetHalf2Point2() + *current.GetPoint0());// the midpoint between point0 and point2
	double half2_length_point1_to_mid = (half2_midpiont - *current.GetHalf2Point1()).norm();
	double half2_mid_max_curvature = (half2_length_point1_to_mid * half2_length_point1_to_mid * half2_length_point1_to_mid) / (half2_triangle_area * half2_triangle_area);
	double half2_start_max_curvature = abs(half2_triangle_area) /
		(current.GetHalf2Vector01Length() * current.GetHalf2Vector01Length() * current.GetHalf2Vector01Length());
	double half2_end_max_curvature = abs(half2_triangle_area) /
		(current.GetHalf2Vector12Length() * current.GetHalf2Vector12Length() * current.GetHalf2Vector12Length());

	return half2_end_max_curvature;
}

void CrossSectionCurvature::compareDifference(const Eigen::VectorXd& _g, const Eigen::VectorXd& _g_fd) const
{
	double average_difference = 0.0;
	double max_difference = -FLT_MAX;
	double min_difference = FLT_MAX;
	for (int iRow = 0; iRow < _g.rows(); ++iRow)
	{

			double difference = abs(_g(iRow) - _g_fd(iRow));
			average_difference += difference;
			if (difference > max_difference)
			{
				max_difference = difference;
			}
			if (difference < min_difference)
			{
				min_difference = difference;
			}
			if (difference > 1e-10)
			{
				cout << "(" << iRow << ") " << "_g1 =" << _g(iRow) << " " << "g2_ =" << _g_fd(iRow) << endl;
			}
		
	}
	average_difference /= double(_g.rows());

	cout << "average_difference=" << average_difference << endl;
	cout << "max_difference=" << max_difference << endl;
	cout << " min_difference=" << min_difference << endl;
}

double CrossSectionCurvature::TriangleAreaHalf1(const Eigen::VectorXd& _x) const
{
	double l1 = _x(0);
	double theta = _x(1);
	double l2 = _x(2);
	std::vector<Eigen::Vector2d> half1_rest(1, Eigen::Vector2d(_x(3), _x(4)));
	Eigen::Vector2d point0(_x(5), _x(6));

	CrossSection current;
	current = *cross_section_;
	current.setControlPoints_SameDegree(half1_rest, l1, theta, l2, std::vector<Eigen::Vector2d>(), point0);
	current.UpdateMoreInfo();
	return ComputeControlPolygonArea(-(*current.GetHalf1Vector01()), *current.GetHalf1Vector12());
}

double CrossSectionCurvature::TriangleAreaHalf2(const Eigen::VectorXd& _x) const
{
	double l1 = _x(0);
	double theta = _x(1);
	double l2 = _x(2);
	std::vector<Eigen::Vector2d> half1_rest(1, Eigen::Vector2d(_x(3), _x(4)));
	Eigen::Vector2d point0(_x(5), _x(6));

	CrossSection current;
	current = *cross_section_;
	current.setControlPoints_SameDegree(half1_rest, l1, theta, l2, std::vector<Eigen::Vector2d>(), point0);
	current.UpdateMoreInfo();
	return ComputeControlPolygonArea(*current.GetHalf2Vector12(), -(*current.GetHalf2Vector01()));
}
