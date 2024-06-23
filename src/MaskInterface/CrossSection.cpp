#include "CrossSection.h"
#include "CrossSectionHelpr.h"
#include "../CurveNSurface/ChangeFrame.h"
#include "../CurveNSurface/BSpline.h"
#include "../CurveNSurface/BezierCurve.h"
#include <fstream>
#include <iostream>
using std::cout;
using std::endl;

CrossSection::CrossSection()
{
	half1_vector01_ = nullptr;
	half1_point1_ = nullptr;
	half2_vector01_ = nullptr;
	half2_point1_ = nullptr;
	half1_vector12_ = nullptr;
	half2_vector12_ = nullptr;
	half1_vector23_ = nullptr;
	half2_vector23_ = nullptr;
	connect_point_.setZero(2);
}

CrossSection::CrossSection(const CrossSection& _other):
	half1_rest_(_other.half1_rest_), half1_point1_length_(_other.half1_point1_length_), tangent_(_other.tangent_),
	tangent_angle_(_other.tangent_angle_), half2_point1_length_(_other.half2_point1_length_), half2_rest_(_other.half2_rest_),
	connect_point_(_other.connect_point_), half1_degree_(_other.half1_degree_), half2_degree_(_other.half2_degree_)
{
	half1_vector01_ = nullptr;
	half1_point1_ = nullptr;
	half2_vector01_ = nullptr;
	half2_point1_ = nullptr;
	half1_vector12_ = nullptr;
	half2_vector12_ = nullptr;
	half1_vector23_ = nullptr;
	half2_vector23_ = nullptr;
}

CrossSection& CrossSection::operator=(const CrossSection& _other)
{
	half1_rest_ = _other.half1_rest_;
	half1_point1_length_ = _other.half1_point1_length_; 
	tangent_ = _other.tangent_;
	tangent_angle_ = _other.tangent_angle_;
	half2_point1_length_ = _other.half2_point1_length_; 
	half2_rest_ = _other.half2_rest_;
	connect_point_ = _other.connect_point_; 
	half1_degree_ = _other.half1_degree_;
	half2_degree_ = _other.half2_degree_;

	half1_vector01_ = nullptr;
	half1_point1_ = nullptr;
	half2_vector01_ = nullptr;
	half2_point1_ = nullptr;
	half1_vector12_ = nullptr;
	half2_vector12_ = nullptr;
	half1_vector23_ = nullptr;
	half2_vector23_ = nullptr;

	return *this;

}

CrossSection::~CrossSection()
{
	if (half1_point1_ != nullptr)
	{
		delete half1_point1_;
		half1_point1_ = nullptr;
	}
	if (half1_vector01_ != nullptr)
	{
		delete half1_vector01_;
		half1_vector01_ = nullptr;
	}
	if (half2_point1_ != nullptr)
	{
		delete half2_point1_;
		half2_point1_ = nullptr;

	}
	if (half2_vector01_ != nullptr)
	{
		delete half2_vector01_;
		half2_vector01_ = nullptr;

	}
	if (half1_vector12_ != nullptr)
	{
		delete half1_vector12_;
		half1_vector12_ = nullptr;

	}
	if (half2_vector12_ != nullptr)
	{
		delete half2_vector12_;
		half2_vector12_ = nullptr;

	}
	if (half1_vector23_ != nullptr)
	{
		delete half1_vector23_;
		half1_vector23_ = nullptr;

	}
	if (half2_vector23_ != nullptr)
	{
		delete half2_vector23_;
		half2_vector23_ = nullptr;

	}


}

void CrossSection::FitOrigianl3DPolyline(const std::vector<Mesh::Point>& _polyline, const Eigen::Matrix4d& _local_frame, int _degree)
{
	half1_degree_ = _degree;
	half2_degree_ = _degree;

	std::vector<Eigen::VectorXd> polyline_2d;
	compute2DLocalCoordinateOf3DPolyline(_polyline, _local_frame, polyline_2d);

	int half2_begin_idx = SeparatePolyineTo2Halfs(polyline_2d);

	Fit2DPolyline(polyline_2d, half2_begin_idx, _degree);

	TranslateConnectPointToOrigin();

	//for visulaization
	//getBezierCurve(bezier_, bezier_points_);
}

void CrossSection::FitOrigianl3DPolyline(const std::vector<Mesh::Point>& _polyline, const int _half2_begin_idx, const Eigen::Matrix4d& _local_frame, int _degree)
{
	half1_degree_ = _degree;
	half2_degree_ = _degree;

	std::vector<Eigen::VectorXd> polyline_2d;
	compute2DLocalCoordinateOf3DPolyline(_polyline, _local_frame, polyline_2d);

	Fit2DPolyline(polyline_2d, _half2_begin_idx, _degree);

	TranslateConnectPointToOrigin();
}

void CrossSection::FitOrigianl3DPolyline_ChangeableDegree(const std::vector<Mesh::Point>& _polyline, const Eigen::Matrix4d& _local_frame, int _degree_half1, int _degree_half2)
{
	half1_degree_ = _degree_half1;
	half2_degree_ = _degree_half2;

	std::vector<Eigen::VectorXd> polyline_2d;
	compute2DLocalCoordinateOf3DPolyline(_polyline, _local_frame, polyline_2d);

	// get the two separate polylines
	int half2_begin_idx = SeparatePolyineTo2Halfs(polyline_2d);
	Fit2DPolyline_ChangebaleDegree(polyline_2d, half2_begin_idx, Eigen::Vector2d::Zero(), half1_degree_, half2_degree_);
}

void CrossSection::FitOrigianl3DPolyline_ChangeableDegree(const std::vector<Mesh::Point>& _polyline, const int _half2_begin_idx, const Eigen::Matrix4d& _local_frame, int _degree_half1, int _degree_half2)
{
	half1_degree_ = _degree_half1;
	half2_degree_ = _degree_half2;

	std::vector<Eigen::VectorXd> polyline_2d;
	compute2DLocalCoordinateOf3DPolyline(_polyline, _local_frame, polyline_2d);

	// get the two separate polylines
	Fit2DPolyline_ChangebaleDegree(polyline_2d, _half2_begin_idx, polyline_2d[_half2_begin_idx], half1_degree_, half2_degree_);
}

void CrossSection::Fit3DPolyline_FixBoundary(const std::vector<Mesh::Point>& _polyline, const int _half2_begin_idx, const Eigen::Matrix4d& _local_frame, int _degree)
{
	half1_degree_ = _degree;
	half2_degree_ = _degree;

	std::vector<Eigen::VectorXd> polyline_2d;
	compute2DLocalCoordinateOf3DPolyline(_polyline, _local_frame, polyline_2d);

	Fit2DPolyline(polyline_2d, _half2_begin_idx, _degree);

	TranslateConnectPointToOrigin_FixHalf2EndPoint();

	//for visulaization
	//getBezierCurve(bezier_, bezier_points_);
}

void CrossSection::Fit3DPolyline_FixBoundary_ChangeableDegree(const std::vector<Mesh::Point>& _polyline, const int _half2_begin_idx, const Eigen::Matrix4d& _local_frame, int _degree_half1, int _degree_half2)
{
	half1_degree_ = _degree_half1;
	half2_degree_ = _degree_half2;

	std::vector<Eigen::VectorXd> polyline_2d;
	compute2DLocalCoordinateOf3DPolyline(_polyline, _local_frame, polyline_2d);

	// get the two separate polylines
	//int half2_begin_idx = SeparatePolyineTo2Halfs(polyline_2d);
	if (polyline_2d.size() > 0)
	{
		Fit2DPolyline_ChangebaleDegree(polyline_2d, _half2_begin_idx, polyline_2d[_half2_begin_idx], half1_degree_, half2_degree_);
	}
	TranslateConnectPointToOrigin_FixHalf2EndPoint();

}


void CrossSection::setHalf1EndPoint(const Eigen::Vector2d& _end_point)
{
	half1_rest_.back() = _end_point;
}

void CrossSection::setHalf2EndPoint(const Eigen::Vector2d& _end_point)
{
	half2_rest_.back() = _end_point;
}

void CrossSection::setTangentAngle(double _angle)
{
	tangent_angle_ = _angle;
}

void CrossSection::setHalf2Point1Length(double _length)
{
	half2_point1_length_ = _length;
}

void CrossSection::setControlPoints(const std::vector<Eigen::Vector2d>& _half1_rest, double _half1_point1_length, 
	double _tangent_angle, double _half2_point1_length, const std::vector<Eigen::Vector2d>& _half2_rest)
{
	half1_degree_ = _half1_rest.size() + 1;
	half2_degree_ = _half2_rest.size() + 1;

	half1_rest_ = _half1_rest;
	half1_point1_length_ = _half1_point1_length;
	tangent_angle_ = _tangent_angle;
	half2_point1_length_ = _half2_point1_length;
	half2_rest_ = _half2_rest;
}

void CrossSection::setPoint0(const Eigen::Vector2d& _point0)
{
	connect_point_ = _point0;
}

void CrossSection::setControlPoints_SameDegree(const std::vector<Eigen::Vector2d>& _half1_rest, double _half1_point1_length, double _tangent_angle, double _half2_point1_length, const std::vector<Eigen::Vector2d>& _half2_rest_except_end, const Eigen::Vector2d& _connect_point)
{
	setControlPoints_SameDegree(_half1_rest, _half1_point1_length, _tangent_angle, _half2_point1_length, _half2_rest_except_end);
	connect_point_ = _connect_point;
}

void CrossSection::setControlPoints_SameDegree(const std::vector<Eigen::Vector2d>& _half1_rest, double _half1_point1_length, double _tangent_angle, double _half2_point1_length, const std::vector<Eigen::Vector2d>& _half2_rest_except_end)
{
	half1_rest_ = _half1_rest;
	half1_point1_length_ = _half1_point1_length;
	tangent_angle_ = _tangent_angle;
	half2_point1_length_ = _half2_point1_length;
	for (int iRest2 = 0; iRest2 < _half2_rest_except_end.size(); ++iRest2)
	{
		half2_rest_[iRest2] = _half2_rest_except_end[iRest2];
	}
}

void CrossSection::setControlPoints_SameDegree(const std::vector<Eigen::Vector2d>& _half1_rest, double _half1_point1_length, double _tangent_angle,
	double _half2_point1_length, const Eigen::Vector2d& _half2_point2)
{
	half1_rest_ = _half1_rest;
	half1_point1_length_ = _half1_point1_length;
	tangent_angle_ = _tangent_angle;
	half2_point1_length_ = _half2_point1_length;
	half2_rest_[0] = _half2_point2;
}

void CrossSection::getCurveMeshIn3D(Mesh& _curve, const std::vector<Eigen::VectorXd>& _polyline_2d) const
{
	_curve.clear();
	for (int iPoint = 0; iPoint < _polyline_2d.size(); ++iPoint)
	{
		_curve.add_vertex(Mesh::Point(0.0, _polyline_2d[iPoint][0], _polyline_2d[iPoint][1]));
	}
}

void CrossSection::getHalf1EndPoint(Eigen::Vector2d& _end_point) const
{
	_end_point = half1_rest_.back();
}

void CrossSection::getHalf2EndPoint(Eigen::Vector2d& _end_point) const
{
	_end_point = half2_rest_.back();

}

void CrossSection::getControlPoints(std::vector<Eigen::Vector2d>& _half1_rest, double& _half1_point1_length, Eigen::Vector2d& _tangent, double& _half2_point1_length, std::vector<Eigen::Vector2d>& _half2_rest) const
{
	_half1_rest = half1_rest_;
	_half1_point1_length = half1_point1_length_;
	_tangent = tangent_;
	_half2_point1_length = half2_point1_length_;
	_half2_rest = half2_rest_;
}

void CrossSection::getControlPoints(std::vector<Eigen::Vector2d>& _half1_rest, double& _half1_point1_length, double& _tangent_angle, double& _half2_point1_length, std::vector<Eigen::Vector2d>& _half2_rest) const
{
	_half1_rest = half1_rest_;
	_half1_point1_length = half1_point1_length_;
	_tangent_angle = tangent_angle_;
	_half2_point1_length = half2_point1_length_;
	_half2_rest = half2_rest_;
}

void CrossSection::getControlPoints(const std::vector<Eigen::Vector2d>*& _half1_rest, const double*& _half1_point1_length, const double*& _tangent_angle, const double*& _half2_point1_length, const std::vector<Eigen::Vector2d>*& _half2_rest) const
{
	_half1_rest = &half1_rest_;
	_half1_point1_length = &half1_point1_length_;
	_tangent_angle = &tangent_angle_;
	_half2_point1_length = &half2_point1_length_;
	_half2_rest = &half2_rest_;
}

void CrossSection::getControlPoints(const std::vector<Eigen::Vector2d>*& _half1_rest, const double*& _half1_point1_length, const double*& _tangent_angle, const double*& _half2_point1_length, const std::vector<Eigen::Vector2d>*& _half2_rest, const Eigen::Vector2d*& _point0) const
{
	_half1_rest = &half1_rest_;
	_half1_point1_length = &half1_point1_length_;
	_tangent_angle = &tangent_angle_;
	_half2_point1_length = &half2_point1_length_;
	_half2_rest = &half2_rest_;
	_point0 = &connect_point_;
}

void CrossSection::getControlPoints(const std::vector<Eigen::Vector2d>*& _half1_rest, const double*& _half1_point1_length, const double*& _tangent_angle, 
	const double*& _half2_point1_length, const Eigen::Vector2d*& _half2_point2) const
{
	_half1_rest = &half1_rest_;
	_half1_point1_length = &half1_point1_length_;
	_tangent_angle = &tangent_angle_;
	_half2_point1_length = &half2_point1_length_;
	_half2_point2 = &(half2_rest_[0]);
}

double CrossSection::getTangentAngle() const
{
	return tangent_angle_;
}

double CrossSection::getConnectPointX() const
{
	return connect_point_[0];
}

double CrossSection::getHalf1Point2X() const
{
	return half1_rest_[0].x();
}

int CrossSection::getHalf1Degree() const
{
	return half1_degree_;
}

int CrossSection::getHalf2Degree() const
{
	return half2_degree_;
}

void CrossSection::saveToFile(const std::string& _file) const
{
	std::ofstream file(_file);
	file << "point0: " << connect_point_.transpose() << "\n";
	file << "half1_point1_length_: " << half1_point1_length_ << "\n";
	file << "tangent_angle_: " << tangent_angle_ << "\n";
	file << "half2_point1_length_: " << half2_point1_length_ << "\n";
	file << "half1_degree_: " << half1_degree_ << "\n";
	file << "half2_degree_: " << half2_degree_ << "\n";
	file << "half1_rest_size: " << half1_rest_.size() << "\n";
	file << "half1_rest_: ";
	for (const auto& point : half1_rest_) {
		file << point.transpose() << "\n";
	}
	file << "half2_rest_size: " << half2_rest_.size() << "\n";
	file << "half2_rest_: ";
	for (const auto& point : half2_rest_) {
		file << point.transpose() << "\n";
	}
	file.close();
}

void CrossSection::loadFromFile(const std::string& _file)
{
	std::ifstream file(_file);
	std::string dummy;

	file >> dummy;
	if (dummy == "point0:")
	{
		file >> connect_point_[0] >> connect_point_[1];
		file >> dummy >> half1_point1_length_;
	}
	else
	{
		connect_point_.setZero(2);
		file >> half1_point1_length_;
	}

	file >> dummy >> tangent_angle_;
	file >> dummy >> half2_point1_length_;
	file >> dummy >> half1_degree_;
	file >> dummy >> half2_degree_;
	file >> dummy;
	int half1_rest_size;
	file >> half1_rest_size;
	file >> dummy;
	double value1, value2;
	for (int i = 0; i < half1_rest_size; ++i) {
		file >> value1 >> value2;
		Eigen::Vector2d point(value1, value2);
		half1_rest_.push_back(point);
	}
	file >> dummy;
	int half2_rest_size;
	file >> half2_rest_size;
	file >> dummy;
	for (int i = 0; i < half2_rest_size; ++i) {
		file >> value1 >> value2;
		Eigen::Vector2d point(value1, value2);
		half2_rest_.push_back(point);
	}
	file.close();
}

void CrossSection::UpdateMoreInfo()
{
	if (half1_point1_ == nullptr)
	{
		half1_point1_ = new Eigen::Vector2d;
	}
	if (half2_point1_ == nullptr)
	{
		half2_point1_ = new Eigen::Vector2d;
	}
	if (half1_vector01_ == nullptr)
	{
		half1_vector01_ = new Eigen::Vector2d;
	}
	if (half2_vector01_ == nullptr)
	{
		half2_vector01_ = new Eigen::Vector2d;
	}

	// data around the half1 point1
	CrossSectionHelpr convert_tangent;
	convert_tangent.TangentAngle2Vector(tangent_angle_, tangent_);
	*half1_vector01_ = -tangent_ * half1_point1_length_;
	*half1_point1_ = *half1_vector01_ + connect_point_;

	if (half1_rest_.size() > 0)
	{
		half1_end_norm_ = half1_rest_.back().norm();

		if (half1_vector12_ == nullptr)
		{
			half1_vector12_ = new Eigen::Vector2d;
		}
		*half1_vector12_ = half1_rest_[0] - *half1_point1_;

		if (half1_rest_.size() > 1)
		{
			if (half1_vector23_ == nullptr)
			{
				half1_vector23_ = new  Eigen::Vector2d;
			}
			*half1_vector23_ = half1_rest_[1] - half1_rest_[0];
		}
	}


	// data around the half2 point1
	*half2_vector01_ = tangent_ * half2_point1_length_;
	*half2_point1_ = *half2_vector01_ + connect_point_;

	if (half2_rest_.size() > 0)
	{
		if (half2_vector12_ == nullptr)
		{
			half2_vector12_ = new Eigen::Vector2d;
		}
		*half2_vector12_ = half2_rest_[0] - *half2_point1_;

		if (half2_rest_.size() > 1)
		{
			if (half2_vector23_ == nullptr)
			{
				half2_vector23_ = new  Eigen::Vector2d;
			}
			*half2_vector23_ = half2_rest_[1] - half2_rest_[0];

		}
	}
}

bool CrossSection::IsConvex_Half1Point1() const
{
	if (half1_vector01_ == nullptr || half1_vector12_ == nullptr)
	{
		return false;
	}
	else
	{
		Eigen::Vector2d half1_vector10 = -*half1_vector01_;

		//// using the sign of cross product of the two vectors to determine if the angle is larger than 180
		//double cross_product = half1_vector10[0] * (*half1_vector12_)[1] - half1_vector10[1] * (*half1_vector12_)[0];
		//bool is_convex = cross_product > 0;// angle < 180

		return IsAngleSmallerThan180(half1_vector10, (*half1_vector12_));
	}
}

bool CrossSection::IsConvex_Half1Point2() const
{
	if (half1_vector12_ == nullptr || half1_vector23_ == nullptr)
	{
		cout << "[WARNING CrossSection::IsConvex_Half1Point2] required data not updated" << endl;
		return false;
	}
	else
	{
		// get the two vector around the point2
		Eigen::Vector2d half1_vector21 = -(*half1_vector12_);

		return IsAngleSmallerThan180(half1_vector21, *half1_vector23_);
	}
}

bool CrossSection::IsConvex_Half2Point1() const
{
	if (half2_vector01_ == nullptr || half2_vector12_ == nullptr)
	{
		return false;
	}
	else
	{
		Eigen::Vector2d half2_vector10 = -*half2_vector01_;

		//// using the sign of cross product of the two vectors to determine if the angle is larger than 180
		//double cross_product = half1_vector10[0] * (*half1_vector12_)[1] - half1_vector10[1] * (*half1_vector12_)[0];
		//bool is_convex = cross_product > 0;// angle < 180

		return IsAngleSmallerThan180((*half2_vector12_), half2_vector10);
	}
}

bool CrossSection::IsConvex_Half2Point2() const
{
	if (half2_vector12_ == nullptr || half2_vector23_ ==nullptr)
	{
		cout << "[WARNING CrossSection::IsConvex_Half2Point2] required data not updated" << endl;
		return false;
	}
	else
	{
		// get the two vector around the point2
		Eigen::Vector2d half2_vector21 = -(*half2_vector12_);

		return IsAngleSmallerThan180(*half2_vector23_, half2_vector21);
	}
}

double CrossSection::AngleAtHalf1Point1() const
{
	return AngleAtTwoVector(-(*half1_vector01_), *half1_vector12_);
}

double CrossSection::AngleAtHalf1Point2() const
{
	return AngleAtTwoVector(-(*half1_vector12_), *half1_vector23_);
}

double CrossSection::AngleAtHalf2Point1() const
{
	return AngleAtTwoVector(-(*half2_vector01_), *half2_vector12_);
}

double CrossSection::AngleAtHalf2Point2() const
{
	return AngleAtTwoVector(-(*half2_vector12_), *half2_vector23_);
}

double CrossSection::AngleBetweenHalf1Vector01And23() const
{
	return AngleAtTwoVector(*half1_vector01_, *half1_vector23_);
}

double CrossSection::GetHalf1Vector01Length() const
{
	return half1_point1_length_;
}

double CrossSection::GetHalf1Vector12Length() const
{
	return (*half1_vector12_).norm();
}

double CrossSection::GetHalf1Vector23Length() const
{
	return (*half1_vector23_).norm();
}

double CrossSection::GetHalf2Vector01Length() const
{
	return half2_point1_length_;
}

double CrossSection::GetHalf2Vector12Length() const
{
	return (*half2_vector12_).norm();
}

double CrossSection::GetHalf2Vector23Length() const
{
	return (*half2_vector23_).norm();
}

double CrossSection::GetHalf1EndNorm() const
{
	return half1_end_norm_;
}

const Eigen::Vector2d* CrossSection::GetHalf1Point1() const
{
	return half1_point1_;
}

const Eigen::Vector2d* CrossSection::GetHalf1Point2() const
{
	return &half1_rest_[0];
}

const Eigen::Vector2d* CrossSection::GetHalf2Point1() const
{
	return half2_point1_;
}

const Eigen::Vector2d* CrossSection::GetHalf2Point2() const
{
	return &half2_rest_[0];
}

double CrossSection::CrossProductHalf1Vector10AndVector12() const
{
	return CrossProduct2D(-(*half1_vector01_), *half1_vector12_);
}

double CrossSection::CrossProductHalf2Vector12AndVector10() const
{
	return CrossProduct2D(*half2_vector12_, -(*half2_vector01_));
}

const Eigen::Vector2d* CrossSection::GetHalf1Vector01() const
{
	return half1_vector01_;
}

const Eigen::Vector2d* CrossSection::GetHalf1Vector12() const
{
	return half1_vector12_;
}

const Eigen::Vector2d* CrossSection::GetHalf2Vector01() const
{
	return half2_vector01_;
}

const Eigen::Vector2d* CrossSection::GetHalf2Vector12() const
{
	return half2_vector12_;
}

const Eigen::Vector2d* CrossSection::GetPoint0() const
{
	return &connect_point_;
}

const Eigen::Vector2d* CrossSection::GetTangent() const
{
	return &tangent_;
}

void CrossSection::getBezierCurve(Mesh& _bezier_curve, Mesh& _bezier_points) const
{
	Eigen::Vector2d tangent;
	CrossSectionHelpr convert_tangent;
	convert_tangent.TangentAngle2Vector(tangent_angle_, tangent);

	// half1
	std::vector<Eigen::VectorXd> bezier_points1(half1_rest_.size() + 2);
	for (int iRest = 0; iRest < half1_rest_.size(); ++iRest)
	{
		bezier_points1[iRest] = half1_rest_[half1_rest_.size() - 1 - iRest];
	}
	bezier_points1[bezier_points1.size() - 2] = -tangent * half1_point1_length_ + connect_point_;
	*(bezier_points1.end() - 1) = connect_point_;

	//half2

	std::vector<Eigen::VectorXd> bezier_points2(half2_rest_.size() + 2);
	bezier_points2[0] = connect_point_;
	bezier_points2[1] = tangent * half2_point1_length_ + connect_point_;
	for (int iRest = 0; iRest < half2_rest_.size(); ++iRest)
	{
		bezier_points2[iRest + 2] = half2_rest_[iRest];
	}

	// get beizer point as mesh
	_bezier_points.clear();
	for (int iControl = 0; iControl < bezier_points1.size(); ++iControl)
	{
		_bezier_points.add_vertex(Mesh::Point(0.0, bezier_points1[iControl](0), bezier_points1[iControl](1)));
	}
	for (int iControl = 1; iControl < bezier_points2.size(); ++iControl)
	{
		_bezier_points.add_vertex(Mesh::Point(0.0, bezier_points2[iControl](0), bezier_points2[iControl](1)));
	}

	// get bezier curve
	BezierCurve bezier1;
	bezier1.setControlPoints(std::move(bezier_points1));
	BezierCurve bezier2;
	bezier2.setControlPoints(std::move(bezier_points2));
	std::vector<Eigen::VectorXd> beizer_sample;

	for (int i = 0; i < 100; ++i)
	{
		Eigen::VectorXd value;
		bezier1.getValue(double(i) / 99, value);
		beizer_sample.push_back(Eigen::Vector2d(value));
	}
	for (int i = 0; i < 100; ++i)
	{
		Eigen::VectorXd value;
		bezier2.getValue(double(i) / 99, value);
		beizer_sample.push_back(Eigen::Vector2d(value));
	}
	getCurveMeshIn3D(_bezier_curve, beizer_sample);
}

void CrossSection::FlattenHalf1()
{
	// after flatten, rest control points are represented as it's length to origin, direction is the tangent
	// find the max length, then set all the rest points to this length
	double max_length = 0.0;// no less than 0
	for (int iRest = 0; iRest < half1_rest_.size(); ++iRest)
	{
		// project to tangent line
		double length = -half1_rest_[iRest].dot(tangent_);//(tangent point to half2, so dot product should (tend to) be negative)
		assert(abs(tangent_.norm() - 1.0) < 1e-5);
		if (max_length < length)
		{
			max_length = length;
		}
	}

	// set all the rest points to max length
	Eigen::Vector2d max_length_point = max_length * ( - tangent_);
	for (int iRest = 0; iRest < half1_rest_.size(); ++iRest)
	{
		half1_rest_[iRest] = max_length_point;
	}
}

void CrossSection::FlattenEndPoint_Half1()
{
	// project to the line of two mid points
	Eigen::Vector2d mid_point1 = -tangent_ * half1_point1_length_ + connect_point_;
	Eigen::Vector2d line_direction = half1_rest_[0] - mid_point1;
	Eigen::Vector2d line_normal(-line_direction[1], line_direction[0]);
	line_normal.normalize();
	Eigen::Vector2d projected = half1_rest_.back() + (mid_point1 - half1_rest_.back()).dot(line_normal) * line_normal;

	//check if projected point is lied outside of the two mid points
	// if not, set it as the second mid point
	double lambda = projected[0] / line_direction[0];
	if (abs(line_direction[0]) < 1e-5)
	{
		lambda = projected[1] / line_direction[1];
	}
	if (lambda < 1 - 1e-5)
	{
		projected = half1_rest_[0];
	}
	
	half1_rest_.back() = std::move(projected);
}

void CrossSection::FlattenHalf1_ParallelToAxisX()
{

	// after flatten, rest control points are represented as it's length to first mid point
	// find the max length, then set all the rest points to this length

	Eigen::Vector2d mid_point1 = -tangent_ * half1_point1_length_ + connect_point_;
	Eigen::Vector2d line_direction(-1.0, 0.0);
	Eigen::Vector2d line_normal(0.0, 1.0);
	double max_length = 0.0;// no less than 0
	for (int iRest = 0; iRest < half1_rest_.size(); ++iRest)
	{
		// project to the two mid points on the line
		Eigen::Vector2d projected = half1_rest_[iRest] + (mid_point1 - half1_rest_[iRest]).dot(line_normal) * line_normal;

		double length = projected[0] / line_direction[0];
		if (max_length < length)
		{
			max_length = length;
		}
	}

	// set all the rest points to max length
	Eigen::Vector2d max_length_point = max_length * line_direction + mid_point1;
	for (int iRest = 0; iRest < half1_rest_.size(); ++iRest)
	{
		half1_rest_[iRest] = max_length_point;
	}
}

void CrossSection::compute2DLocalCoordinateOf3DPolyline(const std::vector<Mesh::Point>& _3d_polyline, const Eigen::Matrix4d& _local_frame, std::vector<Eigen::VectorXd>& _2d_polyline) const
{
	// get transformation matrix
	Eigen::Matrix4d world_frame(Eigen::Matrix4d::Identity());
	Eigen::Matrix4d transform_to_local;
	ChangeFrame change_frame(world_frame, _local_frame, transform_to_local);

	// transform
	int num_point = _3d_polyline.size();
	_2d_polyline.resize(num_point);
	for (int iPoint = 0; iPoint < num_point; ++iPoint)
	{
		Eigen::Vector4d point_world(_3d_polyline[iPoint][0], _3d_polyline[iPoint][1], _3d_polyline[iPoint][2], 1.0);
		Eigen::Vector4d point_local = transform_to_local * point_world;
		_2d_polyline[iPoint] = Eigen::Vector2d(point_local[1], point_local[2]);
	}
}

int CrossSection::SeparatePolyineTo2Halfs(std::vector<Eigen::VectorXd>& _2d_polyline) const
{
	if (_2d_polyline.size() > 0)
	{	
		// make half1 at the front, half2 at the back
		if (_2d_polyline[0][0] > 0.0)
		{
			std::reverse(_2d_polyline.begin(), _2d_polyline.end());
		}

		int half2_begin_idx = 0;
		while (_2d_polyline[half2_begin_idx][0] <= 0.0)
		{
			++half2_begin_idx;
			if (half2_begin_idx == _2d_polyline.size())
			{
				cout << "CrossSection::SeparatePolyineTo2Halfs unexpected case, polyline not properly in two halfs" << endl;
				break;
			}
		}
		return half2_begin_idx;
	}
	
	return -1;
}

void CrossSection::ComputeBezierPointsFromBSpline(const BSpline& _b_spline)
{
	// get Bezier points (then translate them, i.e. a translation that move the one connect two parts to origin)
	std::vector<Eigen::VectorXd> bezier_points;// size = 2(degree+1)
	_b_spline.ConvertToBeizerCurves(bezier_points);
	int degree = _b_spline.getDegree();

	connect_point_ = bezier_points[degree];
	half1_rest_.clear();
	half2_rest_.clear();
	int num_rest_control_point = degree - 1;
	half1_rest_.reserve(num_rest_control_point);
	half2_rest_.reserve(num_rest_control_point);
	for (int iRest = 0; iRest < num_rest_control_point; ++iRest)
	{
		half1_rest_.push_back(Eigen::Vector2d(bezier_points[(degree - 2) - iRest]));
		half2_rest_.push_back(Eigen::Vector2d(bezier_points[(degree + 3) + iRest]));
	}
	Eigen::Vector2d half1_point1 = bezier_points[degree - 1];
	Eigen::Vector2d half2_point1 = bezier_points[degree + 2];
	tangent_ = half2_point1 - connect_point_;
	half2_point1_length_ = tangent_.norm();
	half1_point1_length_ = (half1_point1 - connect_point_).norm();
	tangent_.normalize();
	CrossSectionHelpr convert_tangent;
	tangent_angle_ = convert_tangent.TangentVector2Angle(tangent_);
}

void CrossSection::Fit2DPolyline(const std::vector<Eigen::VectorXd>& _polyline_2d, const int _half2_begin_idx, int _degree)
{
	//cout << "_half2_begin_idx=" << _half2_begin_idx << endl;
	//cout << "_polyline_2d size()=" << _polyline_2d.size() << endl;
	int smoothness = 1;
	BSpline spline;
	spline.LeastSquareFitting(_polyline_2d, std::vector<int>{_half2_begin_idx}, _degree, smoothness);

	ComputeBezierPointsFromBSpline(spline);

	/*******************************************************************************************************************/
	//// for visualization
	// render bspline and original polyline in 2d
	if (1)
	{
		//cout << "_polyline_2d.size() =" << _polyline_2d.size() << endl;
		getCurveMeshIn3D(polyline_2d_, _polyline_2d);


		std::vector<Eigen::VectorXd> spline_2d;
		for (int i = 0; i < 100; ++i)
		{
			Eigen::VectorXd value;
			spline.getValue(double(i) / 99, value);
			spline_2d.push_back(value);
		}
		getCurveMeshIn3D(spline_, spline_2d);
	}
	/*******************************************************************************************************************/

}

void CrossSection::Fit2DPolyline_ChangebaleDegree(const std::vector<Eigen::VectorXd>& _polyline_2d, const int _half2_begin_idx,
	const Eigen::Vector2d& _connect_point, int _degree_half1, int _degree_half2)
{
	std::vector<Eigen::VectorXd> polyline_2d_half1(_polyline_2d.begin(), _polyline_2d.begin() + _half2_begin_idx);
	std::vector<Eigen::VectorXd> polyline_2d_half2(_polyline_2d.begin() + _half2_begin_idx, _polyline_2d.end());

	// fit each half
	BSpline half1_polynomial;
	// last control point set as origin
	half1_polynomial.PolynomialLeastSquareFitting_WithSomeControlPointsKnown(
		std::vector<int>(1, _degree_half1), std::vector<Eigen::VectorXd>(1, _connect_point), polyline_2d_half1, _degree_half1);
	BSpline half2_polynomial;
	// first control point set as origin
	half2_polynomial.PolynomialLeastSquareFitting_WithSomeControlPointsKnown(
		std::vector<int>(1, 0), std::vector<Eigen::VectorXd>(1, _connect_point), polyline_2d_half2, _degree_half2);

	// notice since the spline degenerates to polynomial, the b spline basis becomes bezier basis

	//1.set the half1 control points as the fitted 
	// note for half1, the point order is reversed
	const std::vector<Eigen::VectorXd>* half1_control_points = nullptr;
	half1_polynomial.getControlPoints(half1_control_points);
	int num_control_half1 = (*half1_control_points).size();
	assert(num_control_half1 == _degree_half1 + 1);
	if (num_control_half1 > 1)
	{
		tangent_[0] = _connect_point(0) - (*half1_control_points)[num_control_half1 - 2](0);
		tangent_[1] = _connect_point(1) - (*half1_control_points)[num_control_half1 - 2](1);
		half1_point1_length_ = tangent_.norm();
		tangent_.normalize();
		CrossSectionHelpr converter;
		tangent_angle_ = converter.TangentVector2Angle(tangent_);

		int num_rest = num_control_half1 - 2;
		half1_rest_.resize(num_rest);
		for (int iRest = 0; iRest < num_rest; ++iRest)
		{
			half1_rest_[iRest] = (*half1_control_points)[num_control_half1 - 3 - iRest];
		}
	}
	//2.set the half2 control points as the fitted; yet for the second control point, project it to the tangent line
	const std::vector<Eigen::VectorXd>* half2_control_points = nullptr;
	half2_polynomial.getControlPoints(half2_control_points);
	int num_control_half2 = (*half2_control_points).size();
	assert(num_control_half2 == _degree_half2 + 1);
	if (num_control_half2 > 1)
	{
		half2_point1_length_ = ((*half2_control_points)[1] - _connect_point).dot(tangent_);
		if (half2_point1_length_ < 0)
		{
			cout << "[WARNING CrossSection::Fit2DPolyline_ChangebaleDegree] half2_point1_length_ < 0, half2_point1_length_=" << half2_point1_length_ << endl;
			//half2_point1_length_ = 0.1;
		}

		int num_rest = num_control_half2 - 2;
		half2_rest_.resize(num_rest);
		for (int iRest = 0; iRest < num_rest; ++iRest)
		{
			half2_rest_[iRest] = (*half2_control_points)[iRest + 2];
		}
	}
	connect_point_ = _connect_point;
}

void CrossSection::TranslateConnectPointToOrigin()
{
	assert(half1_rest_.size() == half2_rest_.size());
	for (int iRest = 0; iRest < half1_rest_.size(); ++iRest)
	{
		half1_rest_[iRest] -= connect_point_;
		half2_rest_[iRest] -= connect_point_;
	}
	connect_point_ = Eigen::Vector2d::Zero();
}

void CrossSection::TranslateConnectPointToOrigin_FixHalf2EndPoint()
{
	for (int iRest = 0; iRest < half1_rest_.size(); ++iRest)
	{
		half1_rest_[iRest] -= connect_point_;
	}
	for (int iRest = 0; iRest < half2_rest_.size() - 1; ++iRest)// notice not for the last one(end point)
	{
		half2_rest_[iRest] -= connect_point_;
	}
	connect_point_ = Eigen::Vector2d::Zero();
}

bool CrossSection::IsAngleSmallerThan180(const Eigen::Vector2d& _vector1, const Eigen::Vector2d& _vector2) const
{
	// using the sign of cross product of the two vectors to determine if the angle is larger than 180
	// angle < 180 iff cross_product > 0
	//double cross_product = _vector1[0] * _vector2[1] - _vector1[1] * _vector2[0];
	double cross_product = CrossProduct2D(_vector1, _vector2);
	return cross_product > 0;
}

double CrossSection::CrossProduct2D(const Eigen::Vector2d& _vector1, const Eigen::Vector2d& _vector2) const
{
	return _vector1[0] * _vector2[1] - _vector1[1] * _vector2[0];
}

double CrossSection::AngleAtTwoVector(const Eigen::Vector2d& _vector1, const Eigen::Vector2d& _vector2) const
{
	return acos((_vector1.normalized()).dot(_vector2.normalized()));
}





