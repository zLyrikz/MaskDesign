#include "Rail.h"
#include "../CurveNSurface/SamplePolylineFromCurve.h"
#include <iostream>
using std::cout;
using std::endl;

Rail::Rail()
{
}

Rail::~Rail()
{
}

void Rail::setRail(const BSpline& _rail)
{
	rail_ = _rail;
}

void Rail::Fitting(const std::vector<Eigen::VectorXd>& _polyline, int _segments, int _degree, int _smoothness)
{
	rail_.LeastSquareFitting_Periodic(_polyline, _segments, _degree, _smoothness);
}

void Rail::Fitting(const std::vector<Eigen::VectorXd>& _polyline, const std::vector<int>& _parts_idx, int _degree, int _smoothness)
{
	rail_.LeastSquareFitting_Periodic(_polyline, _parts_idx, _degree, _smoothness);
}

void Rail::setControlPoints(const double* _x)
{
	int space_dim = rail_.getSpaceDimension();
	int point_dim = rail_.getControlPointDimension();

	std::vector<Eigen::VectorXd> control_points(space_dim, Eigen::VectorXd(point_dim));
	for (int iPoint = 0; iPoint < space_dim; ++iPoint)
	{
		for (int iDim = 0; iDim < point_dim; ++iDim)
		{
			control_points[iPoint][iDim] = _x[iPoint * point_dim + iDim];
		}
	}

	rail_.setControlPoints(std::move(control_points));
}

void Rail::getControlPoints(double* _x) const
{
	int space_dim = rail_.getSpaceDimension();
	int point_dim = rail_.getControlPointDimension();

	const std::vector<Eigen::VectorXd>* control_points = nullptr;
	rail_.getControlPoints(control_points);

	if (control_points != nullptr)
	{
		for (int iPoint = 0; iPoint < space_dim; ++iPoint)
		{
			for (int iDim = 0; iDim < point_dim; ++iDim)
			{
				_x[iPoint * point_dim + iDim] = (*control_points)[iPoint][iDim];
			}
		}
	}
}

int Rail::getDegreeOfFreedom() const
{
	int space_dim = rail_.getSpaceDimension();
	int point_dim = rail_.getControlPointDimension();
	return space_dim * point_dim;
}

void Rail::SaveResults(std::string _file) const
{
	rail_.saveToFile(_file);
	cout << "rail data saved" << endl;
}

void Rail::ReadResults(std::string _file)
{
	rail_.loadFromFile(_file);
}

void Rail::getValue(double _u, Mesh::Point& _value) const
{
	Eigen::VectorXd value;
	rail_.getValue(_u, value);
	_value[0] = value[0];
	_value[1] = value[1];
	_value[2] = value[2];
}

void Rail::getValue(double _u, Eigen::VectorXd& _value) const
{
	rail_.getValue(_u, _value);
}

void Rail::getSecondDerivative(double _u, Eigen::VectorXd& _value)
{
	rail_.getSecondDerivative(_u, _value);
}

void Rail::getFirstDerivative(double _u, Eigen::VectorXd& _value)
{
	rail_.getFirstDerivative(_u, _value);
}

void Rail::getNormalPlane(double _u, Eigen::VectorXd& _normal, Eigen::VectorXd& _value) const
{
	rail_.getDerivative(_u, _normal);
	rail_.getValue(_u, _value);
}

const BSpline& Rail::getRail() const
{
	return rail_;
}

const std::array<double, 2>& Rail::getDomain() const
{
	return rail_.getDomain();
}

int Rail::getSegments() const
{
	return rail_.getSegments();
}

int Rail::getDegree() const
{
	return rail_.getDegree();
}

void Rail::getPolynomialIntervals(std::vector<double>& _knots) const
{
	rail_.getPolynomialIntervals(_knots);
}

void Rail::getRailSamplePoints(std::vector<Mesh::Point>& _points, int _sample_num) const
{
	SamplePolylineFromCurve sample_points;
	sample_points.getPolylineAsPointVector3D<Mesh::Point>(rail_, _sample_num, _points, getDomain());
}

void Rail::getControlPointMesh(Mesh& _control_points) const
{
	const std::vector<Eigen::VectorXd>* constrol_points = nullptr;
	rail_.getControlPoints(constrol_points);

	int space_dim = rail_.getSpaceDimension();
	if (space_dim > 0)
	{
		for (int iPoint = 0; iPoint < space_dim; ++iPoint)
		{
			_control_points.add_vertex(std::move(
				Mesh::Point((*constrol_points)[iPoint][0], (*constrol_points)[iPoint][1], (*constrol_points)[iPoint][2])));
		}
		// add first control point again to form a polygon
		_control_points.add_vertex(std::move(
			Mesh::Point((*constrol_points)[0][0], (*constrol_points)[0][1], (*constrol_points)[0][2])));
	}
}
