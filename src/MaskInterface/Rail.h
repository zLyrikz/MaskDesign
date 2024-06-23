#pragma once
#include <array>
//#include "../Alglib/interpolation.h"
#include "../MeshViewer/MeshDefinition.h"
#include "../CurveNSurface/BSpline.h"


class Rail
{
public:
	Rail();
	~Rail();

public:
	void setRail(const BSpline& _rail);
	void Fitting(const std::vector<Eigen::VectorXd>& _polyline, int _segments, int _degree, int _smoothness);
	void Fitting(const std::vector<Eigen::VectorXd>& _polyline, const std::vector<int>& _parts_idx, int _degree, int _smoothness);

	// input x size = spline space dim * point dim; (spline space dim = K)
	// control_points[iPoint][iDim] = _x[iPoint * point_dim + iDim]
	void setControlPoints(const double* _x);
	void getControlPoints(double* _x) const;
	int getDegreeOfFreedom() const;

	void SaveResults(std::string _file) const;
	void ReadResults(std::string _file);
public:

	void getValue(double _u, Mesh::Point& _value) const;
	void getValue(double _u, Eigen::VectorXd& _value) const;
	void getSecondDerivative(double _u, Eigen::VectorXd& _value);
	void getFirstDerivative(double _u, Eigen::VectorXd& _value);
	// curve tangent is the normal of the normal plane, the curve value is a point on the normal plane
	// notice the tangent is not normalized
	void getNormalPlane(double _u, Eigen::VectorXd& _normal, Eigen::VectorXd& _value) const;
	const BSpline& getRail() const;
	const std::array<double, 2>& getDomain() const;
	int getSegments() const;
	int getDegree() const;
	// return knots for the spline intervals;
	// first & last knots are domain boundary
	// no duplicated knots, even if it's multiple
	void getPolynomialIntervals(std::vector<double>& _knots) const;

public:
	void getRailSamplePoints(std::vector<Mesh::Point>& _points, int _sample_num = 500) const;

	// for visualization
	
	//// output mesh has the first control point again at the end to form a polygon
	void getControlPointMesh(Mesh& _control_points) const;

private:
	BSpline rail_;



};

