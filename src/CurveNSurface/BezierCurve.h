#pragma once
#include <vector>
#include <Eigen/Core>

// calculate sample points' coordinates on curve
// only support order 3(degree 2) for now
// domain = [0,1]
class BezierCurve
{
public:
	BezierCurve();
	BezierCurve(int _degree);
	BezierCurve(const std::vector<std::vector<double>>& _control_points);

	void setControlPoints(const std::vector<Eigen::VectorXd>& _control_points);

	void getValue_Degree2(double _x, Eigen::VectorXd& _value) const;
	void getValue_Degree3(double _x, Eigen::VectorXd& _value) const;
	void getValue(double _x, Eigen::VectorXd& _value) const;// compute degree from the number of control points
	void getDerivativeCurve(BezierCurve& _derivative);
	double getBasisValue(int id, double _x) const;//id in [0,degree]
	double getBasisDerivative(int id, double _x) const;//id in [0,degree]

	//reference https://stackoverflow.com/questions/67193273/is-there-a-built-in-function-to-calculate-ncr-in-c
	double CombinationNchooseK(double _n, double _k) const;

private:
	int degree_ = 0;
	std::vector<Eigen::VectorXd> control_points_;

/*******************************************************************************************************/
// previous implementation
/*******************************************************************************************************/
public:
	void getCurvePoints(std::vector<std::vector<double>>& _points, int _point_num = 1000);
	void setControlPoints(const std::vector<std::vector<double>>& _control_points);

private:
	void generateCurvePoints(int _point_num);

private:
	std::vector<std::vector<double>> control_points1_; // control_points1_[i] is the i-th control point
	std::vector<std::vector<double>> curve_points_; // curve_points_[i] is the i-th sample point on the curve

	bool changing_controlpoints_ = true;
};

