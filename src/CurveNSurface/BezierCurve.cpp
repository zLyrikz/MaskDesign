#include "BezierCurve.h"

BezierCurve::BezierCurve()
{
}

BezierCurve::BezierCurve(int _degree)
{
	degree_ = _degree;
	control_points_.resize(degree_ + 1, Eigen::VectorXd::Zero(1));
}

BezierCurve::BezierCurve(const std::vector<std::vector<double>>& _control_points)
{
	control_points1_ = _control_points;
}

void BezierCurve::setControlPoints(const std::vector<Eigen::VectorXd>& _control_points)
{
	control_points_ = _control_points;
	degree_ = control_points_.size() - 1;
}

void BezierCurve::getValue_Degree2(double _x, Eigen::VectorXd& _value) const
{
	_value = 
		control_points_[0] * pow((1 - _x), 2) +
		control_points_[1] * 2 * _x * (1 - _x) +
		control_points_[2] * pow(_x, 2);
}

void BezierCurve::getValue_Degree3(double _x, Eigen::VectorXd& _value) const
{
	_value =
		control_points_[0] * pow((1 - _x), 3) +
		control_points_[1] * 3.0 * _x * pow((1 - _x), 2) +
		control_points_[2] * 3.0 * pow(_x, 2) * (1 - _x) +
		control_points_[3] * pow(_x, 3);
}

void BezierCurve::getValue(double _x, Eigen::VectorXd& _value) const
{
	_value = Eigen::VectorXd::Zero(control_points_[0].rows());
	int degree = control_points_.size() - 1;
	double y = 1 - _x;
	for (int i = 0; i < degree + 1; ++i)
	{
		_value += control_points_[i] * CombinationNchooseK(degree, i) * pow(_x, i) * pow(y, (degree - i));
	}
}

void BezierCurve::getDerivativeCurve(BezierCurve& _derivative)
{
	// check formula of derivative https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/bezier-der.html

	// find control points
	double current_degree = control_points_.size() - 1;
	std::vector<Eigen::VectorXd> control_points(current_degree);// control points of the derivative
	for (int iControl = 0; iControl < current_degree; ++iControl)
	{
		control_points[iControl] = double(current_degree) * (control_points_[iControl + 1] - control_points_[iControl]);
	}

	// construct curve
	_derivative.setControlPoints(control_points);
}

double BezierCurve::getBasisValue(int _id, double _x) const
{
	return CombinationNchooseK(degree_, _id) * pow(_x, _id) * pow(1 - _x, degree_ - _id);
}

double BezierCurve::getBasisDerivative(int id, double _x) const
{
	//https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/bezier-der.html#:~:text=A%20Relation%20Between%20the%20Derivative,curves%20of%20degree%20n%2D1.
	double value = 0.0;
	if (degree_ >= 1)
	{
		BezierCurve lower_order(degree_ - 1);
		if (id == 0)
		{
			value = -degree_ * lower_order.getBasisValue(id, _x);
		}
		else if (id == degree_)
		{
			value = degree_ * lower_order.getBasisValue(id - 1, _x);
		}
		else
		{
			value = degree_ * (lower_order.getBasisValue(id - 1, _x) - lower_order.getBasisValue(id, _x));
		}

	}
	return value;
}

double BezierCurve::CombinationNchooseK(double _n, double _k) const
{
	return std::tgamma(_n + 1) / std::tgamma(_k + 1) / std::tgamma(_n - _k + 1);
}

void BezierCurve::generateCurvePoints(int _point_num)
{
	curve_points_.clear();
	curve_points_.resize(_point_num);

	if (control_points1_.size() == 3) // degree 2
	{
		int dimension = control_points1_[0].size();

		// basis: (1-t)^2,  2(1-t)t,  t^2

		for (int i_point = 0; i_point < _point_num; ++i_point)
		{
			double t = (double)i_point / (double)(_point_num - 1);
			curve_points_[i_point].resize(3);

			for (int i_coordinate = 0; i_coordinate < dimension; ++i_coordinate)
			{
				curve_points_[i_point][i_coordinate] = 
					control_points1_[0][i_coordinate] * pow((1 - t), 2) +
					control_points1_[1][i_coordinate] * 2 * t * (1 - t) +
					control_points1_[2][i_coordinate] * pow(t, 2);
			}
		}
	}
}

void BezierCurve::getCurvePoints(std::vector<std::vector<double>>& _points, int _point_num)
{
	if (_point_num != curve_points_.size() || changing_controlpoints_)
	{
		generateCurvePoints(_point_num);
		changing_controlpoints_ = false;
	}

	_points = curve_points_;
}

void BezierCurve::setControlPoints(const std::vector<std::vector<double>>& _control_points)
{
	control_points1_ = _control_points;
	changing_controlpoints_ = true;
}
