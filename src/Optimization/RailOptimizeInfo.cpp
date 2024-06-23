#include "RailOptimizeInfo.h"
#include "../Evaluation/Variance.h"
#include "../Evaluation/GaussianQuadrature.h"
#include <iostream>
using std::cout;
using std::endl;
RailOptimizeInfo::RailOptimizeInfo()
{
	num_control_point_ = 0;
	control_point_dim_ = 0;
	dim_x_ = 0;
	initial_max_y_ = 0.0;
	initial_min_y_ = 0.0;

	initial_x_ = nullptr;

	box_edge_radius_ = 0;

	distance_to_connector_ = 0.0;
	connector_ = nullptr;
	sample_num_ = 500;

	SetWeights();
}

RailOptimizeInfo::~RailOptimizeInfo()
{
	if (initial_x_ != nullptr)
	{
		delete[] initial_x_;
		initial_x_ = nullptr;
	}
}

void RailOptimizeInfo::SetHeadSdf(const std::vector<const tmd::TriangleMeshDistance*>& _head_sdfs)
{
	head_sdfs_ = _head_sdfs;
}

void RailOptimizeInfo::SetInitialRail(const Rail* _initial_rail)
{
	rail_ = *_initial_rail;
	SetRailControlPointInfo(rail_.getRail().getSpaceDimension(), rail_.getRail().getControlPointDimension());
	// find the top and bottom y value of the rail curve
	int sample_num = 1000;
	std::vector<Mesh::Point> rail_points;
	rail_.getRailSamplePoints(rail_points, sample_num);
	int top_id = 0;
	int bottom_id = 0;
	FindRailTopBottomPoint(rail_points, top_id, bottom_id);
	initial_max_y_ = rail_points[top_id][1];
	initial_min_y_ = rail_points[bottom_id][1];

	if (initial_x_ != nullptr)
	{
		delete[] initial_x_;
		initial_x_ = nullptr;
	}
	initial_x_ = new double[dim_x_];
	GetXFromRail(rail_, initial_x_);
}

void RailOptimizeInfo::SetBoxConstraintRadius(double _box_edge_radius)
{
	box_edge_radius_ = _box_edge_radius;

}

void RailOptimizeInfo::SetWeights(double _objective, double _cushion_width, double _strain_energy, double _symmetry_energy, double _size_regularity, double _angle_energy)
{
	objective_ = _objective;
	cushion_width_ = _cushion_width;
	strain_energy_ = _strain_energy;
	symmetry_energy_ = _symmetry_energy;
	size_regularity_ = _size_regularity;
	angle_energy_ = _angle_energy;
}

void RailOptimizeInfo::SetConnectorDistanceTree(const tmd::TriangleMeshDistance* _connector)
{
	connector_ = _connector;
}

void RailOptimizeInfo::FindBmcRailToConnectorDistance(const Rail& _bmc_rail)
{
	std::vector<Mesh::Point> rail_points;
	_bmc_rail.getRailSamplePoints(rail_points, sample_num_);
	distance_to_connector_ = DistanceToConnector(rail_points);
}

void RailOptimizeInfo::GetRailFromX(Rail& _rail, const double* _x) const
{
	_rail.setControlPoints(_x);
}

void RailOptimizeInfo::GetXFromRail(const Rail& _rail, double* _x) const
{
	_rail.getControlPoints(_x);
}

int RailOptimizeInfo::GetDimX() const
{
	return dim_x_;
}

const double* RailOptimizeInfo::GetInitialX() const
{
	return initial_x_;
}

double RailOptimizeInfo::ObjectiveFunction(const double* _x)
{
	// get rail curve
	GetRailFromX(rail_, _x);

	// sample points on rail, sample num initialized to 500
	std::vector<Mesh::Point> rail_points;
	rail_.getRailSamplePoints(rail_points, sample_num_);

	// ==========================================objective==========================================
	double objective = 0.0;
	double sum_distance = 0.0;

	// average face expressions
	for (int iExpression = 0; iExpression < head_sdfs_.size(); ++iExpression)
	{
		std::vector<double> distances_to_head(sample_num_);
		for (int iPoint = 0; iPoint < sample_num_; ++iPoint)
		{
			tmd::Result result = head_sdfs_[iExpression]->signed_distance({ rail_points[iPoint][0], rail_points[iPoint][1], rail_points[iPoint][2] });
			distances_to_head[iPoint] = result.distance;
			//sum_distance += abs(distances_to_head[iPoint]);// try not use this objective now
		}
		Variance compute_variance;
		objective += compute_variance.computeVariance(distances_to_head);
	}
	objective /= double(head_sdfs_.size());

	double strain_energy = StrainEnergy(rail_);

	int top_id = 0;
	int bottom_id = 0;
	FindRailTopBottomPoint(rail_points, top_id, bottom_id);
	double symmetry_energy = SymmetryEnergy(rail_points, top_id, bottom_id, 600);

	double top_y = rail_points[top_id][1];
	double bottom_y = rail_points[bottom_id][1];
	double size_regularity = SizeRegularity(top_y, bottom_y);
	double angle_energy = WideTopAngleEnergy(_x);
	double cushion_width = CushionWidth(rail_points);

	double value = objective_ * objective + cushion_width_ * cushion_width + strain_energy_ * strain_energy + symmetry_energy_ * symmetry_energy + size_regularity_ * size_regularity + angle_energy_ * angle_energy;
	//double value = 1e3 * objective + 2e-4 * cushion_width + 5.0 * 1e-6 * strain_energy + 0.01 * symmetry_energy + /*1e8*/ 1e8 * size_regularity + 120.0 * angle_energy;
	//double value = 1e3 * objective + 2e-4 * cushion_width + 10.0 * 1e-6 * strain_energy + 0.01 * symmetry_energy + /*1e8*/ 1e8 * size_regularity + 80.0 * angle_energy;
	return value;
}

bool RailOptimizeInfo::IsFeasible(const double* _x) const
{
	return true;

}

double RailOptimizeInfo::StrainEnergy(Rail& _rail) const
{
	//since on D^2 r is just a spline, we can use Guassian Quadrature to accurately compute the strain energy on each polynomial interval

	// 1 get rail intervals
	std::vector<double> knots;
	_rail.getPolynomialIntervals(knots);
	int interval_number = knots.size() - 1;

	// 2 Gaussian Quadrature of integrand (D^2 r)^2 on each interval (degree 2 polynomial)
	std::function<double(double)> rail_second_derivative_inner_product = std::bind(&RailOptimizeInfo::StrainEnergyIntegrand, this, std::placeholders::_1, _rail);
	double strain_energy = 0.0;
	for (int iInterval = 0; iInterval < interval_number; ++iInterval)
	{
		GaussianQuadrature numerical_intergrate;
		strain_energy += numerical_intergrate.TwoPointQuadrature(knots[iInterval], knots[iInterval + 1], rail_second_derivative_inner_product);
	}

	return strain_energy;
}

double RailOptimizeInfo::StrainEnergyIntegrand(double _x, Rail& _rail) const
{
	Eigen::VectorXd rail_second_derivative;
	_rail.getSecondDerivative(_x, rail_second_derivative);

	return rail_second_derivative.dot(rail_second_derivative);
}

void RailOptimizeInfo::FindRailTopBottomPoint(const std::vector<Mesh::Point>& _sample_rail_points, int& _top_id, int& _bottom_id) const
{
	int num_point = _sample_rail_points.size();
	_top_id = 0;   // max y value
	_bottom_id = 0;// min y value
	for (int iPoint = 1; iPoint < num_point; ++iPoint)
	{
		if (_sample_rail_points[iPoint][1] < _sample_rail_points[_bottom_id][1])
		{
			_bottom_id = iPoint;
		}
		if (_sample_rail_points[iPoint][1] > _sample_rail_points[_top_id][1])
		{
			_top_id = iPoint;
		}
	}
	assert(_top_id != _bottom_id);
}

double RailOptimizeInfo::SymmetryEnergy(const std::vector<Mesh::Point>& _sample_rail_points, int _top_id, int _bottom_id, int _sample_y) const
{
	double symmetry_energy = 0.0;
	double width_energy = 0.0;
	if (_sample_rail_points.size() > 1)
	{
		int num_point = _sample_rail_points.size();

		//1 separate the points to two lists: one from top to bottom, the other from bottom to top
		std::vector<Mesh::Point> top2bottom;
		std::vector<Mesh::Point> bottom2top;
		if (_top_id < _bottom_id)
		{
			std::copy(_sample_rail_points.begin() + _top_id, _sample_rail_points.begin() + _bottom_id, std::back_inserter(top2bottom));
			int size_bottom2top = num_point - top2bottom.size() + 2;//+2 is to include the top and bottom points
			bottom2top.reserve(size_bottom2top);
			for (int i = 0; i < size_bottom2top; ++i)
			{
				bottom2top.push_back(_sample_rail_points[(_bottom_id + i) % num_point]);
			}

		}
		else // top_id > bottom_id
		{
			std::copy(_sample_rail_points.begin() + _bottom_id, _sample_rail_points.begin() + _top_id, std::back_inserter(bottom2top));
			int size_top2bottom = num_point - bottom2top.size() + 2;//+2 is to include the top and bottom points
			top2bottom.reserve(size_top2bottom);
			for (int i = 0; i < size_top2bottom; ++i)
			{
				top2bottom.push_back(_sample_rail_points[(_top_id + i) % num_point]);
			}
		}

		//2 scan the two lists with many different y values
		// at each y, acquire the intersect points(a linear combination between two points) with the two lists
		double top_y = _sample_rail_points[_top_id][1];
		double bottom_y = _sample_rail_points[_bottom_id][1];
		double interval = (top_y - bottom_y) / double(_sample_y + 1);
		int scan_top2bottom_id = 1;
		int scan_bottom2top_id = 1;
		std::vector<double> top2bottom_intersects(_sample_y);// x values
		std::vector<double> bottom2top_intersects(_sample_y);// x values
		// scan from bottom to top
		for (int iY = 0; iY < _sample_y; ++iY)
		{			
			double y = bottom_y + double(iY + 1) * interval;

			// check top2bottom list from back to front
			while ((scan_top2bottom_id < top2bottom.size() - 1) && top2bottom[top2bottom.size() - 1 - scan_top2bottom_id][1] < y)
			{
				++scan_top2bottom_id;
			}
			// note this doesn't accurately work for non-monotone y values, but program should not crash
			int top2bottom_last = top2bottom.size() - 1 - scan_top2bottom_id + 1;
			int top2bottom_next = top2bottom.size() - 1 - scan_top2bottom_id;
			double top2bottom_last_y = top2bottom[top2bottom_last][1];
			double top2bottom_next_y = top2bottom[top2bottom_next][1];
			assert(y > top2bottom_last_y && y < top2bottom_next_y);
			double top2bottom_gap = top2bottom_next_y - top2bottom_last_y;
			top2bottom_intersects[iY] = (y - top2bottom_last_y) / top2bottom_gap * top2bottom[top2bottom_next][0]
				+ (top2bottom_next_y - y) / top2bottom_gap * top2bottom[top2bottom_last][0];

			// check bottom2top list from front to back
			while ((scan_bottom2top_id < bottom2top.size() - 1) && bottom2top[scan_bottom2top_id][1] < y)
			{
				++scan_bottom2top_id;
			}
			// note this doesn't accurately work for non-monotone y values, but program should not crash
			int bottom2top_last = scan_bottom2top_id - 1;
			int bottom2top_next = scan_bottom2top_id;
			double bottom2top_last_y = bottom2top[bottom2top_last][1];
			double bottom2top_next_y = bottom2top[bottom2top_next][1];
			assert(y > bottom2top_last_y && y < bottom2top_next_y);
			double bottom2top_gap = bottom2top_next_y - bottom2top_last_y;
			bottom2top_intersects[iY] = (y - bottom2top_last_y) / bottom2top_gap * bottom2top[bottom2top_next][0]
				+ (bottom2top_next_y - y) / bottom2top_gap * bottom2top[bottom2top_last][0];
		}

		//3 compute symmetry
		for (int iY = 0; iY < _sample_y; ++iY)
		{
			// reflect one list w.r.t. YOZ plane
			double reflect_x = -top2bottom_intersects[iY];
			double distance = reflect_x - bottom2top_intersects[iY];
			symmetry_energy += distance * distance;
		}
	}

	return symmetry_energy;
}

double RailOptimizeInfo::SizeRegularity(double _max_y, double _min_y) const
{
	double difference_top_y = _max_y - initial_max_y_;
	return difference_top_y * difference_top_y;
}

double RailOptimizeInfo::WideTopAngleEnergy(const double* _x) const
{
	// only consider XY plane
	Eigen::Vector2d top_control_point(_x[2 * 3], _x[2 * 3 + 1]);
	Eigen::Vector2d left_of_top_control_point(_x[1 * 3], _x[1 * 3 + 1]);
	Eigen::Vector2d right_of_top_control_point(_x[3 * 3], _x[3 * 3 + 1]);
	Eigen::Vector2d vector1 = left_of_top_control_point - top_control_point;
	Eigen::Vector2d vector2 = right_of_top_control_point - top_control_point;
	vector1.normalize();
	vector2.normalize();
	double angle = acos(vector1.dot(vector2));
	return -angle;
}

double RailOptimizeInfo::WideTopEnergy(const double* _x) const
{
	double top_control_point_y = _x[2 * 3 + 1];
	double left_of_top_control_point_y = _x[1 * 3 + 1];
	double right_of_top_control_point_y = _x[3 * 3 + 1];

	double left_y_distance_to_top = top_control_point_y - left_of_top_control_point_y;
	double right_y_distance_to_top = top_control_point_y - right_of_top_control_point_y;

	return left_y_distance_to_top * left_y_distance_to_top + right_y_distance_to_top * right_y_distance_to_top - 276;
}

double RailOptimizeInfo::WideUpperCurve(const double* _x) const
{
	double distance_control13 = _x[1 * 3] - _x[3 * 3];
	double distance_control04 = _x[0 * 3] - _x[4 * 3];
	double distance_control_last5 = _x[(dim_x_/3 - 1) * 3] - _x[5 * 3];//last control point id = (dim_x_/3 - 1)
	return -distance_control13 * distance_control13 - distance_control04 * distance_control04 - distance_control_last5 * distance_control_last5;
}

double RailOptimizeInfo::DistanceToConnector(const std::vector<Mesh::Point>& _sample_rail_points) const
{
	double sum_distance = 0.0;
	for (int iPoint = 0; iPoint < sample_num_; ++iPoint)
	{
		tmd::Result result = connector_->signed_distance({ _sample_rail_points[iPoint][0], _sample_rail_points[iPoint][1], _sample_rail_points[iPoint][2] });
		sum_distance += result.distance;
	}
	return sum_distance;
}

double RailOptimizeInfo::CushionWidth(const std::vector<Mesh::Point>& _sample_rail_points) const
{
	return pow((DistanceToConnector(_sample_rail_points) - distance_to_connector_), 2.0);
}

bool RailOptimizeInfo::IsFeasible_BoxConstraint(const double* _initial_x, const double* _current_x, const double& _box_edge_radius) const
{
	for (int iX = 0; iX < dim_x_; ++iX)
	{
		double difference = _current_x[iX] - _initial_x[iX];
		if ((difference < -_box_edge_radius) ||
			(difference > _box_edge_radius))
		{
			return false;
		}
	}
	return true;
}

void RailOptimizeInfo::SetRailControlPointInfo(int _point_num, int _point_dim)
{
	num_control_point_ = _point_num;
	control_point_dim_ = _point_dim;
	dim_x_ = num_control_point_ * control_point_dim_;
}
