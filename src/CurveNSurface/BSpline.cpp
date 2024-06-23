#include "BSpline.h"
#include "../Utility/VectorOperations.h"
#include "../Utility/DivisionWithReminder.h"
//#include <deque>
#include <Eigen/SVD>
#include<Eigen/SparseLU>
#include <fstream>
#include <iostream>
#include <numeric>
using std::cout;
using std::endl;

BSpline::BSpline()
{
}

BSpline::~BSpline()
{
}

BSpline::BSpline(const BSpline& _other):
	domain_(_other.domain_), degree_(_other.degree_), order_(_other.order_), dim_space_(_other.dim_space_), 
	all_knots_(_other.all_knots_), is_periodic_(_other.is_periodic_), dim_point_(_other.dim_point_), control_points_(_other.control_points_),
	is_derivative_control_points_updated_(_other.is_derivative_control_points_updated_), derivative_control_points_(_other.derivative_control_points_)
{
	GroupKnots();
}

BSpline& BSpline::operator=(const BSpline& _other)
{
	domain_ = _other.domain_;
	degree_ = _other.degree_;
	order_ = _other.order_;
	dim_space_ = _other.dim_space_;
	all_knots_ = _other.all_knots_;
	is_periodic_ = _other.is_periodic_;
	dim_point_ = _other.dim_point_;
	control_points_ = _other.control_points_;
	is_derivative_control_points_updated_ = _other.is_derivative_control_points_updated_;
	derivative_control_points_ = _other.derivative_control_points_;

	GroupKnots();

	return *this;
}


void BSpline::defineSpace(
	const std::vector<std::pair<double, int>>& _inner_knots,
	const std::vector<std::pair<double, int>>& _left_extended_knots,
	const std::vector<std::pair<double, int>>& _right_extended_knots,
	int _degree, 
	const std::array<double, 2>& _domain)
{	
	is_derivative_control_points_updated_ = false;
	is_periodic_ = false;
	degree_ = _degree;
	order_ = degree_ + 1;
	domain_ = _domain;
	

	all_knots_.clear();
	inner_knots_.clear();
	left_extended_knots_.clear();
	right_extended_knots_.clear();
	inner_knots_.reserve(_inner_knots.size());
	left_extended_knots_.reserve(_left_extended_knots.size());
	right_extended_knots_.reserve(_right_extended_knots.size());

	// record first index of the knot in all_knots, then from this index, construct the group of knots with iterator
	std::vector<int> inner_knots;
	inner_knots.reserve(_inner_knots.size());
	std::vector<int> left_extended_knots;
	left_extended_knots.reserve(_left_extended_knots.size());
	std::vector<int> right_extended_knots;
	right_extended_knots.reserve(_right_extended_knots.size());

	// construct left extended knots
	for (int iLeft = 0; iLeft < _left_extended_knots.size(); ++iLeft)
	{
		int multiplicity = _left_extended_knots[iLeft].second;
		all_knots_.reserve(all_knots_.size() + multiplicity);
		for (int iMultiplicity = 0; iMultiplicity < multiplicity; ++iMultiplicity)
		{
			all_knots_.push_back(_left_extended_knots[iLeft].first);
		}

		left_extended_knots.push_back(std::distance(all_knots_.begin(), all_knots_.end()) - multiplicity);
	}

	// construct inner knots	
	int K = 0;// for finding space dimension
	for (int iInner = 0; iInner < _inner_knots.size(); ++iInner)
	{
		int multiplicity = _inner_knots[iInner].second;
		K += multiplicity;
		all_knots_.reserve(all_knots_.size() + multiplicity);
		for (int iMultiplicity = 0; iMultiplicity < multiplicity; ++iMultiplicity)
		{
			all_knots_.push_back(_inner_knots[iInner].first);
		}

		inner_knots.push_back(std::distance(all_knots_.begin(), all_knots_.end()) - multiplicity);
	}
	dim_space_ = order_ + K;

	// construct right extended knots
	for (int iRight = 0; iRight < _right_extended_knots.size(); ++iRight)
	{
		int multiplicity = _right_extended_knots[iRight].second;
		all_knots_.reserve(all_knots_.size() + multiplicity);
		for (int iMultiplicity = 0; iMultiplicity < multiplicity; ++iMultiplicity)
		{
			all_knots_.push_back(_right_extended_knots[iRight].first);
		}

		right_extended_knots.push_back(std::distance(all_knots_.begin(), all_knots_.end()) - multiplicity);
	}

	for (int iLeft = 0; iLeft < _left_extended_knots.size(); ++iLeft)
	{
		left_extended_knots_.push_back(std::make_pair(all_knots_.begin() + left_extended_knots[iLeft], _left_extended_knots[iLeft].second));
	}
	for (int iInner = 0; iInner < _inner_knots.size(); ++iInner)
	{
		inner_knots_.push_back(std::make_pair(all_knots_.begin() + inner_knots[iInner], _inner_knots[iInner].second));
	}
	for (int iRight = 0; iRight < _right_extended_knots.size(); ++iRight)
	{
		right_extended_knots_.push_back(std::make_pair(all_knots_.begin() + right_extended_knots[iRight], _right_extended_knots[iRight].second));
	}


}

void BSpline::defineSpace(int _inner_knots_number, int _multiplicity, int _degree, const std::array<double, 2>& _domain)
{
	is_derivative_control_points_updated_ = false;

	is_periodic_ = false;
	int K = _inner_knots_number * _multiplicity;
	degree_ = _degree;
	order_ = degree_ + 1;
	dim_space_ = order_ + K;
	domain_ = _domain;

	all_knots_.clear();
	inner_knots_.clear();
	left_extended_knots_.clear();
	right_extended_knots_.clear();
	all_knots_.reserve(2 * order_ + K);
	left_extended_knots_.reserve(1);
	right_extended_knots_.reserve(1);
	inner_knots_.reserve(_inner_knots_number);

	// left extened knots
	for (int iExtend = 0; iExtend < order_; ++iExtend)
	{
		all_knots_.push_back(domain_[0]);
	}
	left_extended_knots_.push_back(std::make_pair((all_knots_.end() - order_), order_));

	// inner knots
	double partition_length = (domain_[1] - domain_[0]) / ((double)_inner_knots_number + 1.0);
	for (int iInnerKnot = 0; iInnerKnot < _inner_knots_number; ++iInnerKnot)
	{	
		double x_i = ((double)iInnerKnot + 1.0) * partition_length + domain_[0];

		for (int iMultiplicity = 0; iMultiplicity < _multiplicity; ++iMultiplicity)
		{
			all_knots_.push_back(x_i);
		}

		inner_knots_.push_back(std::make_pair((all_knots_.end() - _multiplicity), _multiplicity));
	}

	// right extened knots
	for (int iExtend = 0; iExtend < order_; ++iExtend)
	{
		all_knots_.push_back(domain_[1]);
	}
	right_extended_knots_.push_back(std::make_pair((all_knots_.end() - order_), order_));


}

void BSpline::definePeriodicSpace(const std::vector<double>& _inner_knots, int _degree, const std::array<double, 2>& _domain)
{
	is_derivative_control_points_updated_ = false;

	is_periodic_ = true;
	degree_ = _degree;
	order_ = degree_ + 1;
	domain_ = _domain;
	dim_space_ = _inner_knots.size();
	double domain_length = domain_[1] - domain_[0];

	if (dim_space_ <= order_)
	{
		cout << "WARNING from BSpline::definePeriodicSpace, K <= m; behavior undefined; stop construction" << endl;
	}
	else
	{
		std::vector<double> left_extended_knots(order_);
		std::vector<double> right_extended_knots(order_);

		for (int iMove = 0; iMove < order_; ++iMove)
		{
			// left_extended knots are a domain length translation of last m inner knots
			int correspond_inner_knot_id = (_inner_knots.size() - 1) - iMove;
			left_extended_knots[(left_extended_knots.size() - 1) - iMove] = _inner_knots[correspond_inner_knot_id] - domain_length;
		}
		for (int iMove = 0; iMove < order_; ++iMove)
		{
			// right_extended knots are a domain length translation of first m inner knots
			int correspond_inner_knot_id = iMove;
			right_extended_knots[iMove] = _inner_knots[correspond_inner_knot_id] + domain_length;
		}

		//construct all_knots_
		all_knots_.clear();
		VectorOperations concatenate;
		concatenate.AppendBack(all_knots_, left_extended_knots);
		concatenate.AppendBack(all_knots_, _inner_knots);
		concatenate.AppendBack(all_knots_, right_extended_knots);

		GroupKnots();
	}
}

void BSpline::definePeriodicSpace(int _inner_knots_number, int _multiplicity, int _degree, const std::array<double, 2>& _domain)
{
	// get all inner knots including multiplicity
	std::vector<double> inner_knots;
	inner_knots.reserve(_inner_knots_number * _multiplicity);
	double domain_length = _domain[1] - _domain[0];
	for (int iInner = 0; iInner < _inner_knots_number; ++iInner)
	{
		// reference about data cast http://c.biancheng.net/view/1329.html
		double y = (iInner + 1.0) / (_inner_knots_number + 1.0) * domain_length + _domain[0];
		for (int iMultilplicity = 0; iMultilplicity < _multiplicity; ++iMultilplicity)
		{
			inner_knots.push_back(y);
		}		
	}

	definePeriodicSpace(inner_knots, _degree, _domain);
}

void BSpline::saveToFile(const std::string& _file) const
{
	std::ofstream file(_file);
	file << "domain_: " << domain_[0] << " " << domain_[1] << "\n";
	file << "degree_: " << degree_ << "\n";
	file << "order_: " << order_ << "\n";
	file << "dim_space_: " << dim_space_ << "\n";
	file << "is_periodic_: " << is_periodic_ << "\n";
	file << "number_of_all_knots_: " << all_knots_.size() << "\n";
	file << "all_knots_: ";
	for (const auto& knot : all_knots_) {
		file << knot << " ";
	}
	file << "\n";
	file << "dim_point_: " << dim_point_ << "\n";
	if (dim_point_ > 0)
	{
		file << "control_points_: ";
		for (const auto& point : control_points_) {
			file << point.transpose() << "\n";
		}
	}
	file.close();
}

void BSpline::loadFromFile(const std::string& _file)
{
	std::ifstream file(_file);
	std::string dummy;
	file >> dummy >> domain_[0] >> domain_[1];
	file >> dummy >> degree_;
	file >> dummy >> order_;
	file >> dummy >> dim_space_;
	file >> dummy >> is_periodic_;
	int all_knot_size = 0;
	file >> dummy >> all_knot_size;
	file >> dummy;
	double knot;
	for (int iKnot = 0; iKnot < all_knot_size; ++iKnot)
	{
		file >> knot;
		all_knots_.push_back(knot);
	}

	file >> dummy >> dim_point_;
	if (dim_point_ > 0)
	{
		file >> dummy;
		double value;
		while (file >> value) {
			Eigen::VectorXd point(dim_point_);
			point[0] = value;
			for (int i = 1; i < dim_point_; ++i) {
				file >> point[i];
			}
			control_points_.push_back(point);
		}
	}

	file.close();

	is_derivative_control_points_updated_ = false;
	GroupKnots();

}

void BSpline::setControlPoints(const std::vector<Eigen::VectorXd>& _control_points)
{
	assert(_control_points.size() == dim_space_);

	if (dim_space_ > 0)
	{
		if (is_periodic_ == false)
		{
			control_points_ = _control_points;
			dim_point_ = control_points_[0].rows();
		}
		else
		{
			// extend m control points, duplicate the first m control points
			control_points_ = _control_points;
			control_points_.reserve(control_points_.size() + order_);
			for (int iPeriod = 0; iPeriod < order_; ++iPeriod)
			{
				control_points_.push_back(control_points_[iPeriod]);
			}

			dim_point_ = control_points_[0].rows();
		}
		is_derivative_control_points_updated_ = false;
	}
}

void BSpline::setControlPointDimension(int _dim_point)
{
	dim_point_ = _dim_point;
}

void BSpline::UpdateSecondOrderDerivativeControlPoints()
{
	is_derivative_control_points_updated_ = true;
	derivative_control_points_.resize(2);
	// initialize to the control points, will use it to updated the rest, and be updated
	std::vector<Eigen::VectorXd> last_derivative_control_points = control_points_;
	for (int iDerivative = 0; iDerivative < 2; ++iDerivative)
	{
		if (last_derivative_control_points.size() > 1)
		{
			std::vector<Eigen::VectorXd>& derivative_control_points = derivative_control_points_[iDerivative];
			int num_control_point = last_derivative_control_points.size() - 1;
			derivative_control_points.resize(num_control_point);

			// denote m as the last order, the current order = m-1
			int order_derivative = order_ - iDerivative - 1;// current order of the derivative spline
			int knot_start = iDerivative + 1;// start knot index of the involved knots for the derivative spline

			for (int iControl = 0; iControl < num_control_point; ++iControl)
			{
				double denominator = all_knots_[iControl + order_derivative + knot_start] - all_knots_[iControl + knot_start];// the two end knots of N_iControl^(m-1)
				if (abs(denominator) < 1e-10)
				{
					derivative_control_points[iControl] = Eigen::VectorXd::Zero(dim_point_);// make it become 0
				}
				else
				{
					derivative_control_points[iControl] =
						order_derivative * (last_derivative_control_points[iControl + 1] - last_derivative_control_points[iControl]) / denominator;
				}
			}

			// updated for the next loop
			last_derivative_control_points = derivative_control_points;
		}
	}
}

void BSpline::LeastSquareFitting_Periodic(const std::vector<Eigen::VectorXd>& _points, int _segments, int _degree, int _smoothness)
{
	// define space
	int multiplicity = _degree - _smoothness;
	definePeriodicSpace(_segments, multiplicity, _degree);
	assert(dim_space_ == multiplicity * _segments);

	// find control points
	LeastSquareFitting_Periodic(_points);

	//int num_data = _points.size();
	//if (num_data > 0)
	//{
	//	int multiplicity = _degree - _smoothness;
	//	definePeriodicSpace(_segments, multiplicity, _degree);
	//	setControlPointDimension(_points[0].rows());
	//	assert(dim_space_ == multiplicity * _segments);
	//	control_points_.resize(dim_space_ + order_, Eigen::VectorXd(dim_point_));

	//	Eigen::MatrixXd coefficient(num_data, dim_space_);
	//	// construct coefficient matrix (N_j(w_i))
	//	for (int iData = 0; iData < num_data; ++iData)
	//	{
	//		// parameter values are uniformly choosen
	//		// since it's for periodic curve, we don't use end domain point as a parameter value 
	//		double w = double(iData) / num_data * (domain_[1] - domain_[0]) + domain_[0];
	//		for (int iBasis = 0; iBasis < dim_space_; ++iBasis)
	//		{
	//			double value = 0.0;
	//			getPeriodicBasisValue(iBasis, w, value);
	//			coefficient(iData, iBasis) = value;
	//		}
	//	}
	//	//cout << "matrix\n " << coefficient << endl;

	//	// sovle data dimension least square linear systems 
	//	for (int iDim = 0; iDim < dim_point_; ++iDim)
	//	{
	//		// fitting points of current dimension
	//		Eigen::VectorXd points(num_data);
	//		for (int iData = 0; iData < num_data; ++iData)
	//		{
	//			points[iData] = _points[iData][iDim];
	//		}

	//		// solve
	//		Eigen::VectorXd control_points_i = coefficient.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(points);
	//		//cout << "result\n " << control_points_i << endl;
	//		
	//		// get control points
	//		for (int iControl = 0; iControl < dim_space_; ++iControl)
	//		{
	//			control_points_[iControl](iDim) = control_points_i(iControl);
	//		}
	//	}
	//	// since it's periodic, last m control points equal to first m
	//	for (int iPeriod = 0; iPeriod < order_; ++iPeriod)
	//	{
	//		control_points_[iPeriod + dim_space_] = control_points_[iPeriod];
	//	}
	//}
}

void BSpline::LeastSquareFitting_Periodic(const std::vector<Eigen::VectorXd>& _points, const std::vector<int>& _parts_idx, int _degree, int _smoothness)
{
	// To construct non-uniform inner knots

	// since we want the knots to be between the parts as well as inside the domain
	// in case first part begin index = 0, then the first knot would be at the left boundary of the domain
	// we periodically reorder the data points, making the first point to be the middle point of the last interval
	int num_points = _points.size();
	int last_part_start_idx = _parts_idx.back();
	int last_part_end_idx = num_points + (_parts_idx[0] - 1);
	int last_interval_middle_idx = ((last_part_end_idx + last_part_start_idx) / 2) % num_points;
	// translate the middle point of last interval to the first point
	// the correspondence is: 
	// old_idx = (new_idx + last_interval_middle_idx) % num_points
	// new_idx = (old_idx - last_interval_middle_idx) % num_points
	std::vector<int> reorderd_parts_idx(_parts_idx.size());
	for (int iPart = 0; iPart < _parts_idx.size(); ++iPart)
	{
		DivisionWithReminder find_period_idx;
		int correspond_new_idx = find_period_idx.getReminder(_parts_idx[iPart] - last_interval_middle_idx, num_points);
		reorderd_parts_idx[iPart] = correspond_new_idx;
	}

	// finally construct non-uniform inner knots
	int multiplicity = _degree - _smoothness;
	int num_segment = _parts_idx.size();
	std::vector<double> inner_knots;
	inner_knots.reserve((num_segment * multiplicity));
	for (int iPart = 0; iPart < num_segment; ++iPart)
	{
		double knot = double(reorderd_parts_idx[iPart]) / num_points;
		assert(knot > 0 && knot < 1);
		for (int iRepeat = 0; iRepeat < multiplicity; ++iRepeat)
		{
			inner_knots.push_back(knot);
		}
	}


	// define space
	definePeriodicSpace(inner_knots, _degree);
	assert(dim_space_ == multiplicity * num_segment);

	// find control points
	std::vector<Eigen::VectorXd> reorderd_points(num_points);
	for (int iPoint = 0; iPoint < num_points; ++iPoint)
	{
		int correspond_original_idx = (iPoint + last_interval_middle_idx) % num_points;
		reorderd_points[iPoint] = _points[correspond_original_idx];
	}
	LeastSquareFitting_Periodic(reorderd_points);
}

void BSpline::GrevilleInterpolation_Periodic(const std::vector<Eigen::VectorXd>& _points, int _degree, const std::array<double, 2>& _domain)
{
	int num_data = _points.size();
	assert(num_data > _degree + 1);// K>m
	std::vector<double> inner_knots;
	inner_knots.reserve(num_data);
	for (int iInner = 0; iInner < num_data; ++iInner)
	{
		double knot = double(iInner + 1) / (num_data + 1) * (_domain[1] - _domain[0]) + _domain[0];
		inner_knots.push_back(knot);
	}
	definePeriodicSpace(inner_knots, _degree, _domain);

	// get Greville abscissae
	std::vector<double> greville_abscissae;
	greville_abscissae.reserve(num_data);
	for (int iData = 0; iData < num_data; ++iData)
	{
		std::vector<double>::const_iterator knot_begin = all_knots_.begin() + (iData + 1);
		double average = std::accumulate(knot_begin, knot_begin + _degree, 0.0);//sum of the invovled knots of the basis id=iData
		average /= _degree;

		// add period to make the greville_abscissae in the domain
		if (average < _domain[0])
		{
			average += _domain[1] - _domain[0];
		}

		greville_abscissae.push_back(average);
	}

	Interpolation_Periodic(_points, greville_abscissae);
}

void BSpline::InterpolationWithParameters_Periodic(const std::vector<Eigen::VectorXd>& _points, 
	const std::vector<double>& _parameter_values, int _degree, const std::array<double, 2>& _domain,
	/*int* _end_index, */Eigen::SparseMatrix<double>* _interpolation_matrix)
{
	int num_data = _points.size();
	assert(num_data > _degree + 1);// K>m
	assert(num_data == _parameter_values.size());


	// set parameter values as knots
	definePeriodicSpace(_parameter_values, _degree, _domain);
	assert(dim_space_ == num_data);

	// now in the current spline space, parameter values are corresponded to knots: m, m+1,..., m+K-1 (i.e. inner knots)
	// perform a periodic index translation, to make each parameter become the middle knot of the periodic basis 
	// because the middle knot approximately reach the max value of the basis
	// (note the first m period basis can move its knot periodically)
	int id_last_middle_knot = getMiddleKnotIdx(dim_space_ - 1);// note last periodic basis id = K-1 = dim_space_ - 1
	assert(id_last_middle_knot >= order_ && id_last_middle_knot <= order_ + dim_space_ - 1);// it's a inner knot 
	int idx_parameter = id_last_middle_knot - order_;// get it's corresponding idx of the parameter 

	// translate this parameter corresponding to the last basis to the last
	// also do this to the data points!
	std::vector<double> new_parameter_values_back(_parameter_values.begin(), _parameter_values.begin() + (idx_parameter + 1));
	std::vector<double> new_parameter_values_front(_parameter_values.begin() + (idx_parameter + 1), _parameter_values.end());
	VectorOperations append;
	append.AppendBack(new_parameter_values_front, new_parameter_values_back);
	std::vector<Eigen::VectorXd> new_points_back(_points.begin(), _points.begin() + (idx_parameter + 1));
	std::vector<Eigen::VectorXd> new_points_front(_points.begin() + (idx_parameter + 1), _points.end());
	append.AppendBack(new_points_front, new_points_back);

	Eigen::SparseMatrix<double>* interpolation_matrix_permutated = nullptr;
	if (_interpolation_matrix != nullptr)
	{
		interpolation_matrix_permutated = new Eigen::SparseMatrix<double>;
	}

	Interpolation_Periodic(new_points_front, new_parameter_values_front, interpolation_matrix_permutated);

	if (interpolation_matrix_permutated != nullptr)
	{
		int size_permuted_front = _parameter_values.size() - (idx_parameter + 1);
		// build the interpolation matrix, if without the parameter index translation
		// so we build from the permuted one by traversing the permuted
		// if permuted.row < size_permuted_front
		//		unpermuted.row = permuted.row + (idx_parameter + 1)
		// else 
		//		unpermuted.row = permuted.row - size_permuted_front

		std::vector<Eigen::Triplet<double>> interpolation_matrix_triplets;
		for (int k = 0; k < interpolation_matrix_permutated->outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(*interpolation_matrix_permutated, k); it; ++it)
			{
				int unpermuted_row = 0;
				if (it.row() < size_permuted_front)
				{
					unpermuted_row = it.row() + (idx_parameter + 1);
				}
				else
				{
					unpermuted_row = it.row() - size_permuted_front;
				}

				interpolation_matrix_triplets.emplace_back(Eigen::Triplet<double>(unpermuted_row, it.col(), it.value()));
			}
		}

		_interpolation_matrix->resize(interpolation_matrix_permutated->rows(), interpolation_matrix_permutated->cols());
		_interpolation_matrix->setFromTriplets(interpolation_matrix_triplets.begin(), interpolation_matrix_triplets.end());

		delete interpolation_matrix_permutated;
		interpolation_matrix_permutated = nullptr;
	}

}

void BSpline::Interpolation_Periodic(const std::vector<Eigen::VectorXd>& _points, const std::vector<double>& _parameter_values,
	Eigen::SparseMatrix<double>* _interpolation_matrix)
{
	int num_data = _points.size();
	if (num_data > 0)
	{
		assert(num_data == _parameter_values.size() && num_data == dim_space_);

		std::vector<Eigen::Triplet<double>> coefficients_builder;
		coefficients_builder.reserve(num_data * (2 * order_ - 1));
		for (int iRow = 0; iRow < num_data; ++iRow)
		{
			double value = 0.0;
			int state = getPeriodicBasisValue(iRow, _parameter_values[iRow], value);
			coefficients_builder.push_back(Eigen::Triplet<double>(iRow, iRow, value));

			// for this banded coefficient matrix
			// fill the columns near center element, semi band width = degree
			// find the normal basis index that current periodic basis behaves like
			// then get column values near this center
			int corresponding_normal_basis_id = iRow;// this is the actual center id
			if (state == 2)
			{
				corresponding_normal_basis_id += dim_space_;
			}
			for (int iNeighbor = 1; iNeighbor <= order_-1; ++iNeighbor)
			{
				int left_neighbor_id = corresponding_normal_basis_id - iNeighbor;
				int right_neighbor_id = corresponding_normal_basis_id + iNeighbor;
				int left_col = transferToPeriodBasisId(left_neighbor_id, dim_space_, order_-1);
				int right_col = transferToPeriodBasisId(right_neighbor_id, dim_space_, order_-1);
				if (left_col != -1)
				{
					getPeriodicBasisValue(left_col, _parameter_values[iRow], value);
					coefficients_builder.push_back(Eigen::Triplet<double>(iRow, left_col, value));
				}
				if (right_col != -1)
				{
					getPeriodicBasisValue(right_col, _parameter_values[iRow], value);
					coefficients_builder.push_back(Eigen::Triplet<double>(iRow, right_col, value));
				}
			}
		}

		// compute matrix and build solver
		Eigen::SparseMatrix<double> coefficients(num_data, num_data);
		coefficients.setFromTriplets(coefficients_builder.begin(), coefficients_builder.end());
		if (_interpolation_matrix != nullptr)
		{
			*_interpolation_matrix = coefficients;
		}
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		solver.compute(coefficients);
		if (solver.info() != Eigen::Success) {
			cout << "BSpline::Interpolation_Periodic Eigen decomposition failed" << endl;
			return;
		}

		// get unknown control points and the right hand set data points
		setControlPointDimension(_points[0].rows());
		control_points_.resize(dim_space_ + order_, Eigen::VectorXd(dim_point_));
		for (int iDim = 0; iDim < dim_point_; ++iDim)
		{
			Eigen::VectorXd control_points(num_data);// the unknown control points of current dimension
			Eigen::VectorXd points(num_data);// right hand side data of current dimension
			for (int iData = 0; iData < num_data; ++iData)
			{
				points(iData) = _points[iData](iDim);
			}
			// solve
			control_points = solver.solve(points);
			if (solver.info() != Eigen::Success) {
				cout << "BSpline::Interpolation_Periodic Eigen solving failed" << endl;
				return;
			}

			//// test linear system hold
			//Eigen::VectorXd y = coefficients * control_points;
			//for (int iData = 0; iData < num_data; ++iData)
			//{
			//	if (abs(y(iData) - points(iData)) > 1e-10)
			//	{
			//		cout << "abs(y(iData) - points(iData))=" << abs(y(iData) - points(iData)) << endl;
			//	}
			//}

			// get control points
			for (int iPoint = 0; iPoint < num_data; ++iPoint)
			{
				control_points_[iPoint](iDim) = control_points(iPoint);
			}
		}
		// since it's periodic, last m control points equal to first m
		for (int iPeriod = 0; iPeriod < order_; ++iPeriod)
		{
			control_points_[iPeriod + dim_space_] = control_points_[iPeriod];
		}

		//// test interpolation
		//for (int iData = 0; iData < num_data; ++iData)
		//{
		//	Eigen::VectorXd interpolated;
		//	getValue(_parameter_values[iData], interpolated);
		//	Eigen::VectorXd difference = interpolated - _points[iData];
		//	for (int iDim = 0; iDim < dim_point_; ++iDim)
		//	{
		//		if (abs(difference(iDim)) > 1e-10)
		//		{
		//			cout << "iData=" << iData << ", difference(" << iDim << ") = " << difference(iDim) << endl;
		//		}
		//	}
		//}
	}
}

int BSpline::getMiddleKnotIdx(int _basis_id) const
{
	int middle = _basis_id + order_ / 2.0;// notice if order is odd, then the middle is taken as the left knot 
	return middle;
}

void BSpline::ConstructLeastSquareFittingCoefficientMatrix(int _num_data, const std::vector<int>& _parts_idx, Eigen::MatrixXd& _coefficient) const
{
	// construct coefficient matrix (N_j(w_i))
	_coefficient.resize(_num_data, dim_space_);
	int num_part = _parts_idx.size() + 1;
	for (int iPart = 0; iPart < num_part; ++iPart)
	{
		double local_domain_length = (domain_[1] - domain_[0]) / num_part;
		double local_domain[2] = { domain_[0] + local_domain_length * iPart, domain_[0] + local_domain_length * (iPart + 1) };

		// for the first part, parameters reach local domain left end value, and not the right end value
		// for the middle parts, parameters not reach any of the 2 local domain values
		// for the last part, parameters reach local domain right end value, and not the left end value
		if (iPart == 0 && iPart != num_part - 1)// first part but not the last part, so _parts_idx[0] is safe to use
		{
			int num_point = _parts_idx[0];
			for (int iPoint = 0; iPoint < num_point; ++iPoint)
			{
				// parameter values are uniformly choosen
				double w = double(iPoint) / num_point * (local_domain[1] - local_domain[0]) + local_domain[0];
				for (int iBasis = 0; iBasis < dim_space_; ++iBasis)
				{
					double value = 0.0;
					getBasisValue(iBasis, w, value);
					_coefficient(iPoint, iBasis) = value;
				}
			}
		}
		else if (iPart < (num_part - 1))// middle parts
		{
			assert(iPart > 0);
			int num_point = _parts_idx[iPart] - _parts_idx[iPart - 1];
			for (int iPoint = 0; iPoint < num_point; ++iPoint)
			{
				// parameter values are uniformly choosen
				double w = double(iPoint + 1) / (num_point + 1) * (local_domain[1] - local_domain[0]) + local_domain[0];
				for (int iBasis = 0; iBasis < dim_space_; ++iBasis)
				{
					double value = 0.0;
					getBasisValue(iBasis, w, value);
					_coefficient(_parts_idx[iPart - 1] + iPoint, iBasis) = value;
				}
			}
		}
		else // iPart == (num_part - 1)
		{
			assert(iPart == (num_part - 1));
			int last_part_begin_index = 0;
			if (iPart > 0)// else, the part_idx is empty
			{
				last_part_begin_index = _parts_idx[iPart - 1];
			}
			else
			{
				assert(_parts_idx.empty());
			}

			int num_point = _num_data - last_part_begin_index;
			for (int iPoint = 0; iPoint < num_point; ++iPoint)
			{
				// parameter values are uniformly choosen
				double w = double(iPoint + 1) / num_point * (local_domain[1] - local_domain[0]) + local_domain[0];
				for (int iBasis = 0; iBasis < dim_space_; ++iBasis)
				{
					double value = 0.0;
					getBasisValue(iBasis, w, value);
					_coefficient(last_part_begin_index + iPoint, iBasis) = value;
				}
			}
		}
	}
}

void BSpline::LeastSquareFittingSpaceSetUp(int _data_dim, int _num_different_inner_knots, int _degree, int _smoothness)
{
	int multiplicity = _degree - _smoothness;
	defineSpace(_num_different_inner_knots, multiplicity, _degree);

	setControlPointDimension(_data_dim);
	assert(dim_space_ == multiplicity * _num_different_inner_knots + order_);
	control_points_.resize(dim_space_, Eigen::VectorXd(dim_point_));
}

void BSpline::LeastSquareFitting_Periodic(const std::vector<Eigen::VectorXd>& _points)
{
	// about Eigen 
	// least square linear systems https://eigen.tuxfamily.org/dox/group__LeastSquares.html#:~:text=An%20overdetermined%20system%20of%20equations,the%20Euclidean%20norm%20is%20used).
	// solving linear systems https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
	//
	// for BSpline, Schoenberg-Whitney theorem: tells when coefficient matrix M^T*M is singular

	int num_data = _points.size();
	if (num_data > 0 && dim_space_ > 0)
	{
		setControlPointDimension(_points[0].rows());
		control_points_.resize(dim_space_ + order_, Eigen::VectorXd(dim_point_));

		Eigen::MatrixXd coefficient(num_data, dim_space_);
		// construct coefficient matrix (N_j(w_i))
		for (int iData = 0; iData < num_data; ++iData)
		{
			// parameter values are uniformly choosen
			// since it's for periodic curve, we don't use end domain point as a parameter value 
			double w = double(iData) / num_data * (domain_[1] - domain_[0]) + domain_[0];
			for (int iBasis = 0; iBasis < dim_space_; ++iBasis)
			{
				double value = 0.0;
				getPeriodicBasisValue(iBasis, w, value);
				coefficient(iData, iBasis) = value;
			}
		}
		//cout << "matrix\n " << coefficient << endl;

		// sovle data dimension least square linear systems 
		for (int iDim = 0; iDim < dim_point_; ++iDim)
		{
			// fitting points of current dimension
			Eigen::VectorXd points(num_data);
			for (int iData = 0; iData < num_data; ++iData)
			{
				points[iData] = _points[iData][iDim];
			}

			// solve
			Eigen::VectorXd control_points_i = coefficient.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(points);
			//cout << "result\n " << control_points_i << endl;

			// get control points
			for (int iControl = 0; iControl < dim_space_; ++iControl)
			{
				control_points_[iControl](iDim) = control_points_i(iControl);
			}
		}
		// since it's periodic, last m control points equal to first m
		for (int iPeriod = 0; iPeriod < order_; ++iPeriod)
		{
			control_points_[iPeriod + dim_space_] = control_points_[iPeriod];
		}
	}
}

void BSpline::getSecondDerivative_WithoutUpdateCheck(double _x, Eigen::VectorXd& _value) const
{
	if (order_ > 2)
	{
		int I = getIntervalContainX(_x);

		// for the Bpline of the i-th derivative(order = m-i)
		// for second derivative i=2
		// last associated Bspline index is I, and there are m-i BSpline associated in total
		// so the first associated BSpline index is I-(m-i) + 1
		// 
		// since only the control points for the Bspline that have nonzero value in the domain is computed
		// the first associated control point index is (I-(m-i) + 1) - i = I-m+1, last is I - i; here i means i-th derivative
		

		// get involved control points
		// notice that when I = m-1(the first index of interval in the domain), control point index start at 0
		std::vector<Eigen::VectorXd> control_points(derivative_control_points_[1].begin() + (I - order_ + 1), derivative_control_points_[1].begin() + (I - 1));
		// not include the first and last associated knots 
		std::vector<double> knots(all_knots_.begin() + (I - (order_ - 2) + 2), all_knots_.begin() + (I + (order_ - 2)));
		deBoorAlgorithm_LocalIndex(_x, control_points, knots, _value);
	}
	else
	{
		_value = Eigen::VectorXd::Zero(dim_point_);
	}
}

void BSpline::getFirstDerivative_WithoutUpdateCheck(double _x, Eigen::VectorXd& _value) const
{
	if (order_ > 1)
	{
		int I = getIntervalContainX(_x);

		// check comment in getSecondDerivative_WithoutUpdateCheck

		// get involved control points
		// notice that when I = m-1(the first index of interval in the domain), control point index start at 0
		std::vector<Eigen::VectorXd> control_points(derivative_control_points_[0].begin() + (I - order_ + 1), derivative_control_points_[0].begin() + I);
		// not include the first and last associated knots 
		std::vector<double> knots(all_knots_.begin() + (I - (order_ - 1) + 2), all_knots_.begin() + (I + (order_ - 1)));
		deBoorAlgorithm_LocalIndex(_x, control_points, knots, _value);
	}
	else
	{
		_value = Eigen::VectorXd::Zero(dim_point_);
	}
}

void BSpline::LeastSquareFitting(const std::vector<Eigen::VectorXd>& _points, const std::vector<int>& _parts_idx, int _degree, int _smoothness)
{
	int num_data = _points.size();
	if (num_data > 0)
	{
		//int multiplicity = _degree - _smoothness;
		//defineSpace(_parts_idx.size(), multiplicity, _degree);
		//
		//setControlPointDimension(_points[0].rows());
		//assert(dim_space_ == multiplicity * _parts_idx.size() + order_);
		//control_points_.resize(dim_space_, Eigen::VectorXd(dim_point_));


		//// construct coefficient matrix (N_j(w_i))
		//Eigen::MatrixXd coefficient(num_data, dim_space_);	
		//int num_part = _parts_idx.size() + 1;
		//for (int iPart = 0; iPart < num_part; ++iPart)
		//{
		//	double local_domain_length = (domain_[1] - domain_[0]) / num_part;
		//	double local_domain[2] = { domain_[0] + local_domain_length * iPart, domain_[0] + local_domain_length * (iPart + 1) };
		//
		//	// for the first part, parameters reach local domain left end value, and not the right end value
		//	// for the middle parts, parameters not reach any of the 2 local domain values
		//	// for the last part, parameters reach local domain right end value, and not the left end value
		//	if (iPart == 0 && iPart != num_part - 1)// first part but not the last part, so _parts_idx[0] is safe to use
		//	{
		//		int num_point = _parts_idx[0];
		//		for (int iPoint = 0; iPoint < num_point; ++iPoint)
		//		{
		//			// parameter values are uniformly choosen
		//			double w = double(iPoint) / num_point * (local_domain[1] - local_domain[0]) + local_domain[0];
		//			for (int iBasis = 0; iBasis < dim_space_; ++iBasis)
		//			{
		//				double value = 0.0;
		//				getBasisValue(iBasis, w, value);
		//				coefficient(iPoint, iBasis) = value;
		//			}
		//		}
		//	}
		//	else if (iPart < (num_part - 1))// middle parts
		//	{
		//		assert(iPart > 0);
		//		int num_point = _parts_idx[iPart] - _parts_idx[iPart - 1];
		//		for (int iPoint = 0; iPoint < num_point; ++iPoint)
		//		{
		//			// parameter values are uniformly choosen
		//			double w = double(iPoint + 1) / (num_point + 1) * (local_domain[1] - local_domain[0]) + local_domain[0];
		//			for (int iBasis = 0; iBasis < dim_space_; ++iBasis)
		//			{
		//				double value = 0.0;
		//				getBasisValue(iBasis, w, value);
		//				coefficient(_parts_idx[iPart - 1] + iPoint, iBasis) = value;
		//			}
		//		}
		//	}
		//	else // iPart == (num_part - 1)
		//	{
		//		assert(iPart == (num_part - 1));
		//		int num_point = num_data - _parts_idx[iPart - 1];
		//		for (int iPoint = 0; iPoint < num_point; ++iPoint)
		//		{
		//			// parameter values are uniformly choosen
		//			double w = double(iPoint + 1) / num_point * (local_domain[1] - local_domain[0]) + local_domain[0];
		//			for (int iBasis = 0; iBasis < dim_space_; ++iBasis)
		//			{
		//				double value = 0.0;
		//				getBasisValue(iBasis, w, value);
		//				coefficient(_parts_idx[iPart - 1] + iPoint, iBasis) = value;
		//			}
		//		}
		//	}
		//}
		////cout << "matrix\n " << coefficient << endl;

		LeastSquareFittingSpaceSetUp(_points[0].rows(), _parts_idx.size(), _degree, _smoothness);

		Eigen::MatrixXd coefficient;
		ConstructLeastSquareFittingCoefficientMatrix(num_data, _parts_idx, coefficient);


		// sovle data dimension many least square linear systems 
		for (int iDim = 0; iDim < dim_point_; ++iDim)
		{
			// fitting points of current dimension
			Eigen::VectorXd points(num_data);
			for (int iData = 0; iData < num_data; ++iData)
			{
				points[iData] = _points[iData][iDim];
			}

			// solve
			Eigen::VectorXd control_points_i = coefficient.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(points);
			//cout << "result\n " << control_points_i << endl;

			// get control points
			for (int iControl = 0; iControl < dim_space_; ++iControl)
			{
				control_points_[iControl](iDim) = control_points_i(iControl);
			}
		}
	}
}

void BSpline::PolynomialLeastSquareFitting(const std::vector<Eigen::VectorXd>& _points, int _degree)
{
	std::vector<int> part_idx;// an empty vector
	LeastSquareFitting(_points, part_idx, _degree, 0);
}

void BSpline::LeastSquareFitting_WithSomeControlPointsKnown(const std::vector<int>& _control_points_idx, const std::vector<Eigen::VectorXd>& _control_points,
	const std::vector<Eigen::VectorXd>& _points, const std::vector<int>& _parts_idx, int _degree, int _smoothness)
{
	const int num_data = _points.size();
	if (num_data > 0)
	{
		LeastSquareFittingSpaceSetUp(_points[0].rows(), _parts_idx.size(), _degree, _smoothness);

		// get a general coefficent matrix, then edit it to fit the current problem
		Eigen::MatrixXd general_coefficient;
		ConstructLeastSquareFittingCoefficientMatrix(num_data, _parts_idx, general_coefficient);

		// get unknown control point index (find the rest of the index not in the [_control_points_idx])
		int num_unknown_control_points = dim_space_ - _control_points.size();
		assert(num_unknown_control_points);
		assert(_control_points_idx.size() == _control_points.size());
		std::vector<int> unknown_control_points_idx;
		unknown_control_points_idx.reserve(num_unknown_control_points);
		{
			// use the increasing order of the [_control_points_idx], find rest of the index
			int counter_known = 0;
			for (int iControl = 0; iControl < dim_space_; ++iControl)
			{
				if (counter_known < _control_points_idx.size() && iControl == _control_points_idx[counter_known])
				{
					// current iControl index is in [_control_points_idx]
					++counter_known;
				}
				else
				{
					unknown_control_points_idx.push_back(iControl);
				}
			}
		}

		Eigen::MatrixXd coefficient(num_data, num_unknown_control_points);
		coefficient = general_coefficient(Eigen::placeholders::all, unknown_control_points_idx);// see matrix indexing document https://eigen.tuxfamily.org/dox-devel/group__TutorialSlicingIndexing.html
		//cout << "matrix\n " << coefficient << endl;
		//cout << "general_coefficient.rows()*general_coefficient.cols()=" << general_coefficient.rows() << "," << general_coefficient.cols() << endl;
		//cout << "coefficient.rows()*coefficient.cols()=" << coefficient.rows() << "," << coefficient.cols() << endl;

		// sovle data dimension many least square linear systems 
		for (int iDim = 0; iDim < dim_point_; ++iDim)
		{
			// fitting points of current dimension
			Eigen::VectorXd points(num_data);
			for (int iData = 0; iData < num_data; ++iData)
			{
				points[iData] = _points[iData][iDim];
				// minus the known control points
				for (int iKnown = 0; iKnown < _control_points.size(); ++iKnown)
				{
					points[iData] -= general_coefficient(iData, _control_points_idx[iKnown]) * _control_points[iKnown](iDim);
				}
			}

			// solve
			Eigen::VectorXd control_points_i = coefficient.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(points);
			//cout << "result\n " << control_points_i << endl;

			// get control points
			for (int iUnknown = 0; iUnknown < unknown_control_points_idx.size(); ++iUnknown)
			{
				control_points_[unknown_control_points_idx[iUnknown]](iDim) = control_points_i(iUnknown);
			}
		}
		for (int iKnown = 0; iKnown < _control_points_idx.size(); ++iKnown)
		{
			control_points_[_control_points_idx[iKnown]] = _control_points[iKnown];
		}
	}
}

void BSpline::PolynomialLeastSquareFitting_WithSomeControlPointsKnown(const std::vector<int>& _control_points_idx, const std::vector<Eigen::VectorXd>& _control_points, const std::vector<Eigen::VectorXd>& _points, int _degree)
{
	LeastSquareFitting_WithSomeControlPointsKnown(_control_points_idx, _control_points,
		_points, std::vector<int>(), _degree, 0);
}

void BSpline::KnotInsertion(std::vector<double>& _all_knots, std::vector<Eigen::VectorXd>& _control_points, int _I, double _knot) const
{
	assert(_I >= order_ - 1 && _I <= _all_knots.size() - order_);

	// duplicate the control_point at I-m+1
	_control_points.insert(_control_points.begin() + (_I - order_ + 1), _control_points[_I - order_ + 1]);
	// one time deBoor at knot t
	for (int i = 0; i < order_ - 1; ++i)
	{
		int current_knot_idx = _I - order_ + 2 + i;
		double a = (_knot - _all_knots[current_knot_idx]) / (_all_knots[current_knot_idx + order_ - 1] - _all_knots[current_knot_idx]);
		_control_points[current_knot_idx] = (1.0 - a) * _control_points[current_knot_idx] + a * _control_points[current_knot_idx + 1];
	}
	// notice now, _control_points[I+1] is previously the _control_points[I](the index before de Boor operations)

	_all_knots.insert(_all_knots.begin() + (_I + 1), _knot);
}

void BSpline::ConvertToBeizerCurves(std::vector<Eigen::VectorXd>& _control_points) const
{
	// get knot index and mutiplicity
	std::vector<std::pair<int, int>> knots_info;// pair := (index, multiplicity)
	knots_info.reserve(inner_knots_.size() + 2);

	// the left extended one
	knots_info.push_back(std::make_pair(order_ - 1, left_extended_knots_.back().second));
	// the inners
	for (int iInner = 0; iInner < inner_knots_.size(); ++iInner)
	{
		int multiplicity = inner_knots_[iInner].second;
		int last_idx = std::distance(all_knots_.begin(), inner_knots_[iInner].first)/* + multiplicity - 1*/;// TODO: see if using first index works
		knots_info.push_back(std::make_pair(last_idx, multiplicity));
	}
	// the right extended one
	knots_info.push_back(std::make_pair(all_knots_.size() - order_, right_extended_knots_[0].second));

	// for each knot, insert the same knots to reach multiplicity m
	std::vector<double> all_knots = all_knots_;
	std::vector<Eigen::VectorXd> control_points = control_points_;
	// insert sequence from right to left, so the index isn't affected
	for (int iKnot = knots_info.size() - 1 ; iKnot >= 0; --iKnot)
	{
		
		int num_to_insert = order_ - knots_info[iKnot].second;
		for (int iInsert = 0; iInsert < num_to_insert; ++iInsert)
		{
			KnotInsertion(all_knots, control_points, knots_info[iKnot].first, all_knots[knots_info[iKnot].first]);
		}
	}

	// get Bezier Points
	int first_point_idx = std::distance(all_knots_.begin(), left_extended_knots_.back().first);
	_control_points.clear();
	int num_bezier_point = (inner_knots_.size() + 1) * order_;// the points at the joint of two Beizer curves counted twice
	_control_points.reserve((inner_knots_.size() + 1) * order_);// notice we store two same points at the joint of two Beizer curves
	for (int iPoint = 0; iPoint < num_bezier_point; ++iPoint)
	{
		_control_points.push_back(control_points[first_point_idx + iPoint]);
	}
}

int BSpline::getIntervalContainX(double _x) const
{
	int id_interval = 0;

	while (id_interval < inner_knots_.size() && * (inner_knots_[id_interval].first) <= _x)
	{
		++id_interval;
	}
	// now inner_knots_[id_interval] > _x; inner_knots_[id_interval - 1] <= _x (overshoot one inner knot to pass the multiplicity)

	int I = 0;
	if (inner_knots_.size() == 0)
	{
		I = order_ - 1;
	}
	else if (id_interval == inner_knots_.size())
	{
		I = std::distance(all_knots_.begin(), inner_knots_.back().first) + (inner_knots_.back().second - 1);// end index of the inner knot in all knots		
	}
	else
	{
		I = std::distance(all_knots_.begin(), inner_knots_[id_interval].first) - 1;
	}
	return I;
}

// TODO: consider when _x becomes one knot
void BSpline::getValue(double _x, Eigen::VectorXd& _value) const
{
	// de Boor Algorithm (refer to CAGD book by Farin)
	// find the interval, _x in [y_I, y_(I+1)); closed on the left, open on the right
	// 
	// there are m Bspline associated with the interval, therefore m control points (same with Bspline order)
	// need m-1 iteration; N^m -> N^(m-1) -> ... -> N^1

	int I = getIntervalContainX(_x);

	//    the last associated BSpline index is I, and there are m BSpline associated in total
	// so the first associated BSpline index is I-(m-1)
	std::vector<Eigen::VectorXd> control_points(control_points_.begin() + (I - order_ + 1), control_points_.begin() + (I + 1));// used for iteration
	
	std::vector<double> knots(all_knots_.cbegin() + (I - order_ + 2), all_knots_.cbegin() + (I + order_));
	deBoorAlgorithm_LocalIndex(_x, control_points, knots, _value);

	//// explicit way of implementation
	//for (int i = 0; i < degree_; ++i)
	//{
	//	// m-i control points left
	//	// control_points[j] belong to N_j^(m-i), j=0,...,m-i-1; (notice here j is some kind of local index for the B Spline)
	//	// the (m-i-j)th knot of N_j^(m-i) is y_I, so N_j^(m-i) involving knots: y_[I-(m-i-j)+1],...,y_(I+j+1)
	//	// control_points[j] = s[_x<i>, y_[I-(m-i-j)+2],...,y_(I+j)]; (this is the blossom representaion of the spline in interval [y_I, y_(I+1)))
	//	// need (m-i)-1 affine combination operations to get next set of control points
	//	for (int j = 0; j < degree_ - i; ++j)
	//	{
	//		// affine combination of y_[I-(m-i-j)+2] and y_(I+j+1)
	//		// _x = (1-a) * y_[I-(m-i-j)+2] + a * y_(I+j+1)
	//		// notice denominator is never 0 because y_I < y_(I+1)
	//		double a = (_x - all_knots_[I - (order_ - i - j) + 2]) / (all_knots_[I + j + 1] - all_knots_[I - (order_ - i - j) + 2]);
	//		control_points[j] = (1.0 - a) * control_points[j] + a * control_points[j + 1];
	//	}
	//}
	//_value = control_points[0];
}

void BSpline::getDerivative(double _x, Eigen::VectorXd& _value) const
{
	// reference Spline Functions Basic Theory (Larry Schumaker)
	int I = getIntervalContainX(_x);
	int order_derivative = order_ - 1;

	//    the last associated BSpline index is I, and there are m BSpline associated in total
	// so the first associated BSpline index is I-(m-1)
	// get involved knots: the knots of all lower order Bspline N^(m-1), size = 2*(order_ - 1)
	std::vector<double> knots(all_knots_.begin() + (I - (order_ - 1) + 1), all_knots_.begin() + (I + (order_ - 1) + 1));

	// get involved control points of N^m; size = order
	std::vector<Eigen::VectorXd> control_points(control_points_.begin() + (I - (order_ - 1)), control_points_.begin() + (I + 1));

	// compute control points of the derivative
	for (int iControl = 0; iControl < order_derivative; ++iControl)
	{
		double denominator = knots[iControl + order_derivative] - knots[iControl];// the two end knots of N_iControl^(m-1)
		if (abs(denominator) < 1e-10)
		{
			control_points[iControl] = Eigen::VectorXd::Zero(dim_point_);// make it become 0
		}
		else
		{
			control_points[iControl] = order_derivative * (control_points[iControl + 1] - control_points[iControl]) / denominator;
		}
	}

	control_points.pop_back();
	knots.erase(knots.begin());
	knots.pop_back();
	deBoorAlgorithm_LocalIndex(_x, control_points, knots, _value);
}

void BSpline::getSecondDerivative(double _x, Eigen::VectorXd& _value)
{
	if (!is_derivative_control_points_updated_)
	{
		UpdateSecondOrderDerivativeControlPoints();
	}

	getSecondDerivative_WithoutUpdateCheck(_x, _value);
}

void BSpline::getFirstDerivative(double _x, Eigen::VectorXd& _value)
{
	if (!is_derivative_control_points_updated_)
	{
		UpdateSecondOrderDerivativeControlPoints();
	}

	getFirstDerivative_WithoutUpdateCheck(_x, _value);
}

const std::array<double, 2>& BSpline::getDomain() const
{
	return domain_;
}

int BSpline::getSegments() const
{
	int segments = inner_knots_.size();
	if (is_periodic_ == false)
	{
		segments += 1;
	}
	return segments;
}

int BSpline::getDegree() const
{
	return degree_;
}

void BSpline::getPolynomialIntervals(std::vector<double>& _knots) const
{
	_knots.resize(inner_knots_.size() + 2);
	_knots[0] = domain_[0];
	for (int iKnot = 0; iKnot < inner_knots_.size(); ++iKnot)
	{
		_knots[iKnot + 1] = *(inner_knots_[iKnot].first);
	}
	_knots.back() = domain_[1];
}

void BSpline::getControlPoints(const std::vector<Eigen::VectorXd>*& _control_points) const
{
	_control_points = &control_points_;
}

int BSpline::getControlPointDimension() const
{
	return dim_point_;
}

int BSpline::getSpaceDimension() const
{
	return dim_space_;
}

void BSpline::GroupKnots()
{
	if (all_knots_.size() > 0)
	{
		int current_idx = 0;// record progress
		const double EPSILON = 1e-10;
		inner_knots_.clear();
		left_extended_knots_.clear();
		right_extended_knots_.clear();

		// get left knots
		while (all_knots_[current_idx] < domain_[0] + EPSILON)
		{
			// find next set of points with same values 
			int multiplicity = 1;
			while (abs(all_knots_[current_idx + multiplicity - 1] - all_knots_[current_idx + multiplicity]) < EPSILON)
			{
				++multiplicity;
			}
			// note [current_idx] hasn't been updated

			std::vector<double>::const_iterator first_iter = all_knots_.begin() + current_idx;
			left_extended_knots_.push_back(std::make_pair(first_iter, multiplicity));

			current_idx += multiplicity;// update current position, now it's the index of the beginning of a new set of same value points
		}
		assert(current_idx == order_);

		// get inner knots
		while (all_knots_[current_idx] < domain_[1] - EPSILON)
		{
			// find next set of points with same values 
			int multiplicity = 1;
			while (abs(all_knots_[current_idx + multiplicity - 1] - all_knots_[current_idx + multiplicity]) < EPSILON)
			{
				++multiplicity;
			}
			// note [current_idx] hasn't been updated

			std::vector<double>::const_iterator first_iter = all_knots_.begin() + current_idx;
			inner_knots_.push_back(std::make_pair(first_iter, multiplicity));

			current_idx += multiplicity;// update current position, now it's the index of the beginning of a new set of same value points
		}
		if (is_periodic_)
		{
			assert(current_idx == order_ + dim_space_);
		}
		else
		{
			assert(current_idx == dim_space_);
		}

		// get right_extended_knots
		while (current_idx < all_knots_.size())
		{
			// find next set of points with same values 
			int multiplicity = 1;
			while ((current_idx + multiplicity) < all_knots_.size() &&
				abs(all_knots_[current_idx + multiplicity - 1] - all_knots_[current_idx + multiplicity]) < EPSILON)
			{
				++multiplicity;
			}
			// note [current_idx] hasn't been updated

			std::vector<double>::const_iterator first_iter = all_knots_.begin() + current_idx;
			right_extended_knots_.push_back(std::make_pair(first_iter, multiplicity));

			current_idx += multiplicity;// update current position, now it's the index of the beginning of a new set of same value points
		}
		assert(current_idx == all_knots_.size());
	}
}

void BSpline::deBoorAlgorithm_LocalIndex(double _x, 
	const std::vector<Eigen::VectorXd>& _control_points, 
	const std::vector<double>& _knots, 
	Eigen::VectorXd& _value) const
{
	// TODO: figure out when knot mulpilicity larger than order!!!
	// de Boor Algorithm (refer to CAGD book by Farin)
	// find the interval, _x in [y_I, y_(I+1)); closed on the left, open on the right
	// 
	// there are m Bspline associated with the interval, therefore m control points (same with Bspline order)
	// need m-1 iteration; N^m -> N^(m-1) -> ... -> N^1

	int order = _control_points.size();
	int degree = order - 1;
	if (order >= 1)
	{
		// order = control points number := m
		// s(x) = d(0)*N_0^m +...+ d(m-1)*N_(m-1)^m, x in [y_I, y_(I+1))
		// knots: y_(I-m+2),..., y_(I+m-1), size = 2(m-1)
		// 
		int I = order - 2;// index of input _knots; (remap I-m+2 -> 0); although when order=1, I=-1, this implementation still works
		std::vector<Eigen::VectorXd> control_points(_control_points);

		for (int i = 0; i < degree; ++i)
		{
			// m-i control points left
			// control_points[j] belong to N_j^(m-i), j=0,...,m-i-1; (notice here j is some kind of local index for the B Spline)
			// the (m-i-j)th knot of N_j^(m-i) is y_I, so N_j^(m-i) involving knots: y_[I-(m-i-j)+1],...,y_(I+j+1)
			// control_points[j] = s[_x<i>, y_[I-(m-i-j)+2],...,y_(I+j)]; (this is the blossom representaion of the spline in interval [y_I, y_(I+1)))
			// need (m-i)-1 affine combination operations to get next set of control points
			for (int j = 0; j < degree - i; ++j)
			{
				// affine combination of y_[I-(m-i-j)+2] and y_(I+j+1)
				// _x = (1-a) * y_[I-(m-i-j)+2] + a * y_(I+j+1)
				// notice denominator is never 0 because y_I < y_(I+1)
				double a = (_x - _knots[i + j]) / (_knots[I + j + 1] - _knots[i + j]);
				control_points[j] = (1.0 - a) * control_points[j] + a * control_points[j + 1];
			}
		}

		_value = control_points[0];
	}
	else
	{
		cout << "input control points size = 0" << endl;
	}
}

bool BSpline::getBasisValue(int _basis_id, double _x, double& _value) const
{
	bool is_in_interval = false;
	_value = 0.0;

	// x in [y_I, y_(I+1))
	int I = getIntervalContainX(_x);

	int basis_id_local = _basis_id - (I - order_ + 1);
	if (basis_id_local >= 0 && basis_id_local < order_)// current basis is non zero in [y_I, y_(I+1))
	{
		_value = getBasisValue(basis_id_local, _x, I);
		is_in_interval = true;
	}
	
	return is_in_interval;
}

int BSpline::getPeriodicBasisValue(int _basis_id, double _x, double& _value) const
{
	int state = 0;
	_value = 0.0;

	if (is_periodic_)
	{
		// x in [y_I, y_(I+1))
		int I = getIntervalContainX(_x);

		// [y_I, y_(I+1)) involve basis: N_j, j = I-m+1,...,I
		// remap them to local index j -> j-(I-m+1) (i.e. N_j, j = 0,...,m-1)
		int basis_id_local = _basis_id - (I - order_ + 1);
		if (basis_id_local >= 0 && basis_id_local < order_)// current basis is non zero in [y_I, y_(I+1))
		{
			_value = getBasisValue(basis_id_local, _x, I);
			state = 1;
		}
		else if (_basis_id <= order_ && basis_id_local < 0) // these basis have more definition like it's index + K
		{
			basis_id_local += dim_space_;
			if (basis_id_local >= 0 && basis_id_local < order_)// current basis is non zero in [y_I, y_(I+1))
			{
				_value = getBasisValue(basis_id_local, _x, I);
				state = 2;
			}
		}

	}

	return state;
}

double BSpline::getBasisValue(int _local_basis_id, double _x, int _I) const
{
	// [y_I, y_(I+1)) involve basis: N_j, j = I-m+1,...,I
	// remap them to local index j -> j-(I-m+1) (i.e. N_j, j = 0,...,m-1)
	std::vector<Eigen::VectorXd> control_points(order_, Eigen::VectorXd::Zero(1));
	control_points[_local_basis_id] = Eigen::VectorXd::Ones(1);

	std::vector<double> knots_local(all_knots_.cbegin() + (_I - order_ + 2), all_knots_.cbegin() + (_I + order_));

	Eigen::VectorXd value;
	deBoorAlgorithm_LocalIndex(_x, control_points, knots_local, value);

	return value[0];
}

int BSpline::transferToPeriodBasisId(int _normal_basis_id, int _K, int _degree) const
{
	int new_id = -1;

	if (_normal_basis_id >= 0)
	{
		if (_normal_basis_id <= _K - 1)
		{
			new_id = _normal_basis_id;
		}
		else if (_normal_basis_id <= _degree + _K)
		{
			new_id = _normal_basis_id - _K;// this periodic basis got nonzero value
		}
	}

	return new_id;
}
