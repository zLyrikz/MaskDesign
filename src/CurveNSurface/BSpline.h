#pragma once
#include <array>
#include <vector>
#include <Eigen/Core>
#include <Eigen/SparseCore>
// consider all the splines in a spline space
// using B spline as basis 
// model the curve by its control points(n-dimensional)
// 
// some notation: 
// order = m
// domain = [a,b]
// inner knots = x_1,...,x_k;  (a < x_1 <...< x_k < b)
// multiplicity  m_1,...,m_k;  (1 <= m_i <= m)
// K := m_1 +...+ m_k
// all knots = y_0 <= y_1 <= ... <= y_(2m+K-1), including multiplicity
// Bspline: N_i^m(x); 
// a B spline's support covers m intervals, having m+1 knots
// 
// spline has at least C^(m - 1 - m_i) smoothness at knot x_i
// 
// space dim = m + K (in the periodic case, sapce dim = K (because we require 0 to m-1 order smoothness at domain boundary))
// to get a basis of the space, need extend 2m knots to get m+K B splines: m knots <= a and m knots >= b
// so total knots number = 2m + K
// 
// control points size = m+K; in periodic case, last m points equal to first m points (check THEOREM8.3 on P299 of Spline Functions Basic Theory (Larry Schumaker))
// 
// in the implementation of Rhino, the very first and very last knots are omitted (https://www.rhino3d.com/features/nurbs/#:~:text=Rhino%20%2D%20What%20are%20NURBS%3F,free%2Dform%20surface%20or%20solid.)
	
// defineSpace functions decide the space parameters and knot information
// fitting, interpolation functions decide space as well as the control points


//TODO: get some validation check functions, like m_i <= m, etc.
class BSpline
{
public:
	BSpline();
	~BSpline();
	BSpline(const BSpline& _other);// have to define copy constructor because there is iterator as members
	BSpline& operator=(const BSpline& _other);

public:
	// for all define(Periodic)Space function, will clear previous space info if it's defined before, and refresh to get the new one
	

	// given 3 groups of knots, 
	// each pair give the knot and its multiplicity;
	void defineSpace(
		const std::vector<std::pair<double, int>>& _inner_knots,
		const std::vector<std::pair<double, int>>& _left_extended_knots,
		const std::vector<std::pair<double, int>>& _right_extended_knots,
		int _degree, 
		const std::array<double, 2>& _domain = std::array<double, 2>{0.0, 1.0});
	
	// a frequent case, inner knots are equally spaced with same multiplicity, extended knots are choosen to be domain points a and b
	//_inner_knots_number is number of different inner knots
	void defineSpace(int _inner_knots_number, int _multiplicity, int _degree, const std::array<double, 2>& _domain = std::array<double, 2>{0.0, 1.0});
	
	// a frequent case, extened knots are choosen to be domain points a and b
	// TODO: completely this additional space definition
	//void defineSpace(const std::vector<std::pair<double, int>>& _inner_knots, int _degree, const std::array<double, 2>& _domain = std::array<double, 2>{0.0, 1.0});
	
	// input inner knots including multiplicity
	// define a periodic spline space, reference Spline Functions Basic Theory (Larry Schumaker) Chapter8
	// make sure inner knots strictly in the domain
	// extened knots are periodicly choosen from the inner knot
	// spline will be infinitely smooth at domain end points
	// notice a hard requirement K>m; 
	void definePeriodicSpace(const std::vector<double>& _inner_knots, int _degree, const std::array<double, 2>& _domain = std::array<double, 2>{0.0, 1.0});
	// uniform inner knots, same multiplicity
	//_inner_knots_number is number of different inner knots
	void definePeriodicSpace(int _inner_knots_number, int _multiplicity, int _degree, const std::array<double, 2>& _domain = std::array<double, 2>{0.0, 1.0});

	void saveToFile(const std::string& _file) const;
	void loadFromFile(const std::string& _file);
public:
	// input points size = space dim (= m + K or K for period)
	void setControlPoints(const std::vector<Eigen::VectorXd>& _control_points);
	void setControlPointDimension(int _dim_point);

	// reference Spline Functions Basic Theory (Larry Schumaker) Ch5.2
	// NOTE: need to set derivative control points not unpdated if control points or space is changed!
	void UpdateSecondOrderDerivativeControlPoints();

public:
	// TODO: need more consideration on the choosen parameters relation to the Periodic Basis!
	// given data points and polynomial segments, using uniform knots to fit a periodic spline, domain is [0,1]; TODO: this is not uniform actually, make it truly uniform?
	// different inner knot number = segments
	// parameter values are uniformly choosen
	// assign min smoothness required (continuous derivatives up to this number); 
	// max possible smoothness <= degree-1; smoothness >= -1; and knot multiplicity = degree - smoothness
	// note also make sure K > m; here K = multiplicity * knot number = (degree - smoothness) * segments; e.g. degree = 3, smoothness = 2, then knot number >= 5
	void LeastSquareFitting_Periodic(const std::vector<Eigen::VectorXd>& _points, int _segments, int _degree, int _smoothness);

	// fitting points are in several parts, each part is fitted by a polynomial, and connect them in a spline
	// different from the non-periodic implementation, parameters are uniformly selected, we select non-uniform knots, to separate different parts
	// input _parts_idx give the begin indexs of the first to last part, must be increasing order
	// (for most parts, start index < end index; yet for the last part, end index may periodicly come back to a smaller index, right before the start index of the first part)
	void LeastSquareFitting_Periodic(const std::vector<Eigen::VectorXd>& _points, const std::vector<int>& _parts_idx, int _degree, int _smoothness);


	// using Greville Interpolation from CAGD by Farin Chapter 9.1 
	// also the so called the universal method in this slides http://www.cad.zju.edu.cn/home/zhx/GM/009/00-bsia.pdf from State Key Lab of CAD&CG Zhejiang University	
	// using uniform knots in a given domain (mutiplicity = 1); adjust smoothness by changing degree, smoothness = degree - multiplicity
	// fit given data points with Greville abscissae as parameters, 
	// we consider the 0...K-1 periodic basis Greville abscissae as the 0...K-1 normal basis, 
	// so it may happen that the 0...m-1 Greville abscissae out of the domain, then simply add a period domain length to the abscissae, 
	// and the periodic basis by definition will periodicly take value as the id+K normal basis 
	// notice periodic space require K>m i.e. data point number > degree + 1
	void GrevilleInterpolation_Periodic(const std::vector<Eigen::VectorXd>& _points, int _degree, const std::array<double, 2>& _domain = std::array<double, 2>{0.0, 1.0});
	
	// given the interpolate the data points with parameter values (in strict increasing order in the domain)
	// will set inner knot exact as the parameter values
	// notice periodic space require K>m i.e. input data point number > degree + 1
	void InterpolationWithParameters_Periodic(const std::vector<Eigen::VectorXd>& _points, 
		const std::vector<double>& _parameter_values, int _degree, const std::array<double, 2>& _domain = std::array<double, 2>{0.0, 1.0}, 
		/*int* _end_index = nullptr, */Eigen::SparseMatrix<double>* _interpolation_matrix = nullptr);


	// fitting points are in several parts, each part is fitted by a polynomial, and connect them in a spline
	// part intervals equally distributed in domain [0,1]
	// for each part, parameter values are uniformly choosen in its interval
	// input _parts_idx give the begin indexs of the second to last part (note begin index of first part is always 0)
	// so different inner knot number = number of parts -1 = _parts_idx.size()
	// choose smoothness at knots, max possible smoothness <= degree-1; smoothness >= -1; and knot multiplicity = degree - smoothness
	void LeastSquareFitting(const std::vector<Eigen::VectorXd>& _points, const std::vector<int>& _parts_idx, int _degree, int _smoothness);
	// degenerate case, fitting with a polyline;space has no inner knots
	void PolynomialLeastSquareFitting(const std::vector<Eigen::VectorXd>& _points, int _degree);
	// input the known control points by a list of [index,point] pair
	// require _control_points_idx has increasing order
	void LeastSquareFitting_WithSomeControlPointsKnown(const std::vector<int>& _control_points_idx, const std::vector<Eigen::VectorXd>& _control_points,
		const std::vector<Eigen::VectorXd>& _points, const std::vector<int>& _parts_idx, int _degree, int _smoothness);
	void PolynomialLeastSquareFitting_WithSomeControlPointsKnown(const std::vector<int>& _control_points_idx, const std::vector<Eigen::VectorXd>& _control_points, 
		const std::vector<Eigen::VectorXd>& _points, int _degree);


public:
	// insert knot t in interval [y_I, y_(I+1))
	// I >= m-1; I <= m+K to make sure index used not out of range
	// input original knots, control points; knot t; interval index I; 
	// we will update the original knot and control points
	// TODO: understand what if y_I = t = y_(I+1) ?
	void KnotInsertion(std::vector<double>& _all_knots, std::vector<Eigen::VectorXd>& _control_points, int _I, double _knot) const;

	// calling knot insertion to get Bezier control points for each polynomial segment
	// for the first and last segment in the domain, we insert knot at the first extended knots rather than the domain end points
	// output points size = m * segments
	void ConvertToBeizerCurves(std::vector<Eigen::VectorXd>& _control_points) const;

public:
	// all knots: y_0 <= y_1 <= ... <= y_(2m+K-1)
	// x in interval [y_I, y_(I+1)); closed on the left, open on the right
	// output index I
	int getIntervalContainX(double _x) const;
	void getValue(double _x, Eigen::VectorXd& _value) const;// compute curve value at given parameter, using de Boor Algorithm
	void getDerivative(double _x, Eigen::VectorXd& _value) const;
	void getSecondDerivative(double _x, Eigen::VectorXd& _value);
	void getFirstDerivative(double _x, Eigen::VectorXd& _value);// unlike getDerivative, this uses the stored derivative control points for more efficiency
	// TODO: add a function, get both value and derivative
	const std::array<double, 2>& getDomain() const;
	// for nonperiod case, segments = different inner knot number + 1
	// for period case, segments = different inner knot number
	int getSegments() const;
	int getDegree() const;
	void getPolynomialIntervals(std::vector<double>& _knots) const;// return a,x_1,...,x_k,b

	// just pass the control points pointer
	// note, e.g. [const int*&] is a reference to a pointer to an integer that is const; just read the notation from right to left
	void getControlPoints(const std::vector<Eigen::VectorXd>*& _control_points) const;
	int getControlPointDimension() const;
	int getSpaceDimension() const;
public:
	// utility functions

	// construct members inner_knots_, left_extended_knots_, right_extended_knots_ from all_knots_ and order_
	void GroupKnots();

	// compute spline value in [y_I, y_(I+1))
	// related control points number = order := m
	// s(x) = d(0)*N_0^m +...+ d(m-1)*N_(m-1)^m, x in [y_I, y_(I+1))
	// knots: y_(I-m+2),..., y_(I+m-1); (note N_0^m have knots: y_(I-m+1)...y_(I+1), and N_(m-1)^m have knots: y_I...y_(I+m))
	// // make sure y_I < y_(I+1) strictly
	void deBoorAlgorithm_LocalIndex(double _x, 
		const std::vector<Eigen::VectorXd>& _control_points, 
		const std::vector<double>& _knots, 
		Eigen::VectorXd& _value) const;
	
	// basis_id can only be choosen from 0 to m+K-1
	// return true iff _x is in the basis's corresponding interval [y_I, y_(I+1))
	bool getBasisValue(int _basis_id, double _x, double& _value) const;

	// periodic basis refer to Spline Functions Basic Theory (Larry Schumaker) Chapter8
	// basis_id can only be choosen from 0 to K-1 
	// return state = 0, iff value not in corresponding intervals of the normal basis with same id and the normal basis with id = input id + K
	// return state = 1, iff value in corresponding interval of the normal basis with same id
	// return state = 2, iff value in corresponding interval of the normal basis with id = input id + K
	int getPeriodicBasisValue(int _basis_id, double _x, double& _value) const;



private:
	// 0 <= _local_basis_id < order
	// interval index I is also given (x in [y_I, y_(I+1)))
	double getBasisValue(int _local_basis_id, double _x, int _I) const;

	// id begin from 0
	// return -1 if no valid transfer
	int transferToPeriodBasisId(int _normal_basis_id, int _K, int _degree) const;

	// requires to predefine the space
	// given the interpolate the data points with parameter values (in the domain)
	// make sure the periodic basis value is nonzero at corresponding parameters
	// as for the nonperiodic case, the interpolation matrix is nonsingular, requiring each parameter value must be provided such that the corresponding Basis is nonzero (reference see Spline Functions Basic Theory (Larry Schumaker) Chapter4.8 MATRICES AND DETERMINANTS)
	// notice periodic space require K>m i.e. data point number > degree + 1
	void Interpolation_Periodic(const std::vector<Eigen::VectorXd>& _points, const std::vector<double>& _parameter_values,
		Eigen::SparseMatrix<double>* _interpolation_matrix = nullptr);

	int getMiddleKnotIdx(int _basis_id) const;

	//use in the function [LeastSquareFitting()]
	// need space to be defined already
	// // fitting points are in several parts
	void ConstructLeastSquareFittingCoefficientMatrix(int _num_data, const std::vector<int>& _parts_idx, Eigen::MatrixXd& _coefficient) const;
	void LeastSquareFittingSpaceSetUp(int _data_dim, int _num_different_inner_knots, int _degree, int _smoothness);

	// with space already defined, find control points
	// uniform parameters 
	void LeastSquareFitting_Periodic(const std::vector<Eigen::VectorXd>& _points);

	void getSecondDerivative_WithoutUpdateCheck(double _x, Eigen::VectorXd& _value) const;
	void getFirstDerivative_WithoutUpdateCheck(double _x, Eigen::VectorXd& _value) const;

private:

	// space parameter
	std::array<double, 2> domain_{0.0, 1.0};
	int degree_ = 0;
	int order_ = 0;
	int dim_space_ = 0;
	bool is_periodic_ = false;

	// all the knots we have, including multiplicity
	// extended knots are the first m and last m knots (counting multiplicity); they must be outside of the open domain (a,b)
	// y_0,...,y_(2m+K-1)
	std::vector<double> all_knots_;

	// each pair give the knot and its multiplicity;
	// iterator at the first elements position in all_knots_
	std::vector<std::pair<std::vector<double>::const_iterator, int>> inner_knots_;
	std::vector<std::pair<std::vector<double>::const_iterator, int>> left_extended_knots_;
	std::vector<std::pair<std::vector<double>::const_iterator, int>> right_extended_knots_;



	// curve parameter
	int dim_point_ = 0;
	std::vector<Eigen::VectorXd> control_points_;// size = m+K

	bool is_derivative_control_points_updated_ = false;
	// derivative_control_points_[d-2] (d>=2) is the set of control points of derivative d-1 (i.e. cd(d, *) in the Larry's book)
	// currently only updated to second order: i.e. size = 2;
	std::vector<std::vector<Eigen::VectorXd>> derivative_control_points_;

	// NOTE: change the copy and assign constructor when adding new members

};