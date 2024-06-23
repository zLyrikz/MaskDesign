#include "CrossSectionSurface.h"
#include "CrossSectionHelpr.h"
#include "../CurveNSurface/SamplePolylineFromCurve.h"
#include<Eigen/SparseLU>

#include <fstream>
#include <iostream>
using std::cout;
using std::endl;
CrossSectionSurface::CrossSectionSurface()
{
}

CrossSectionSurface::~CrossSectionSurface()
{
}

void CrossSectionSurface::setControlCurve(
	const std::vector<BSpline>& _half1_rest, const BSpline& _half1_point1, 
	const BSpline& _tangent, 
	const BSpline& _half2_point1, const std::vector<BSpline>& _half2_rest)
{
	half1_rest_ = _half1_rest;
	half1_point1_length_ = _half1_point1;
	tangent_ = _tangent;
	half2_point1_length_ = _half2_point1;
	half2_rest_ = _half2_rest;

	domain_u_ = tangent_.getDomain();// other curve domain should be the same
}

void CrossSectionSurface::setControlCurve_Degree2(const BSpline& _half1_rest, const BSpline& _half1_point1_length, const BSpline& _tangent, const BSpline& _half2_point1_length, const BSpline& _half2_rest)
{
	half1_rest_.resize(1);
	half2_rest_.resize(1);
	half1_rest_[0] = _half1_rest;
	half1_point1_length_ = _half1_point1_length;
	tangent_ = _tangent;
	half2_point1_length_ = _half2_point1_length;
	half2_rest_[0] = _half2_rest;

	domain_u_ = tangent_.getDomain();// other curve domain should be the same

}

void CrossSectionSurface::setControlCurveByInterpolation(const std::vector<const CrossSection*>& _cross_sections, const std::vector<double>& _u, int _degree, const std::array<double, 2>& _domain_u)
{

	int num_section = _cross_sections.size();
	assert(num_section == _u.size());
	if (num_section > 0)
	{
		domain_u_ = _domain_u;

		// get lists of interpolation points from key cross section list
		int num_rest_half1 = _cross_sections[0]->getHalf1Degree() - 1;
		std::vector<std::vector<Eigen::VectorXd>> half1_rest(num_rest_half1, std::vector<Eigen::VectorXd>(num_section, Eigen::VectorXd(2)));
		std::vector<Eigen::VectorXd> half1_point1_length(num_section, Eigen::VectorXd(1));
		std::vector<Eigen::VectorXd> tangent_angle(num_section, Eigen::VectorXd(1));
		std::vector<Eigen::VectorXd> point0(num_section, Eigen::VectorXd(2));
		std::vector<Eigen::VectorXd> half2_point1_length(num_section, Eigen::VectorXd(1));
		int num_rest_half2 = _cross_sections[0]->getHalf2Degree() - 1;
		std::vector<std::vector<Eigen::VectorXd>> half2_rest(num_rest_half2, std::vector<Eigen::VectorXd>(num_section, Eigen::VectorXd(2)));


		for (int iSection = 0; iSection < num_section; ++iSection)
		{
			// get control points
			std::vector<Eigen::Vector2d> current_half1_rest;
			double current_half1_point1_length = 0.0;
			double current_tangent_angle;
			double current_half2_point1_length = 0.0;
			std::vector<Eigen::Vector2d> current_half2_rest;
			_cross_sections[iSection]->getControlPoints(
				current_half1_rest, current_half1_point1_length,
				current_tangent_angle,
				current_half2_point1_length, current_half2_rest);

			// put them in the list
			for (int iRest = 0; iRest < num_rest_half1; ++iRest)
			{
				half1_rest[iRest][iSection] = std::move(current_half1_rest[iRest]);
			}
			for (int iRest = 0; iRest < num_rest_half2; ++iRest)
			{
				half2_rest[iRest][iSection] = std::move(current_half2_rest[iRest]);
			}
			half1_point1_length[iSection](0) = current_half1_point1_length;
			half2_point1_length[iSection](0) = current_half2_point1_length;
			tangent_angle[iSection](0) = current_tangent_angle;
			point0[iSection] = *_cross_sections[iSection]->GetPoint0();
		}

		// interpolate each of them

		// we just use half1_point1_length_ here to get key_interpolation_matrix_ (it's nothing special with other control curves)
		if (key_u_.size() == _u.size()) // for now, this will make sure that each element of the two vectors are the same
		{
			for (int iU = 0; iU < _u.size(); ++iU)
			{
				assert(key_u_[iU] == _u[iU]);
			}
			half1_point1_length_.InterpolationWithParameters_Periodic(half1_point1_length, _u, _degree, _domain_u);
		}
		else
		{
			key_u_ = _u;
			half1_point1_length_.InterpolationWithParameters_Periodic(half1_point1_length, _u, _degree, _domain_u, &key_interpolation_matrix_);

		}
		half2_point1_length_.InterpolationWithParameters_Periodic(half2_point1_length, _u, _degree, _domain_u);
		tangent_angle_.InterpolationWithParameters_Periodic(tangent_angle, _u, _degree, _domain_u);
		point0_.InterpolationWithParameters_Periodic(point0, _u, _degree, _domain_u);

		half1_rest_.resize(num_rest_half1);
		half2_rest_.resize(num_rest_half2);
		for (int iRest = 0; iRest < num_rest_half1; ++iRest)
		{
			half1_rest_[iRest].InterpolationWithParameters_Periodic(half1_rest[iRest], _u, _degree, _domain_u);

		}
		for (int iRest = 0; iRest < num_rest_half2 - 1; ++iRest)// note, not include the last one, it's the boundary
		{
			half2_rest_[iRest].InterpolationWithParameters_Periodic(half2_rest[iRest], _u, _degree, _domain_u);
		}
		///*******************************************************************************************************/
		//// for debug
		///*******************************************************************************************************/
		//SamplePolylineFromCurve get_curve_mesh;
		//get_curve_mesh.getPolylineAsPointCloud2D(half1_rest_[0], 1000, control2_, _domain_u);
		//get_curve_mesh.getPolylineAsPointCloud2D(half1_rest_[1], 1000, control3_, _domain_u);
		//get_curve_mesh.getPolylineAsPointCloud2D(half2_rest_[0], 1000, control22_, _domain_u);
		//get_curve_mesh.getPolylineAsPointCloud2D(half2_rest_[1], 1000, control23_, _domain_u);
	}
}

void CrossSectionSurface::setBoundaryControlCurveByInterpolation(const std::vector<Eigen::VectorXd>& _interpolate_points, const std::vector<double>& _u, int _degree, const std::array<double, 2>& _domain_u)
{
	if (!half2_rest_.empty())
	{
		half2_rest_.back().InterpolationWithParameters_Periodic(_interpolate_points, _u, _degree, _domain_u);
	}
}

void CrossSectionSurface::setHalf1BoundaryControlCurveByInterpolation(const std::vector<Eigen::VectorXd>& _interpolate_points, const std::vector<double>& _u, int _degree, const std::array<double, 2>& _domain_u)
{
	if (!half1_rest_.empty())
	{
		half1_rest_.back().InterpolationWithParameters_Periodic(_interpolate_points, _u, _degree, _domain_u);
	}
}

void CrossSectionSurface::ComputeSampledValues(const std::vector<double>& _u, const std::vector<double>& _v)
{
	cross_section_values_.resize(_u.size());

	for (int iU = 0; iU < _u.size(); ++iU)
	{
		cross_section_values_[iU].resize(_v.size());
		for (int iV = 0; iV < _v.size(); ++iV)
		{
			getValue(_u[iU], _v[iV], cross_section_values_[iU][iV]);
		}
	}

	if (_u.size() != u_values_.size())
	{
		u_values_ = _u;
		is_derivative_updated_ = false;
	}
	if (_v.size() != v_values_.size())
	{
		v_values_ = _v;
		is_derivative_updated_ = false;
	}
}

void CrossSectionSurface::UpdateDerivativeInfo()
{
	findDerivative_ControlCurve_KeyControlPoints(u_values_);
	UpdateSampledDerivatives();
	findDerivative_CrossSectionSurfaceAndNormal_ControlCurves(u_values_, v_values_, cross_section_derivatives_);
	is_derivative_updated_ = true;
}

void CrossSectionSurface::UpdateSampledDerivatives()
{
	cross_section_derivatives_.resize(u_values_.size());
	for (int iU = 0; iU < u_values_.size(); ++iU)
	{
		// at u
		// get the 2d cross section bezier curve
		// then get the 2d derivative curve
		double u = u_values_[iU];
		BezierCurve half1;
		getBezierCurve_Half1(u, half1);
		BezierCurve derivative_half1;
		half1.getDerivativeCurve(derivative_half1);
		BezierCurve half2;
		getBezierCurve_Half2(u, half2);
		BezierCurve derivative_half2;
		half2.getDerivativeCurve(derivative_half2);

		cross_section_derivatives_[iU].resize(v_values_.size());
		for (int iV = 0; iV < v_values_.size(); ++iV)
		{
			Eigen::VectorXd derivative;// 2d vector, in the local frame

			if (v_values_[iV] < 1.0)
			{
				derivative_half1.getValue(v_values_[iV], derivative);
			}
			else
			{
				derivative_half2.getValue(v_values_[iV] - 1.0, derivative);
			}

			cross_section_derivatives_[iU][iV](0) = derivative(0);
			cross_section_derivatives_[iU][iV](1) = derivative(1);
		}
	}

}

void CrossSectionSurface::findDerivative_CrossSectionSurfaceAndNormal_KeyControlPoints(
	std::vector<std::vector<std::vector<std::array<Eigen::Matrix2Xd, 2>>>>& _half1_rest,
	std::vector<std::vector<std::vector<std::array<Eigen::Matrix2Xd, 2>>>>& _half2_rest,
	std::vector<std::vector<Eigen::Matrix2Xd>>& _half1_point1_length,
	std::vector<std::vector<Eigen::Matrix2Xd>>& _tangent_angle,
	std::vector<std::vector<Eigen::Matrix2Xd>>& _half2_point1_length,
	std::vector<std::vector<std::array<Eigen::Matrix2Xd, 2>>>& _point0,
	std::vector<std::vector<std::vector<std::array<Eigen::Matrix2Xd, 2>>>>& _normal_half1_rest,
	std::vector<std::vector<std::vector<std::array<Eigen::Matrix2Xd, 2>>>>& _normal_half2_rest,
	std::vector<std::vector<Eigen::Matrix2Xd>>& _normal_half1_point1_length,
	std::vector<std::vector<Eigen::Matrix2Xd>>& _normal_tangent_angle,
	std::vector<std::vector<Eigen::Matrix2Xd>>& _normal_half2_point1_length,
	std::vector<std::vector<std::array<Eigen::Matrix2Xd, 2>>>& _normal_point0) const
{
	_half1_rest.resize(half1_rest_.size(), std::vector<std::vector<std::array<Eigen::Matrix2Xd, 2>>>(u_values_.size(), std::vector<std::array<Eigen::Matrix2Xd, 2>>(v_values_.size())));
	_half2_rest.resize(half2_rest_.size() - 1, std::vector<std::vector<std::array<Eigen::Matrix2Xd, 2>>>(u_values_.size(), std::vector<std::array<Eigen::Matrix2Xd, 2>>(v_values_.size())));
	_half1_point1_length.resize(u_values_.size(), std::vector<Eigen::Matrix2Xd>(v_values_.size()));
	_tangent_angle.resize(u_values_.size(), std::vector<Eigen::Matrix2Xd>(v_values_.size()));
	_half2_point1_length.resize(u_values_.size(), std::vector<Eigen::Matrix2Xd>(v_values_.size()));
	_point0.resize(u_values_.size(), std::vector<std::array<Eigen::Matrix2Xd, 2>>(v_values_.size()));

	_normal_half1_rest.resize(half1_rest_.size(), std::vector<std::vector<std::array<Eigen::Matrix2Xd, 2>>>(u_values_.size(), std::vector<std::array<Eigen::Matrix2Xd, 2>>(v_values_.size())));
	_normal_half2_rest.resize(half2_rest_.size() - 1, std::vector<std::vector<std::array<Eigen::Matrix2Xd, 2>>>(u_values_.size(), std::vector<std::array<Eigen::Matrix2Xd, 2>>(v_values_.size())));
	_normal_half1_point1_length.resize(u_values_.size(), std::vector<Eigen::Matrix2Xd>(v_values_.size()));
	_normal_tangent_angle.resize(u_values_.size(), std::vector<Eigen::Matrix2Xd>(v_values_.size()));
	_normal_half2_point1_length.resize(u_values_.size(), std::vector<Eigen::Matrix2Xd>(v_values_.size()));
	_normal_point0.resize(u_values_.size(), std::vector<std::array<Eigen::Matrix2Xd, 2>>(v_values_.size()));

	for (int iU = 0; iU < u_values_.size(); ++iU)
	{
		const Eigen::VectorXd& derivative_curve2keys = control_curve_derivative_key_control_points_[iU];
		for (int iV = 0; iV < v_values_.size(); ++iV)
		{
			for (int iRest1 = 0; iRest1 < half1_rest_.size(); ++iRest1)
			{
				CombineDerivative_CrossSectionSurfaceAndNormal_ControlCurves_And_ControlCurve_KeyControlPoints(
					_half1_rest[iRest1][iU][iV][0], cross_section_surface_derivative_half1_rest_[iRest1][iU][iV][0], derivative_curve2keys);
				CombineDerivative_CrossSectionSurfaceAndNormal_ControlCurves_And_ControlCurve_KeyControlPoints(
					_normal_half1_rest[iRest1][iU][iV][0], cross_section_surface_normal_derivative_half1_rest_[iRest1][iU][iV][0], derivative_curve2keys);
				CombineDerivative_CrossSectionSurfaceAndNormal_ControlCurves_And_ControlCurve_KeyControlPoints(
					_half1_rest[iRest1][iU][iV][1], cross_section_surface_derivative_half1_rest_[iRest1][iU][iV][1], derivative_curve2keys);
				CombineDerivative_CrossSectionSurfaceAndNormal_ControlCurves_And_ControlCurve_KeyControlPoints(
					_normal_half1_rest[iRest1][iU][iV][1], cross_section_surface_normal_derivative_half1_rest_[iRest1][iU][iV][1], derivative_curve2keys);

			}
			for (int iRest2 = 0; iRest2 < cross_section_surface_normal_derivative_half2_rest_.size(); ++iRest2)
			{
				CombineDerivative_CrossSectionSurfaceAndNormal_ControlCurves_And_ControlCurve_KeyControlPoints(
					_half2_rest[iRest2][iU][iV][0], cross_section_surface_derivative_half2_rest_[iRest2][iU][iV][0], derivative_curve2keys);
				CombineDerivative_CrossSectionSurfaceAndNormal_ControlCurves_And_ControlCurve_KeyControlPoints(
					_normal_half2_rest[iRest2][iU][iV][0], cross_section_surface_normal_derivative_half2_rest_[iRest2][iU][iV][0], derivative_curve2keys);
				CombineDerivative_CrossSectionSurfaceAndNormal_ControlCurves_And_ControlCurve_KeyControlPoints(
					_half2_rest[iRest2][iU][iV][1], cross_section_surface_derivative_half2_rest_[iRest2][iU][iV][1], derivative_curve2keys);
				CombineDerivative_CrossSectionSurfaceAndNormal_ControlCurves_And_ControlCurve_KeyControlPoints(
					_normal_half2_rest[iRest2][iU][iV][1], cross_section_surface_normal_derivative_half2_rest_[iRest2][iU][iV][1], derivative_curve2keys);
			}
			CombineDerivative_CrossSectionSurfaceAndNormal_ControlCurves_And_ControlCurve_KeyControlPoints(
				_half1_point1_length[iU][iV], cross_section_surface_derivative_half1_point1_length_[iU][iV], derivative_curve2keys);
			CombineDerivative_CrossSectionSurfaceAndNormal_ControlCurves_And_ControlCurve_KeyControlPoints(
				_tangent_angle[iU][iV], cross_section_surface_derivative_tangent_angle_[iU][iV], derivative_curve2keys);
			CombineDerivative_CrossSectionSurfaceAndNormal_ControlCurves_And_ControlCurve_KeyControlPoints(
				_half2_point1_length[iU][iV], cross_section_surface_derivative_half2_point1_length_[iU][iV], derivative_curve2keys);
			CombineDerivative_CrossSectionSurfaceAndNormal_ControlCurves_And_ControlCurve_KeyControlPoints(
				_point0[iU][iV][0], cross_section_surface_derivative_point0_[iU][iV][0], derivative_curve2keys);
			CombineDerivative_CrossSectionSurfaceAndNormal_ControlCurves_And_ControlCurve_KeyControlPoints(
				_point0[iU][iV][1], cross_section_surface_derivative_point0_[iU][iV][1], derivative_curve2keys);

			CombineDerivative_CrossSectionSurfaceAndNormal_ControlCurves_And_ControlCurve_KeyControlPoints(
				_normal_half1_point1_length[iU][iV], cross_section_surface_normal_derivative_half1_point1_length_[iU][iV], derivative_curve2keys);
			CombineDerivative_CrossSectionSurfaceAndNormal_ControlCurves_And_ControlCurve_KeyControlPoints(
				_normal_tangent_angle[iU][iV], cross_section_surface_normal_derivative_tangent_angle_[iU][iV], derivative_curve2keys);
			CombineDerivative_CrossSectionSurfaceAndNormal_ControlCurves_And_ControlCurve_KeyControlPoints(
				_normal_half2_point1_length[iU][iV], cross_section_surface_normal_derivative_half2_point1_length_[iU][iV], derivative_curve2keys);
			CombineDerivative_CrossSectionSurfaceAndNormal_ControlCurves_And_ControlCurve_KeyControlPoints(
				_normal_point0[iU][iV][0], cross_section_surface_normal_derivative_point0_[iU][iV][0], derivative_curve2keys);
			CombineDerivative_CrossSectionSurfaceAndNormal_ControlCurves_And_ControlCurve_KeyControlPoints(
				_normal_point0[iU][iV][1], cross_section_surface_normal_derivative_point0_[iU][iV][1], derivative_curve2keys);
		}
	}
}

const std::vector<Eigen::VectorXd>* CrossSectionSurface::getControlCurveDerivativeKeyControlPoints() const
{
	return &control_curve_derivative_key_control_points_;
}

void CrossSectionSurface::findDerivative_ControlCurve_KeyControlPoints(const std::vector<double>& _u)
{
	if (!is_derivative_updated_)
	{
		Eigen::SparseMatrix<double> basis_vectors;
		// since any control curve has same derivative
		// we take half1_point1_length_, it's nothing special with other control curves
		FindBSplineBasisVector(basis_vectors, half1_point1_length_, _u);

		// [N(u_i)]M^(-1) is the derivative for all control curves
		// solve M^T*X=[N(u_i)]^T=basis_vector at u_i

		// construct solver
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		solver.compute(key_interpolation_matrix_.transpose());
		if (solver.info() != Eigen::Success) {
			cout << "[ERROR!!!!!!!!! CrossSectionSurface::findDerivative_ControlCurve_KeyControlPoints] Eigen decomposition failed" << endl;
			return;
		}

		// solve
		control_curve_derivative_key_control_points_.resize(_u.size());
		for (int iU = 0; iU < _u.size(); ++iU)
		{
			control_curve_derivative_key_control_points_[iU] = solver.solve(basis_vectors.col(iU));
			if (solver.info() != Eigen::Success) {
				cout << "[ERROR!!!!!!!!! CrossSectionSurface::findDerivative_ControlCurve_KeyControlPoints] Eigen solving failed" << endl;
				return;
			}
		}
	}
}

void CrossSectionSurface::findDerivative_CrossSectionSurfaceAndNormal_ControlCurves(const std::vector<double>& _u, const std::vector<double>& _v, 
	const std::vector<std::vector<Eigen::Vector2d>>& _cross_section_derivatives)
{
	// for aquiring basis values
	// degree = half.size() + 1
	BezierCurve half1(half1_rest_.size() + 1);
	BezierCurve half2(half2_rest_.size() + 1);

	// cross-section surface: p(u,v) = X_i(u) * Y_i(v)
	// Dp/Dc = DX_i/Dc * Y_i
	// Dp'/Dc = DX_i/Dc * Y_i'
	// here Y_i is just some beizer curves
	// 

	// Dp/Dc
	cross_section_surface_derivative_half1_rest_.clear();
	cross_section_surface_derivative_half1_rest_.resize(half1_rest_.size(), std::vector<std::vector<std::array<Eigen::Vector2d, 2>>>(_u.size(), std::vector<std::array<Eigen::Vector2d, 2>>(_v.size())));
	cross_section_surface_derivative_half2_rest_.clear();
	cross_section_surface_derivative_half2_rest_.resize(half2_rest_.size() - 1, std::vector<std::vector<std::array<Eigen::Vector2d, 2>>>(_u.size(), std::vector<std::array<Eigen::Vector2d, 2>>(_v.size())));
	cross_section_surface_derivative_point0_.clear();
	cross_section_surface_derivative_point0_.resize(_u.size(), std::vector<std::array<Eigen::Vector2d, 2>>(_v.size()));
	cross_section_surface_derivative_half1_point1_length_.clear();
	cross_section_surface_derivative_half1_point1_length_.resize(_u.size(), std::vector<Eigen::Vector2d>(_v.size()));
	cross_section_surface_derivative_tangent_angle_.clear();
	cross_section_surface_derivative_tangent_angle_.resize(_u.size(), std::vector<Eigen::Vector2d>(_v.size()));
	cross_section_surface_derivative_half2_point1_length_.clear();
	cross_section_surface_derivative_half2_point1_length_.resize(_u.size(), std::vector<Eigen::Vector2d>(_v.size()));

	// normal_derivative = rotation of D(p'/||p'||)/Dc
	// D(p'/||p'||)/Dc = D(x/||x||)/Dx * Dp'/Dc
	cross_section_surface_normal_derivative_half1_rest_.clear();
	cross_section_surface_normal_derivative_half1_rest_.resize(half1_rest_.size(), std::vector<std::vector<std::array<Eigen::Vector2d, 2>>>(_u.size(), std::vector<std::array<Eigen::Vector2d, 2>>(_v.size())));
	cross_section_surface_normal_derivative_half2_rest_.clear();
	cross_section_surface_normal_derivative_half2_rest_.resize(half2_rest_.size() - 1, std::vector<std::vector<std::array<Eigen::Vector2d, 2>>>(_u.size(), std::vector<std::array<Eigen::Vector2d, 2>>(_v.size())));
	cross_section_surface_normal_derivative_half1_point1_length_.clear();
	cross_section_surface_normal_derivative_half1_point1_length_.resize(_u.size(), std::vector<Eigen::Vector2d>(_v.size()));
	cross_section_surface_normal_derivative_tangent_angle_.clear();
	cross_section_surface_normal_derivative_tangent_angle_.resize(_u.size(), std::vector<Eigen::Vector2d>(_v.size()));
	cross_section_surface_normal_derivative_half2_point1_length_.clear();
	cross_section_surface_normal_derivative_half2_point1_length_.resize(_u.size(), std::vector<Eigen::Vector2d>(_v.size()));
	cross_section_surface_normal_derivative_point0_.clear();
	cross_section_surface_normal_derivative_point0_.resize(_u.size(), std::vector<std::array<Eigen::Vector2d, 2>>(_v.size()));

	// DX_i/Dc: dirivative of control curves without multiplying with the basis value
	// first get DX_i/Dc, then multiply them with basis value or its derivative to get the derivative or normal_derivative
	for (int iU = 0; iU < _u.size(); ++iU)
	{
		// get tangle angle curve value at current u
		// note the tangent direction is pointing to half2 point1
		Eigen::VectorXd tangent_angle;
		tangent_angle_.getValue(_u[iU], tangent_angle);
		Eigen::Vector2d tangent;
		CrossSectionHelpr convert_tangent;
		convert_tangent.TangentAngle2Vector(tangent_angle(0), tangent);

		for (int iV = 0; iV < _v.size(); ++iV)
		{
			// initialize everything to zero
			for (int iRest = 0; iRest < half2_rest_.size() - 1; ++iRest)
			{
				cross_section_surface_derivative_half2_rest_[iRest][iU][iV][0].setZero(2);
				cross_section_surface_derivative_half2_rest_[iRest][iU][iV][1].setZero(2);
				
				cross_section_surface_normal_derivative_half2_rest_[iRest][iU][iV][0].setZero(2);
				cross_section_surface_normal_derivative_half2_rest_[iRest][iU][iV][1].setZero(2);

			}
			for (int iRest = 0; iRest < half1_rest_.size(); ++iRest)
			{
				cross_section_surface_derivative_half1_rest_[iRest][iU][iV][0].setZero(2);
				cross_section_surface_derivative_half1_rest_[iRest][iU][iV][1].setZero(2);
				cross_section_surface_normal_derivative_half1_rest_[iRest][iU][iV][0].setZero(2);
				cross_section_surface_normal_derivative_half1_rest_[iRest][iU][iV][1].setZero(2);
			}
			cross_section_surface_derivative_half1_point1_length_[iU][iV].setZero(2);
			cross_section_surface_derivative_half2_point1_length_[iU][iV].setZero(2);
			cross_section_surface_normal_derivative_half1_point1_length_[iU][iV].setZero(2);
			cross_section_surface_normal_derivative_half2_point1_length_[iU][iV].setZero(2);
			// tangent_angle common coefficients for half1 and 2
			// take derivative of tangent w.r.t. angle
			cross_section_surface_derivative_tangent_angle_[iU][iV](0) = -tangent(1);
			cross_section_surface_derivative_tangent_angle_[iU][iV](1) = tangent(0);
			cross_section_surface_normal_derivative_tangent_angle_[iU][iV] = cross_section_surface_derivative_tangent_angle_[iU][iV];

			// point0 is common for both halfs, initialize to DX_i/Dc
			cross_section_surface_derivative_point0_[iU][iV][0].x() = 1.0;
			cross_section_surface_derivative_point0_[iU][iV][0].y() = 0.0;
			cross_section_surface_derivative_point0_[iU][iV][1].y() = 1.0;
			cross_section_surface_derivative_point0_[iU][iV][1].x() = 0.0;
			cross_section_surface_normal_derivative_point0_[iU][iV] = cross_section_surface_derivative_point0_[iU][iV];

			// get jacobian of x/||x||, take x=p'(u,v)
			Eigen::Matrix2d jacobian_normed_x;
			double temp = -_cross_section_derivatives[iU][iV].x() * _cross_section_derivatives[iU][iV].y();
			jacobian_normed_x <<
				_cross_section_derivatives[iU][iV].y() * _cross_section_derivatives[iU][iV].y(), temp,
				temp, _cross_section_derivatives[iU][iV].x() * _cross_section_derivatives[iU][iV].x();
			jacobian_normed_x /= pow(_cross_section_derivatives[iU][iV].squaredNorm(), 1.5);

			// find Dp/Dc and D(x/||x||)/Dx * Dp'/Dc(stored in "normal_derivative")
			if (_v[iV] >= 1.0)// notice 1.0 is in the v_values, so it's very different whether 1.0 is included in half1 or half2!!!!!!
			{
				// current point on half2

				double local_v = _v[iV] - 1.0;

				// half2 rest (not include the end point)
				for (int iRest = 0; iRest < half2_rest_.size() - 1; ++iRest)
				{
					int basis_id = iRest + 2;
					double bezier_value = half2.getBasisValue(basis_id, local_v);
					double bezier_derivative = half2.getBasisDerivative(basis_id, local_v);

					// derivative w.r.t. x element
					cross_section_surface_derivative_half2_rest_[iRest][iU][iV][0].x() = bezier_value;
					cross_section_surface_normal_derivative_half2_rest_[iRest][iU][iV][0].x() = bezier_derivative;
					cross_section_surface_normal_derivative_half2_rest_[iRest][iU][iV][0] = jacobian_normed_x * cross_section_surface_normal_derivative_half2_rest_[iRest][iU][iV][0];

					// derivative w.r.t. y element
					cross_section_surface_derivative_half2_rest_[iRest][iU][iV][1].y() = bezier_value;
					cross_section_surface_normal_derivative_half2_rest_[iRest][iU][iV][1].y() = bezier_derivative;
					cross_section_surface_normal_derivative_half2_rest_[iRest][iU][iV][1] = jacobian_normed_x * cross_section_surface_normal_derivative_half2_rest_[iRest][iU][iV][1];

				}

				// half2_point1_length_
				int half2_point1_length_basis_id = 1;
				double half2_point1_length_bezier_value = half2.getBasisValue(half2_point1_length_basis_id, local_v);
				double half2_point1_length_bezier_derivative = half2.getBasisDerivative(half2_point1_length_basis_id, local_v);
				cross_section_surface_derivative_half2_point1_length_[iU][iV] = tangent * half2_point1_length_bezier_value;
				cross_section_surface_normal_derivative_half2_point1_length_[iU][iV] = jacobian_normed_x * tangent * half2_point1_length_bezier_derivative;

				// tangent_angle
				// it is at the same basis as length_
				Eigen::VectorXd half2_point1_length;
				half2_point1_length_.getValue(_u[iU], half2_point1_length);
				cross_section_surface_derivative_tangent_angle_[iU][iV] *= half2_point1_length(0) * half2_point1_length_bezier_value;
				cross_section_surface_normal_derivative_tangent_angle_[iU][iV] *= half2_point1_length(0) * half2_point1_length_bezier_derivative;
				cross_section_surface_normal_derivative_tangent_angle_[iU][iV] = jacobian_normed_x * cross_section_surface_normal_derivative_tangent_angle_[iU][iV];

				// point0
				double basis0 = half2.getBasisValue(0, local_v);
				double basis0_derivative = half2.getBasisDerivative(0, local_v);
				double basis0and1 = basis0 + half2_point1_length_bezier_value;
				double basis0and1_derivative = basis0_derivative + half2_point1_length_bezier_derivative;
				cross_section_surface_derivative_point0_[iU][iV][0] *= basis0and1;
				cross_section_surface_derivative_point0_[iU][iV][1] *= basis0and1;
				cross_section_surface_normal_derivative_point0_[iU][iV][0] *= basis0and1_derivative;
				cross_section_surface_normal_derivative_point0_[iU][iV][1] *= basis0and1_derivative;
				cross_section_surface_normal_derivative_point0_[iU][iV][0] = jacobian_normed_x * cross_section_surface_normal_derivative_point0_[iU][iV][0];
				cross_section_surface_normal_derivative_point0_[iU][iV][1] = jacobian_normed_x * cross_section_surface_normal_derivative_point0_[iU][iV][1];
			}
			else
			{
				// current point on half1

				// this this equivalent to: set local_v = 1-_v[iV], and reverse the basis id
				double local_v = _v[iV];

				// half1 rest
				for (int iRest = 0; iRest < half1_rest_.size(); ++iRest)
				{
					int basis_id = half1_rest_.size() - 1 - iRest;
					double bezier_value = half1.getBasisValue(basis_id, local_v);
					double bezier_derivative = half1.getBasisDerivative(basis_id, local_v);

					// derivative w.r.t. x element of the control curve
					cross_section_surface_derivative_half1_rest_[iRest][iU][iV][0].x() = bezier_value;
					cross_section_surface_normal_derivative_half1_rest_[iRest][iU][iV][0].x() = bezier_derivative;
					cross_section_surface_normal_derivative_half1_rest_[iRest][iU][iV][0] = jacobian_normed_x * cross_section_surface_normal_derivative_half1_rest_[iRest][iU][iV][0];

					// derivative w.r.t. y element of the control curve
					cross_section_surface_derivative_half1_rest_[iRest][iU][iV][1].y() = bezier_value;
					cross_section_surface_normal_derivative_half1_rest_[iRest][iU][iV][1].y() = bezier_derivative;
					cross_section_surface_normal_derivative_half1_rest_[iRest][iU][iV][1] = jacobian_normed_x * cross_section_surface_normal_derivative_half1_rest_[iRest][iU][iV][1];

				}

				// half1_point1_length_
				int half1_point1_length_basis_id = half1_rest_.size();
				double half1_point1_length_bezier_value = half1.getBasisValue(half1_point1_length_basis_id, local_v);
				double half1_point1_length_bezier_derivative = half1.getBasisDerivative(half1_point1_length_basis_id, local_v);
				cross_section_surface_derivative_half1_point1_length_[iU][iV] = -tangent * half1_point1_length_bezier_value;
				cross_section_surface_normal_derivative_half1_point1_length_[iU][iV] = jacobian_normed_x * (- tangent)  * half1_point1_length_bezier_derivative;

				// tangent angle
				Eigen::VectorXd half1_point1_length;
				half1_point1_length_.getValue(_u[iU], half1_point1_length);
				cross_section_surface_derivative_tangent_angle_[iU][iV] *= -half1_point1_length(0) * half1_point1_length_bezier_value;
				cross_section_surface_normal_derivative_tangent_angle_[iU][iV] *= -half1_point1_length(0) * half1_point1_length_bezier_derivative;
				cross_section_surface_normal_derivative_tangent_angle_[iU][iV] = jacobian_normed_x * cross_section_surface_normal_derivative_tangent_angle_[iU][iV];

				// point0
				// notice the basis0 corresponds to the last basis in the profile curve direction
				double basis0 = half1.getBasisValue(half1_rest_.size() + 1, local_v);
				double basis0_derivative = half1.getBasisDerivative(half1_rest_.size() + 1, local_v);
				double basis0and1 = basis0 + half1_point1_length_bezier_value;
				double basis0and1_derivative = basis0_derivative + half1_point1_length_bezier_derivative;
				cross_section_surface_derivative_point0_[iU][iV][0] *= basis0and1;
				cross_section_surface_derivative_point0_[iU][iV][1] *= basis0and1;
				cross_section_surface_normal_derivative_point0_[iU][iV][0] *= basis0and1_derivative;
				cross_section_surface_normal_derivative_point0_[iU][iV][1] *= basis0and1_derivative;
				cross_section_surface_normal_derivative_point0_[iU][iV][0] = jacobian_normed_x * cross_section_surface_normal_derivative_point0_[iU][iV][0];
				cross_section_surface_normal_derivative_point0_[iU][iV][1] = jacobian_normed_x * cross_section_surface_normal_derivative_point0_[iU][iV][1];
			}

			// now cross_section_surface_normal_derivative_point0_ = D(x/||x||)/Dx * Dp'/Dc
			// rotate it to get the normal: (x,y)->(-y,x)
			for (int iRest1 = 0; iRest1 < cross_section_surface_normal_derivative_half1_rest_.size(); ++iRest1)
			{
				RotateTangent(cross_section_surface_normal_derivative_half1_rest_[iRest1][iU][iV][0]);
				RotateTangent(cross_section_surface_normal_derivative_half1_rest_[iRest1][iU][iV][1]);
			}
			for (int iRest2 = 0; iRest2 < cross_section_surface_normal_derivative_half2_rest_.size(); ++iRest2)
			{
				RotateTangent(cross_section_surface_normal_derivative_half2_rest_[iRest2][iU][iV][0]);
				RotateTangent(cross_section_surface_normal_derivative_half2_rest_[iRest2][iU][iV][1]);
			}
			RotateTangent(cross_section_surface_normal_derivative_half1_point1_length_[iU][iV]);
			RotateTangent(cross_section_surface_normal_derivative_tangent_angle_[iU][iV]);
			RotateTangent(cross_section_surface_normal_derivative_half2_point1_length_[iU][iV]);
			RotateTangent(cross_section_surface_normal_derivative_point0_[iU][iV][0]);
			RotateTangent(cross_section_surface_normal_derivative_point0_[iU][iV][1]);
		}
	}
}


void CrossSectionSurface::getValue(double _u, double _v, Eigen::Vector2d& _value) const
{
	if (_v > 1)
	{
		//construct Bezier curve
		BezierCurve half2;
		getBezierCurve_Half2(_u, half2);

		Eigen::VectorXd value(2);
		half2.getValue(_v - 1, value);
		_value(0) = value(0);
		_value(1) = value(1);
		//_value.setZero(2);
	}
	else
	{
		BezierCurve half1;
		getBezierCurve_Half1(_u, half1);

		Eigen::VectorXd value(2);
		half1.getValue(_v, value);
		_value(0) = value(0);
		_value(1) = value(1);
	}
}

void CrossSectionSurface::getValue(int _u, int _v, Eigen::Vector2d& _value) const
{
	assert(_u < cross_section_values_.size());
	assert(cross_section_values_.size() > 0);
	assert(_v < cross_section_values_[0].size());
	_value = cross_section_values_[_u][_v];
}

void CrossSectionSurface::getValue_Degree2(double _u, double _v, Eigen::Vector3d& _value) const
{
	assert(half1_rest_.size() == 1);
	assert(half2_rest_.size() == 1);

	if (_v > 1)
	{

		Eigen::VectorXd half2_point1_length;
		half2_point1_length_.getValue(_u, half2_point1_length);
		Eigen::VectorXd tangent;
		tangent_.getValue(_u, tangent);
		tangent.normalize();
		//cout << "check tangent.normalize() CrossSectionSurface::getValue_Degree2(double _u, double _v, Eigen::Vector3d& _value) const" << endl;
		Eigen::VectorXd point0;
		if (point0_.getSpaceDimension() == 0)
		{
			point0.setZero(2);
		}
		else
		{
			point0_.getValue(_u, point0);
		}

		Eigen::VectorXd half2_point1;
		half2_point1 = half2_point1_length[0] * tangent + point0;// note we will need to add point0 here
		Eigen::VectorXd half2_point2;
		half2_rest_[0].getValue(_u, half2_point2);

		BezierCurve half2;
		half2.setControlPoints(std::vector<Eigen::VectorXd>{ point0, half2_point1, half2_point2 });

		Eigen::VectorXd value(2);
		half2.getValue_Degree2(_v - 1, value);
		_value(1) = value(0);
		_value(2) = value(1);
		_value(0) = 0.0;
	}
	else
	{

		Eigen::VectorXd half1_point1_length;
		half1_point1_length_.getValue(_u, half1_point1_length);

		Eigen::VectorXd tangent;
		tangent_.getValue(_u, tangent);
		tangent.normalize();
		//cout << "check tangent.normalize() CrossSectionSurface::getValue_Degree2(double _u, double _v, Eigen::Vector3d& _value) const" << endl;

		Eigen::VectorXd half1_point2;
		half1_rest_[0].getValue(_u, half1_point2);
		Eigen::VectorXd point0;
		if (point0_.getSpaceDimension() == 0)
		{
			point0.setZero(2);
		}
		else
		{
			point0_.getValue(_u, point0);
		}
		Eigen::VectorXd half1_point1 = half1_point1_length[0] * (-tangent) + point0;


		BezierCurve half1;
		half1.setControlPoints(std::vector<Eigen::VectorXd>{half1_point2, half1_point1, point0});

		Eigen::VectorXd value(2);
		half1.getValue_Degree2(_v, value);
		_value(1) = value(0);
		_value(2) = value(1);
		_value(0) = 0.0;
	}
}

void CrossSectionSurface::getBezierCurve_Half1(double _u, BezierCurve& _half1) const
{
	// notice the order of the control points arranged here,
	// the origin is set as the last control point
	std::vector<Eigen::VectorXd> control_points(2 + half1_rest_.size());
	int num_control = control_points.size();

	Eigen::VectorXd point0;
	if (point0_.getSpaceDimension() == 0)
	{
		point0.setZero(2);
	}
	else
	{
		point0_.getValue(_u, point0);
	}
	control_points.back() = point0;

	// get second to last control point
	Eigen::VectorXd half1_point1_length;
	half1_point1_length_.getValue(_u, half1_point1_length);

	Eigen::VectorXd tangent_angle;
	tangent_angle_.getValue(_u, tangent_angle);
	Eigen::Vector2d tangent;
	CrossSectionHelpr convert_tangent;
	convert_tangent.TangentAngle2Vector(tangent_angle(0), tangent);
	control_points[num_control - 2] = half1_point1_length[0] * (-tangent) + point0;

	for (int iRest = 0; iRest < half1_rest_.size(); ++iRest)
	{
		half1_rest_[iRest].getValue(_u, control_points[num_control - 3 - iRest]);
	}

	_half1.setControlPoints(control_points);
}

void CrossSectionSurface::getBezierCurve_Half2(double _u, BezierCurve& _half2) const
{
	std::vector<Eigen::VectorXd> control_points(2 + half2_rest_.size());

	Eigen::VectorXd point0;
	if (point0_.getSpaceDimension() == 0)
	{
		point0.setZero(2);
	}
	else
	{
		point0_.getValue(_u, point0);
	}
	control_points[0] = point0;

	// get second control point
	Eigen::VectorXd half2_point1_length;
	half2_point1_length_.getValue(_u, half2_point1_length);
	//Eigen::VectorXd tangent;
	//tangent_.getValue(_u, tangent);
	//tangent.normalize();
	//control_points[1] = half2_point1_length[0] * tangent;

	Eigen::VectorXd tangent_angle;
	tangent_angle_.getValue(_u, tangent_angle);

	////cout << "u=" << _u << ", tangent_angle_=" << tangent_angle[0] << endl;

	Eigen::Vector2d tangent;
	CrossSectionHelpr convert_tangent;
	convert_tangent.TangentAngle2Vector(tangent_angle(0), tangent);
	control_points[1] = half2_point1_length[0] * tangent + point0;

	// get rest of the control point
	for (int iRest = 0; iRest < half2_rest_.size(); ++iRest)
	{
		half2_rest_[iRest].getValue(_u, control_points[iRest + 2]);
	}

	//construct Bezier curve
	_half2.setControlPoints(control_points);
}

void CrossSectionSurface::getControlParameters(double _u, std::vector<Eigen::Vector2d>& _half1_rest, double& _half1_point1_length,
	double& _tangent_angle, double& _half2_point1_length, std::vector<Eigen::Vector2d>& _half2_rest) const
{
	_half1_rest.resize(half1_rest_.size());
	for (int iRest1 = 0; iRest1 < half1_rest_.size(); ++iRest1)
	{
		Eigen::VectorXd point;
		half1_rest_[iRest1].getValue(_u, point);
		_half1_rest[iRest1](0) = point(0);
		_half1_rest[iRest1](1) = point(1);
	}

	Eigen::VectorXd half1_point1_length;
	half1_point1_length_.getValue(_u, half1_point1_length);
	_half1_point1_length = half1_point1_length(0);

	Eigen::VectorXd tangent_angle;
	tangent_angle_.getValue(_u, tangent_angle);
	_tangent_angle = tangent_angle(0);

	Eigen::VectorXd half2_point1_length;
	half2_point1_length_.getValue(_u, half2_point1_length);
	_half2_point1_length = half2_point1_length(0);

	_half2_rest.resize(half2_rest_.size());
	for (int iRest = 0; iRest < half2_rest_.size(); ++iRest)
	{
		Eigen::VectorXd point;
		half2_rest_[iRest].getValue(_u, point);
		_half2_rest[iRest](0) = point(0);
		_half2_rest[iRest](1) = point(1);
	}
}

void CrossSectionSurface::getPoint0(double _u, Eigen::Vector2d& _point0) const
{
	if (point0_.getSpaceDimension() == 0)
	{
		_point0.setZero(2);
	}
	else
	{
		Eigen::VectorXd point;
		point0_.getValue(_u, point);
		_point0(0) = point(0);
		_point0(1) = point(1);
	}
}

const std::array<double, 2>& CrossSectionSurface::getUDomain() const
{
	return domain_u_;
}

void CrossSectionSurface::FindBSplineBasisVector(Eigen::SparseMatrix<double>& _basis_vectors, const BSpline& _spline, const std::vector<double>& _u) const
{
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(_u.size() * (_spline.getDegree() + 1));
	for (int iU = 0; iU < _u.size(); ++iU)
	{
		// construct column id = iU
		for (int iBasis = 0; iBasis < _spline.getSpaceDimension(); ++iBasis)
		{
			double basis_value = 0.0;
			int basis_state = _spline.getPeriodicBasisValue(iBasis, _u[iU], basis_value);
			if (basis_state != 0) // has nonzero value
			{
				triplets.emplace_back(Eigen::Triplet<double>(iBasis, iU, basis_value));
			}
		}
	}

	_basis_vectors.resize(_spline.getSpaceDimension(), _u.size());
	_basis_vectors.setFromTriplets(triplets.begin(), triplets.end());


}

void CrossSectionSurface::FindBezierBasisValues(std::vector<std::vector<double>>& _values, int _order, const std::vector<double>& _t) const
{
	_values.resize(_order, std::vector<double>(_t.size()));
	BezierCurve basis(_order - 1);
	for (int iBasis = 0; iBasis < _order; ++iBasis)
	{
		for (int i = 0; i < _t.size(); ++i)
		{
			_values[iBasis][i] = basis.getBasisValue(iBasis, _t[i]);
		}
	}
}

void CrossSectionSurface::RotateTangent(Eigen::Vector2d& _tangent)
{
	Eigen::Vector2d normal(-_tangent.y(), _tangent.x());
	_tangent = normal;
}

void CrossSectionSurface::CombineDerivative_CrossSectionSurfaceAndNormal_ControlCurves_And_ControlCurve_KeyControlPoints(
	Eigen::Matrix2Xd& _combined, const Eigen::Vector2d& _surface2curve, const Eigen::VectorXd& _curve2keys) const
{
	_combined.resize(Eigen::NoChange, _curve2keys.rows());
	_combined.row(0) = _surface2curve(0) * _curve2keys.transpose();
	_combined.row(1) = _surface2curve(1) * _curve2keys.transpose();
}

Eigen::VectorXd CrossSectionSurface::FunctionNormedX(const Eigen::VectorXd& _x) const
{
	return _x / _x.norm();
}