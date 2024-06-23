#pragma once
#include <array>
//#include "./src/Alglib/interpolation.h"
#include "../CurveNSurface/BSpline.h"
#include "../CurveNSurface/BezierCurve.h"
#include "../MaskInterface/CrossSection.h"
#include <Eigen/Core>
#include <vector>

// each control point of the bezier should be 2D
// XOY local frame is the YOZ plane in the world coordinate
class CrossSectionSurface
{
public:
	CrossSectionSurface();
	~CrossSectionSurface();

	void setControlCurve(
		const std::vector<BSpline>& _half1_rest, const BSpline& _half1_point1_length,
		const BSpline& _tangent, 
		const BSpline& _half2_point1_length, const std::vector<BSpline>& _half2_rest);

	// the case when each cross section is a degree 2 Bezier curve (so 3 control points for each half )
	void setControlCurve_Degree2(
		const BSpline& _half1_rest, const BSpline& _half1_point1_length,
		const BSpline& _tangent,
		const BSpline& _half2_point1_length, const BSpline& _half2_rest);

	//void setControlCurveByInterpolation_Interpolate2DTangent(const std::vector<const CrossSection*>& _cross_sections, const std::vector<double>& _u, int _degree, const std::array<double, 2>& _domain_u = std::array<double, 2>{0.0, 1.0});
	// not set the boundary control curve
	void setControlCurveByInterpolation(const std::vector<const CrossSection*>& _cross_sections, const std::vector<double>& _u, int _degree, const std::array<double, 2>& _domain_u = std::array<double, 2>{0.0, 1.0});
	// set half2 end curve by interpolation
	// note, not changing the [half2_rest_] size, only replace the last element
	void setBoundaryControlCurveByInterpolation(const std::vector<Eigen::VectorXd>& _interpolate_points, const std::vector<double>& _u, int _degree = 3, const std::array<double, 2>& _domain_u = std::array<double, 2>{0.0, 1.0});
	void setHalf1BoundaryControlCurveByInterpolation(const std::vector<Eigen::VectorXd>& _interpolate_points, const std::vector<double>& _u, int _degree = 3, const std::array<double, 2>& _domain_u = std::array<double, 2>{0.0, 1.0});

	// note regard u,v as the same samples if size doesn't change!!!
	void ComputeSampledValues(const std::vector<double>& _u, const std::vector<double>& _v);

	// pre-request: call ComputeSampledValues()
	void UpdateDerivativeInfo();

	// pre-request: call UpdateDerivativeInfo()
	// note input u, v must correspond to cross-section values sampled u,v
	// note each row of the Eigen::Matrix2Xd is the derivative w.r.t. the parameters of key cross-sections, the order of the parameters is the same with key cross-section order
	void findDerivative_CrossSectionSurfaceAndNormal_KeyControlPoints(
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
		std::vector<std::vector<std::array<Eigen::Matrix2Xd, 2>>>& _normal_point0) const;
	const std::vector<Eigen::VectorXd>* getControlCurveDerivativeKeyControlPoints() const;

	// when v in [0,1] it's half1, when v in [1,2] it's half2
	// return 3d point, x=0
	void getValue(double _u, double _v, Eigen::Vector2d& _value) const;
	void getValue(int _u, int _v, Eigen::Vector2d& _value) const;
	// a special case, this is when each side is a degree 2 bezier curve
	void getValue_Degree2(double _u, double _v, Eigen::Vector3d& _value) const;
	void getBezierCurve_Half1(double _u, BezierCurve& _half1) const;
	void getBezierCurve_Half2(double _u, BezierCurve& _half2) const;
	void getControlParameters(double _u, std::vector<Eigen::Vector2d>& _half1_rest, double& _half1_point1_length,
		double& _tangent_angle,
		double& _half2_point1_length, std::vector<Eigen::Vector2d>& _half2_rest) const;
	void getPoint0(double _u, Eigen::Vector2d& _point0) const;

	const std::array<double, 2>& getUDomain() const;

	// save and read result
	//void saveToFile(const std::string& _file) const;// input full file path, but file name do not include suffix(.txt)

public:
	/*******************************************************************************************************/
	// for debug
	/*******************************************************************************************************/
	Mesh control2_;
	Mesh control3_;
	Mesh control22_;
	Mesh control23_;

private:
	// sparse matrix column number = _u.size
	// each column of the sparse matrix is a vector of basis at the given point u
	// for periodic bspline only
	void FindBSplineBasisVector(Eigen::SparseMatrix<double>& _basis_vectors, const BSpline& _spline, const std::vector<double>& _u) const;
	// _values[i][j] = B_i(t_j)
	void FindBezierBasisValues(std::vector<std::vector<double>>& _values, int _order, const std::vector<double>& _t) const;
	// (x,y)->(-y,x)
	void RotateTangent(Eigen::Vector2d& _tangent);
	void CombineDerivative_CrossSectionSurfaceAndNormal_ControlCurves_And_ControlCurve_KeyControlPoints(
		Eigen::Matrix2Xd& _combined, const Eigen::Vector2d& _surface2curve, const Eigen::VectorXd& _curve2keys) const;

	// at a given u_i of a control curve, take derivative w.r.t. its control points(which is the key control points of the key cross sections)
	// [N(u_i)]M^(-1) is the derivative for all control curves
	// due to the linear relation, this derivative doesn't depend on the key control points(only related to the interpolation technique)
	// pre-request: setControlCurveByInterpolation is called(to get the interpolation matrix)
	void findDerivative_ControlCurve_KeyControlPoints(const std::vector<double>& _u);
	// derivative of cross-section surface sample points w.r.t. control curves (at sampled u and v)
	// and derivative of cross-section surface normal sample points w.r.t. control curves (at sampled u and v)
	// the input include cross-section values at u,v; for computing the jacobian x/||x|| in the normal derivative
	void findDerivative_CrossSectionSurfaceAndNormal_ControlCurves(const std::vector<double>& _u, const std::vector<double>& _v, 
		const std::vector<std::vector<Eigen::Vector2d>>& _cross_section_derivatives);
	//void findDerivative_CrossSectionSurfaceNormal_ControlCurves(const std::vector<double>& _u, const std::vector<double>& _v);
	void UpdateSampledDerivatives();

	// for derivative test
	void TestJacobianNormedX();
	Eigen::VectorXd FunctionNormedX(const Eigen::VectorXd& _x) const;

private:
	// at each u of the curves below, these form control points of two bezier curve
	// NOTE now the connect point may not be 0, so point1 = point0 +- length*tangent
	std::vector<BSpline> half1_rest_;// set of 2d curves
	BSpline half1_point1_length_;	 // 1d curve
	BSpline tangent_;				 // note, this may not be updated accordingly! 2d curve on sphere (i.e. point norm doesn't matter), by default tangent point to half2 point1
	BSpline tangent_angle_;			 // 1d curve, represent the angle between tangent and unit x; angle range in [-pi, pi]
	BSpline half2_point1_length_;	 // 1d curve
	std::vector<BSpline> half2_rest_;// set of 2d curve
	BSpline point0_;// the point that connect two halfs, allowed to be nonzero now; if space dim==0, then regard it as 0 function


	std::array<double, 2> domain_u_ = std::array<double, 2>{0.0, 0.1};// domain of the control curves

	// values at sampled u,v points
	std::vector<std::vector<Eigen::Vector2d>> cross_section_values_;
	std::vector<std::vector<Eigen::Vector2d>> cross_section_derivatives_;// derivative w.r.t. v
	// u and v values in the sampled discrete cushion
	// will be updated in ComputeSampledValues
	std::vector<double> u_values_;
	std::vector<double> v_values_;
	bool is_derivative_updated_ = false;// update derivative when there is sample v or u parameters are changed

	// for computing deriatives

	std::vector<double> key_u_;
	Eigen::SparseMatrix<double> key_interpolation_matrix_;
	// i-th Eigen::VectorXd(column major) is the derivative when take u_i at control curve(i.e. [N(u_i)]M^(-1))
	std::vector<Eigen::VectorXd> control_curve_derivative_key_control_points_;

	//// store the bezier curve values for faster gradient evaluation
	//std::vector<double> bezier_basis1_values_at_half1_;
	//std::vector<double> bezier_basis1_values_at_half2_;

	// at sampled u,v(v major) point of the cross-section surface,
	// take derivative of the sampled point w.r.t. the control curve
	// cross_section_surface_derivative_half1_rest_[i][u][v][d] is for d-th coordinate of i-th point(i.e. half1_rest_[i][d]), at sampled u,v of cross_section_surface
	std::vector<std::vector<std::vector<std::array<Eigen::Vector2d, 2>>>> cross_section_surface_derivative_half1_rest_;
	std::vector<std::vector<Eigen::Vector2d>> cross_section_surface_derivative_half1_point1_length_;
	std::vector<std::vector<Eigen::Vector2d>> cross_section_surface_derivative_tangent_angle_;
	std::vector<std::vector<Eigen::Vector2d>> cross_section_surface_derivative_half2_point1_length_;
	std::vector<std::vector<std::vector<std::array<Eigen::Vector2d, 2>>>> cross_section_surface_derivative_half2_rest_;
	std::vector<std::vector<std::array<Eigen::Vector2d, 2>>> cross_section_surface_derivative_point0_;

	// normal is (-y,x), where (x,y) is the normlized tangent
	std::vector<std::vector<std::vector<std::array<Eigen::Vector2d, 2>>>> cross_section_surface_normal_derivative_half1_rest_;
	std::vector<std::vector<Eigen::Vector2d>> cross_section_surface_normal_derivative_half1_point1_length_;
	std::vector<std::vector<Eigen::Vector2d>> cross_section_surface_normal_derivative_tangent_angle_;
	std::vector<std::vector<Eigen::Vector2d>> cross_section_surface_normal_derivative_half2_point1_length_;
	std::vector<std::vector<std::vector<std::array<Eigen::Vector2d, 2>>>> cross_section_surface_normal_derivative_half2_rest_;
	std::vector<std::vector<std::array<Eigen::Vector2d, 2>>> cross_section_surface_normal_derivative_point0_;


};

