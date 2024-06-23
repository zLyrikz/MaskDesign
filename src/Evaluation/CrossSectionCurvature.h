#pragma once
#include "../MaskInterface/CrossSection.h"

class CrossSectionCurvature
{
public:
	CrossSectionCurvature();
	~CrossSectionCurvature();
	
	//NOTE only suitable for degree 2 cross-section
	// formula follow the paper https://inria.hal.science/inria-00072572/PDF/RR-4064.pdf Interpolation with Curvature Constraints, Hafsa Deddi, Hazel Everett, Sylvain Lazard
	// require MoreInfoUpdated in CrossSection
	void ComputeMaxCurvatureCandidates(const CrossSection* _cross_section, 
		double& _half1_start_max_curvature, double& _half1_mid_max_curvature, double& _half1_end_max_curvature
		, double& _half2_start_max_curvature, double& _half2_mid_max_curvature, double& _half2_end_max_curvature);

	void ComputeMaxCurvature(const CrossSection* _cross_section,
		double& _half1, double& _half2);

	// each gradient is w.r.t. control points(7 parameters) in this order: l1,theta,l2,half1_rest,point0 (no half2 rest because of degree2)
	// note we already include point0 at last
	void ComputeMaxCurvatureCandidates_WithGradient(const CrossSection* _cross_section,
		double& _half1_start_max_curvature, double& _half1_mid_max_curvature, double& _half1_end_max_curvature
		, double& _half2_start_max_curvature, double& _half2_mid_max_curvature, double& _half2_end_max_curvature,
		Eigen::VectorXd& _gradient_half1_start_max_curvature, Eigen::VectorXd& _gradient_half1_mid_max_curvature, Eigen::VectorXd& _gradient_half1_end_max_curvature
		, Eigen::VectorXd& _gradient_half2_start_max_curvature, Eigen::VectorXd& _gradient_half2_mid_max_curvature, Eigen::VectorXd& _gradient_half2_end_max_curvature);
	

private:
	void SetCrossSection(const CrossSection* _cross_section);
	double ComputeControlPolygonArea(const Eigen::Vector2d& _point1, const Eigen::Vector2d& _point2) const;// point0 is origin

	void MidMaxCurvatureDerivative(Eigen::VectorXd& _gradient, double _area, double _pow2_area, const Eigen::VectorXd& _derivative_area,
		const Eigen::VectorXd& _derivative_pow3_length_point1_to_mid, double _pow3_length_point1_to_mid) const;
	void TwoEndsCurvatuveDerivative(Eigen::VectorXd& _gradient, double _curvatuve, const Eigen::VectorXd& _derivative_area,
		const Eigen::VectorXd& _derivative_pow3_vector_length, double _pow3_vector_length) const;

	// triange area = 0.5 * (_v1 x _v2)
	void AreaGradient(Eigen::VectorXd& _gradient, double _signed_area, 
		const Eigen::Matrix2Xd& _v1_graident, const Eigen::Matrix2Xd& _v2_graident,
		const Eigen::Vector2d& _v1, const Eigen::Vector2d& _v2) const;
	// D||x||^3/Dx
	void NormXPow3Gradient(Eigen::Vector2d& _gradient, const Eigen::Vector2d& _x, double _norm_x) const;

	// for debug derivative only
	double Half1StartMaxCurvatureFunction(const Eigen::VectorXd& _x) const;
	double Half1MidMaxCurvatureFunction(const Eigen::VectorXd& _x) const;
	double Half1EndMaxCurvatureFunction(const Eigen::VectorXd& _x) const;
	double Half2StartMaxCurvatureFunction(const Eigen::VectorXd& _x) const;
	double Half2MidMaxCurvatureFunction(const Eigen::VectorXd& _x) const;
	double Half2EndMaxCurvatureFunction(const Eigen::VectorXd& _x) const;
	void compareDifference(const Eigen::VectorXd& _g, const Eigen::VectorXd& _g_fd) const;
	double TriangleAreaHalf1(const Eigen::VectorXd& _x) const;
	double TriangleAreaHalf2(const Eigen::VectorXd& _x) const;

private:
	const CrossSection* cross_section_;
	// cross product of its two edges, possibly be negative
	double half1_triangle_area_;
	double half2_triangle_area_;
};

