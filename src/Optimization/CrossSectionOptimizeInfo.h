#pragma once
#include "../MeshViewer/MeshDefinition.h"
#include "../MeshProcessing/AabbTree.h"
#include "../MaskInterface/CushionSurface.h"
#include "../Evaluation/LinearTetrahedralFemEvaluation.h"
#include <Discregrid/All>

// TODO prvide deep copy functions
// we have dynamic allocate pointer member

//	// input x size = dim_x_
	// the correspondence is x0,x1... = half1_point1_length, tangent_angle_, half2_point1_length_, half1_rest, half2_rest_(except end point), point0;
	// note if half1 is degree 3, then half1_rest = (half1_point2, half1 point3) (point2 is in front of point3)

// NOTE: 
// 1. curvature barrier and curveture terms only for quadratic half 1 and half 2 
// 2. height field and length variance term is only for quadratic half1
// 3. length variance is only for quadratic half1
//  max curvature is not the max of all three (DO NOT use it as a barrier term, only use it as a max curvature check)
class CrossSectionOptimizeInfo
{
public:
	// input initial cushion must have an initial key cross sections defined
	CrossSectionOptimizeInfo(const Mesh* _head, const AabbTree* _head_tree, const CushionSurface* _initial_cushion, int _sample_v, int _half1_degree, int _half2_degree);
	CrossSectionOptimizeInfo(/*FEM parameters*/double _E, double _nu, std::vector<int> _sdf_function_index, std::vector < const Discregrid::CubicLagrangeDiscreteGrid*> _real_head_sdf, std::vector < const tmd::TriangleMeshDistance*> _head_sdf, double _collision_weight,
		int _max_iteration, double _epsilon_gradient, double _epsilon_fucntion, double _epsilon_step,
		/*evaluation weight*/double _weight_average_pressure, double _weight_force_distribution, double _weight_area_distribution,
		double _wieght_tetrahedral_inversion, 
		/*sigmoid parameters (force threshold)*/double _scale, double _epsilon_force,
		/*tetrahedralization parameters*/int _num_layer, double _height,
		const CushionSurface* _initial_cushion, int _sample_v, int _half1_degree, int _half2_degree);

	~CrossSectionOptimizeInfo();

	//void SetHeadMesh(const Mesh* _head);
	//void SetHeadTree(const AabbTree* _head_tree);
	void SetInitialCushion(const CushionSurface* _cushion, int _sample_v);
	void SetRegularizationTermWeight(double _weight);
	void SetSharpAnglePenaltyTermWeight(double _weight);
	void SetRegularityPenaltyTermWeight(double _weight);
	void SetGlobalLengthVarianceWeight(double _weight);// variance of half1 endpoint norm
	void SetCurvatureBarrierWeight(double _weight);
	void SetHeightFieldWeight(double _weight_height_field);
	void SetWidthRegularizationWeight(double _weight_width);
	void SetShapeRegularizationWeight(double _weight_shape);
	void SetSizeRegularizationWeight(double _weight_size);
	void SetConvexityWeight(double _weight_convexity);

	// limit half1 end point a bit more to avoid cushion length change too much
	// _half1_end absolute difference with initial < [_half1_end]
	// point absolute difference with initial < [_point]
	// length in [min length, max length]
	// angle in interval [-_angle, _angle]
	void SetConstraint(double _half1_end, double _point, double _max_length, double _min_length, double _angle);
public:

	void GetKeyCrossSectionsFromX(std::vector<CrossSection>& _key_cross_sections, const double* _x) const;
	void GetXFromKeyCrossSections(const std::vector<CrossSection>& _key_cross_sections, double* _x) const;

	void GetCurrentX(double* _x) const;
	void GetXFromCushionSurface(double* _x, const CushionSurface& _cushion) const;// note the dimension should be the same with initial cushion
	void ChangeCushionFromX(const double* _x, CushionSurface& _cushion, int _sample_v = 120);//TODO: note we interpolate beizer point with degree 3 bspline by default!

	int GetDimX() const;
	int GetNumKeyCrossSection() const;
	int GetDimCrossSection() const;
	int GetHalf1Degree() const;
	int GetHalf2Degree() const;
	const double* GetInitialX() const;
	int GetSampleV() const;

	const LinearTetrahedralFemEvaluation* GetFemEvaluator() const;
	void GetCurrentCushion(CushionSurface& _cushion) const;
	void SaveSimulationResultXml(std::string _filepath, std::string _filename) const;

	void SaveCushionSurface(std::string _filepath, std::string _filename) const;

public:
	//double ObjectiveFunction(const double* _x);
	//bool IsFeasible(const double* _x) const;


	double ObjectiveFunction2_WithGradient(const double* _x, Eigen::VectorXd& _gradient, double _inexact_gradient_epsilon = 0.0);
	bool CurvatureBarrier(/*const Eigen::VectorXd& _x*/);// feasibility check
	bool HeightFieldBarrier(/*const Eigen::VectorXd& _x*/);
	bool IsFeasible2(const Eigen::VectorXd& _x);

	void DecreaseCurvatureUntilFeasible();
private:
	void Initialization(const CushionSurface* _initial_cushion, int _sample_v, int _half1_degree, int _half2_degree);
	void FindHeightFieldBarrierFunctionParameterX0();
	// other evaluation terms

	// some regularity term on the control polyogn
	double SharpAnglePenalty(const std::vector<CrossSection>& _key_cross_sections) const;
	double PenaltyOnExtraAnlge(double _extra_angle) const;// used in ConcavePenalty

	// length variance of control polygon
	double Regularity(const std::vector<CrossSection>& _key_cross_sections) const;
	// note, this is replaced by LengthVariance function
	double GlobalLengthVarianceRegularity(const std::vector<CrossSection>& _key_cross_sections) const;

	// this is previous version
	// need update before using
	double Regularization(const double* _x) const;

	double CurvatureBarrier(const std::vector<CrossSection>& _key_cross_sections, bool& _is_violated) const;

	// new versions
	double Convexity(Eigen::VectorXd& _gradient) const;
	// only for quadratic cross-section!!!!!!!
	double CurvatureBarrier_WithGradient(const std::vector<CrossSection>& _cross_sections, Eigen::VectorXd& _gradient) const;
	// note, we regard point2 as the end point of half1 here
	double HeightFieldBarrier_WithGradient(const std::vector<CrossSection>& _cross_sections, Eigen::VectorXd& _gradient) const;
	double BarrierFunction_WithDerivative(double _x, double& _derivative) const;
	double WidthRegularization(const std::vector<CrossSection>& _cross_sections, Eigen::VectorXd& _gradient) const;
	// note, we regard point2 as the end point of half1 here
	double LengthVariance(const std::vector<CrossSection>& _key_cross_sections, Eigen::VectorXd& _gradient) const;
	double ShapeRegularization(const double* _x, Eigen::VectorXd& _gradient) const;
	double SizeRegularization(Eigen::VectorXd& _gradient) const;

private:
	// ultility function used in DecreaseCurvatureUntilFeasible();
	void SetInvolvedKeyCrossSection(std::vector<int>& _involved_key_cross_sections, int _key, int _frame_state);
	void DecreaseHaf1Curvature(int _key, double _rate);// for degree 2 only
	void DecreaseHaf2Curvature(int _key, double _rate);// for degree 2 only
	double FindVectorAngle(const Eigen::Vector2d& _vector) const;
private:
	//const Mesh* head_;
	//const AabbTree* head_tree_;
	int sample_v_;
	const CushionSurface* initial_cushion_;
	double* initial_x_;
	CushionSurface cushion_;// for iteration used in the optimization
	std::vector<CrossSection> key_cross_sections_;
	int half1_degree_;// degree - 1 is the rest control point size
	int half2_degree_;// degree - 1 is the rest control point size

	bool evaluator_indicator_;
	std::vector<LinearTetrahedralFemEvaluation> fem_evaluator_;
	double regularization_term_weight_;
	double sharp_angle_penalty_weight_;
	double regularity_weight_;

	double curvature_barrier_weight_;
	double cushion_thickness_;

	// new constraints
	double weight_height_field_;
	double weight_width_;
	double global_length_variance_weight_;
	double weight_shape_;
	double weight_size_;
	double weight_convexity_;
	std::vector<double> height_field_barrier_parameter_for_term3_;

	// box constraint parameters
	double box_half1_end_radius_;
	double box_point_radius_;
	double max_length_;
	double box_angle_radius_;
	double min_length_;

	int num_key_cross_section_;
	const int dim_each_cross_section_;
	int dim_x_;// num_key_cross_section_ * dim_each_cross_section_

	int function_evaluation_times_;// record totoal function evaluation time


	// for temporary testing
public:
	//CushionSurface record_cushion_;
	int curvature_violate_counter_ = 0;
};

