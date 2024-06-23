#include "CrossSectionOptimizeInfo.h"
#include "../Utility/ReadRhino.h"
#include "../Evaluation/Variance.h"
#include "../Evaluation/CrossSectionCurvature.h"

#include <random>
#include <iostream>
using std::cout;
using std::endl;

CrossSectionOptimizeInfo::CrossSectionOptimizeInfo(const Mesh* _head, const AabbTree* _head_tree, const CushionSurface* _initial_cushion, int _sample_v, int _half1_degree, int _half2_degree):
	dim_each_cross_section_(2 * (_half1_degree - 1) + 2 * (_half2_degree - 2) + 3) // note half rest point number = degree - 1
{
	Initialization(_initial_cushion, _sample_v, _half1_degree, _half2_degree);
	evaluator_indicator_ = false;
}

CrossSectionOptimizeInfo::CrossSectionOptimizeInfo(
	double _E, double _nu, std::vector<int> _sdf_function_index, std::vector < const Discregrid::CubicLagrangeDiscreteGrid*> _real_head_sdf, std::vector < const tmd::TriangleMeshDistance*> _head_sdf, double _collision_weight, int _max_iteration, double _epsilon_gradient, double _epsilon_fucntion, double _epsilon_step,
	double _weight_average_pressure, double _weight_force_distribution, double _weight_area_distribution,
	double _wieght_tetrahedral_inversion,
	double _scale, double _epsilon_force, 
	int _num_layer, double _height, 
	const CushionSurface* _initial_cushion, int _sample_v, int _half1_degree, int _half2_degree):
	//dim_each_cross_section_(2 * (_half1_degree - 1) + 2 * (_half2_degree - 2) + 3) // note half rest point number = degree - 1
	dim_each_cross_section_(2 * (_half1_degree - 1) + 2 * (_half2_degree - 2) + 5) // add point0 now! note half rest point number = degree - 1
{
	Initialization(_initial_cushion, _sample_v, _half1_degree, _half2_degree);

	cushion_thickness_ = _height;

	if (half1_degree_ == 2 && half2_degree_ == 2)
	{
		if (!CurvatureBarrier(/*initial_x*/))
		{
			cout << "Curvature bad, resolve it...";
			DecreaseCurvatureUntilFeasible();
			cout << "resolved" << endl;

		}
	}

	evaluator_indicator_ = true;
	fem_evaluator_.resize(_sdf_function_index.size());
	for (int iExpression = 0; iExpression < _sdf_function_index.size(); ++iExpression)
	{
		fem_evaluator_[iExpression].SetCushion(&cushion_);
		fem_evaluator_[iExpression].SetFemParameters(_E, _nu, _head_sdf[iExpression], _collision_weight, _max_iteration, _epsilon_gradient, _epsilon_fucntion, _epsilon_step, _sdf_function_index[iExpression], _real_head_sdf[iExpression]);
		fem_evaluator_[iExpression].SetObjectiveFunctionWeight(_weight_average_pressure, _weight_force_distribution, _weight_area_distribution);
		fem_evaluator_[iExpression].SetTetrahdedralInversionWeight(_wieght_tetrahedral_inversion);
		fem_evaluator_[iExpression].SetSigmoidParameters(_scale, _epsilon_force);// scale, force epsilon
		fem_evaluator_[iExpression].TetrahedralizeCushion(_num_layer, _height);// num_layer, height
		fem_evaluator_[iExpression].FindBarrierFunctionX0();
	}
	if (half1_degree_ > 1)
	{
		FindHeightFieldBarrierFunctionParameterX0();
	}
}

CrossSectionOptimizeInfo::~CrossSectionOptimizeInfo()
{
	if (initial_x_ != nullptr)
	{
		delete[] initial_x_;
		initial_x_ = nullptr;
	}
}

//void CrossSectionOptimizeInfo::SetHeadMesh(const Mesh* _head)
//{
//	head_ = _head;
//}
//
//void CrossSectionOptimizeInfo::SetHeadTree(const AabbTree* _head_tree)
//{
//	head_tree_ = _head_tree;
//}

void CrossSectionOptimizeInfo::SetInitialCushion(const CushionSurface* _cushion, int _sample_v)
{
	sample_v_ = _sample_v;
	initial_cushion_ = _cushion;
	cushion_ = *initial_cushion_;

	num_key_cross_section_ = _cushion->getKeyCrossSectionsNumber();
	dim_x_ = num_key_cross_section_ * dim_each_cross_section_;

	//cout << "num_key_cross_section_=" << num_key_cross_section_ << ", dim_each_cross_section_=" << dim_each_cross_section_ << ", dim_x =" << dim_x_ << endl;

	_cushion->getKeyCrossSections(key_cross_sections_);

	if (initial_x_ != nullptr)
	{
		delete[] initial_x_;
		initial_x_ = nullptr;
	}
	initial_x_ = new double[dim_x_];
	GetCurrentX(initial_x_);
}

void CrossSectionOptimizeInfo::SetRegularizationTermWeight(double _weight)
{
	regularization_term_weight_ = _weight;
}

void CrossSectionOptimizeInfo::SetSharpAnglePenaltyTermWeight(double _weight)
{
	sharp_angle_penalty_weight_ = _weight;
}

void CrossSectionOptimizeInfo::SetRegularityPenaltyTermWeight(double _weight)
{
	regularity_weight_ = _weight;
}

void CrossSectionOptimizeInfo::SetGlobalLengthVarianceWeight(double _weight)
{
	global_length_variance_weight_ = _weight;
}

void CrossSectionOptimizeInfo::SetCurvatureBarrierWeight(double _weight)
{
	curvature_barrier_weight_ = _weight;
}

void CrossSectionOptimizeInfo::SetHeightFieldWeight(double _weight_height_field)
{
	weight_height_field_ = _weight_height_field;
}

void CrossSectionOptimizeInfo::SetWidthRegularizationWeight(double _weight_width)
{
	weight_width_ = _weight_width;
}

void CrossSectionOptimizeInfo::SetShapeRegularizationWeight(double _weight_shape)
{
	weight_shape_ = _weight_shape;
}

void CrossSectionOptimizeInfo::SetSizeRegularizationWeight(double _weight_size)
{
	weight_size_ = _weight_size;
}

void CrossSectionOptimizeInfo::SetConvexityWeight(double _weight_convexity)
{
	weight_convexity_ = _weight_convexity;
}

void CrossSectionOptimizeInfo::SetConstraint(double _half1_end, double _point, double _max_length, double _min_length, double _angle)
{
	box_half1_end_radius_ = _half1_end;
	box_point_radius_ = _point;
	max_length_ = _max_length;
	box_angle_radius_ = _angle;
	min_length_ = _min_length;
}

void CrossSectionOptimizeInfo::GetKeyCrossSectionsFromX(std::vector<CrossSection>& _key_cross_sections, const double* _x) const
{
	assert(_key_cross_sections.size() == num_key_cross_section_);

	for (int iSection = 0; iSection < num_key_cross_section_; ++iSection)
	{
		std::vector<Eigen::Vector2d> half1_rest(half1_degree_-  1);
		double half1_point1_length = 0.0;
		double tangent_angle = 0.0;
		double half2_point1_length = 0.0;
		std::vector<Eigen::Vector2d> half2_rest(half2_degree_ - 2);//except end point

		int base_index = iSection * dim_each_cross_section_;
		half1_point1_length = _x[base_index];
		tangent_angle = _x[base_index + 1];
		half2_point1_length = _x[base_index + 2];

		int start_index_half1 = base_index + 3;
		for (int iRest1 = 0; iRest1 < half1_rest.size(); ++iRest1)
		{
			half1_rest[iRest1].x() = _x[start_index_half1 + 2 * iRest1];
			half1_rest[iRest1].y() = _x[start_index_half1 + 2 * iRest1 + 1];
		}
		int start_index_half2 = start_index_half1 + 2 * half1_rest.size();
		for (int iRest2 = 0; iRest2 < half2_rest.size(); ++iRest2)
		{
			half2_rest[iRest2].x() = _x[start_index_half2 + 2 * iRest2];
			half2_rest[iRest2].y() = _x[start_index_half2 + 2 * iRest2 + 1];
		}

		// add point0 at the end
		int start_index_point0 = start_index_half2 + 2 * half2_rest.size();
		Eigen::Vector2d point0(_x[start_index_point0], _x[start_index_point0 + 1]);

		_key_cross_sections[iSection].setControlPoints_SameDegree(half1_rest, half1_point1_length, tangent_angle, half2_point1_length, half2_rest, point0);
	}
}

void CrossSectionOptimizeInfo::GetXFromKeyCrossSections(const std::vector<CrossSection>& _key_cross_sections, double* _x) const
{
	assert(_key_cross_sections.size() == num_key_cross_section_);

	for (int iSection = 0; iSection < num_key_cross_section_; ++iSection)
	{
		const std::vector<Eigen::Vector2d>* half1_rest = nullptr;
		const double* half1_point1_length = nullptr;
		const double* tangent_angle = nullptr;
		const double* half2_point1_length = nullptr;
		const std::vector<Eigen::Vector2d>* half2_rest = nullptr;

		_key_cross_sections[iSection].getControlPoints(half1_rest, half1_point1_length, tangent_angle, half2_point1_length, half2_rest);

		int base_index = iSection * dim_each_cross_section_;
		_x[base_index] = *half1_point1_length;
		_x[base_index + 1] = *tangent_angle;
		_x[base_index + 2] = *half2_point1_length;
		int start_index_half1 = base_index + 3;
		for (int iRest1 = 0; iRest1 < half1_degree_ - 1; ++iRest1)
		{
			_x[start_index_half1 + 2 * iRest1] = (*half1_rest)[iRest1].x();
			_x[start_index_half1 + 2 * iRest1 + 1] = (*half1_rest)[iRest1].y();
		}
		int start_index_half2 = start_index_half1 + 2 * (half1_degree_ - 1);
		for (int iRest2 = 0; iRest2 < half2_degree_ - 2; ++iRest2)
		{
			_x[start_index_half2 + 2 * iRest2] = (*half2_rest)[iRest2].x();
			_x[start_index_half2 + 2 * iRest2 + 1] = (*half2_rest)[iRest2].y();
		}

		// add point0
		const Eigen::Vector2d* point0 = _key_cross_sections[iSection].GetPoint0();
		int start_index_point0 = start_index_half2 + 2 * (half2_degree_ - 2);
		_x[start_index_point0] = point0->x();
		_x[start_index_point0 + 1] = point0->y();

	}

}

void CrossSectionOptimizeInfo::GetCurrentX(double* _x) const
{
	GetXFromCushionSurface(_x, cushion_);
}

void CrossSectionOptimizeInfo::GetXFromCushionSurface(double* _x, const CushionSurface& _cushion) const
{
	std::vector<CrossSection> key_cross_sections;
	_cushion.getKeyCrossSections(key_cross_sections);
	GetXFromKeyCrossSections(key_cross_sections, _x);
}

void CrossSectionOptimizeInfo::ChangeCushionFromX(const double* _x, CushionSurface& _cushion, int _sample_v)
{
	GetKeyCrossSectionsFromX(key_cross_sections_, _x);

	_cushion.changeKeyCrossSections(key_cross_sections_);

	_cushion.ConstructCrossSectionSurface(3);

	_cushion.computeDiscreteCrossSections_NewFrames(_sample_v);// TODO this can be speed up by avoid repeatedly constructing bezier curve at the same frame

	_cushion.UpdateDerivativeInfo();

}

int CrossSectionOptimizeInfo::GetDimX() const
{
	return dim_x_;
}

int CrossSectionOptimizeInfo::GetNumKeyCrossSection() const
{
	return num_key_cross_section_;
}

int CrossSectionOptimizeInfo::GetDimCrossSection() const
{
	return dim_each_cross_section_;
}

int CrossSectionOptimizeInfo::GetHalf1Degree() const
{
	return half1_degree_;
}

int CrossSectionOptimizeInfo::GetHalf2Degree() const
{
	return half2_degree_;
}

const double* CrossSectionOptimizeInfo::GetInitialX() const
{
	return initial_x_;
}

int CrossSectionOptimizeInfo::GetSampleV() const
{
	return sample_v_;
}

const LinearTetrahedralFemEvaluation* CrossSectionOptimizeInfo::GetFemEvaluator() const
{
	return &fem_evaluator_[0];
}

void CrossSectionOptimizeInfo::GetCurrentCushion(CushionSurface& _cushion) const
{
	_cushion = cushion_;
}

void CrossSectionOptimizeInfo::SaveSimulationResultXml(std::string _filepath, std::string _filename) const
{
	fem_evaluator_[0].SaveSimulationResultXml(_filepath, _filename);
}

void CrossSectionOptimizeInfo::SaveCushionSurface(std::string _filepath, std::string _filename) const
{
	cushion_.SaveResults(_filepath, _filename);
}

double CrossSectionOptimizeInfo::ObjectiveFunction2_WithGradient(const double* _x, Eigen::VectorXd& _gradient, double _inexact_gradient_epsilon)
{
	ChangeCushionFromX(_x, cushion_, sample_v_);

	//================================================================= update key cross-section info =============================================================
	for (auto& iCrossSection : key_cross_sections_)// here key corss sections should already be updated (in ChangeCushionFromX())
	{
		iCrossSection.UpdateMoreInfo();
	}
	//cout << "global_length_variance_=" << global_length_variance_ << "    " << std::flush;


	//============================================================= height field =============================================================||
	Eigen::VectorXd weighted_height_field_gradient;
	weighted_height_field_gradient.setZero(dim_x_);
	double weighted_height_field = 0.0;
	if (weight_height_field_ > 0.0)
	{
		Eigen::VectorXd height_field_gradient;

		double height_field = HeightFieldBarrier_WithGradient(key_cross_sections_, height_field_gradient);
		// add weight
		weighted_height_field_gradient = weight_height_field_ * height_field_gradient;
		weighted_height_field = weight_height_field_ * height_field;
		//cout << "weighted_height_field=" << weighted_height_field << " weighted_height_field_gradient norm=" << weighted_height_field_gradient.norm() << endl;
	}
	//============================================================= width_regularization =============================================================||
	Eigen::VectorXd width_regularization_gradient;
	double width_regularization = WidthRegularization(key_cross_sections_, width_regularization_gradient);
	// add weight
	Eigen::VectorXd weighted_width_regularization_gradient = weight_width_ * width_regularization_gradient;
	double weighted_width_regularization = weight_width_ * width_regularization;
	//cout << "weighted_width_regularization=" << weighted_width_regularization << " weighted_width_regularization_gradient norm=" << weighted_width_regularization_gradient.norm() << endl;

	//============================================================= length_variance =============================================================||
	Eigen::VectorXd weighted_length_variance_gradient;
	weighted_length_variance_gradient.setZero(dim_x_);
	double weighted_length_variance = 0.0;
	if (global_length_variance_weight_ > 0.0 && half1_degree_ > 1)
	{
		Eigen::VectorXd length_variance_gradient;
		double length_variance = LengthVariance(key_cross_sections_, length_variance_gradient);
		// add weight
		Eigen::VectorXd weighted_length_variance_gradient = global_length_variance_weight_ * length_variance_gradient;
		double weighted_length_variance = global_length_variance_weight_ * length_variance;
		//cout << "weighted_length_variance=" << weighted_length_variance << " weighted_length_variance_gradient norm=" << weighted_length_variance_gradient.norm() << endl;
	}
	//============================================================= shape regularization =============================================================||
	Eigen::VectorXd shape_regularization_gradient;
	double shape_regularization = ShapeRegularization(_x, shape_regularization_gradient);
	// add weight
	Eigen::VectorXd weighted_shape_regularization_gradient = weight_shape_ * shape_regularization_gradient;
	double weighted_shape_regularization = weight_shape_ * shape_regularization;
	//cout << "weighted_shape_regularization=" << weighted_shape_regularization << " weighted_shape_regularization_gradient norm=" << weighted_shape_regularization_gradient.norm() << endl;

	//============================================================= size regularization =============================================================||
	Eigen::VectorXd size_regularization_gradient;
	double size_regularization = SizeRegularization(size_regularization_gradient);
	// add weight
	Eigen::VectorXd weighted_size_regularization_gradient = weight_size_ * size_regularization_gradient;
	double weighted_size_regularization = weight_size_ * size_regularization;
	//cout << "weighted_size_regularization=" << weighted_size_regularization << " weighted_size_regularization_gradient norm=" << weighted_size_regularization_gradient.norm() << endl;

	//============================================================= convexity =============================================================||
	Eigen::VectorXd convexity_gradient;
	double convexity = Convexity(convexity_gradient);
	// add weight
	Eigen::VectorXd weighted_convexity_gradient = weight_convexity_ * convexity_gradient;
	double weighted_convexity = weight_convexity_ * convexity;
	//cout << "weighted_convexity=" << weighted_convexity << " weighted_convexity_gradient norm=" << weighted_convexity_gradient.norm() << endl;

	//=============================================================FEM metric=============================================================||
	double value = 0.0;
	Eigen::VectorXd gradient_fem;
	gradient_fem.setZero(dim_x_);
	std::vector<double> fem_values(fem_evaluator_.size());
	std::vector<Eigen::VectorXd> fem_gradients(fem_evaluator_.size());
	for (int iExpression = 0; iExpression < fem_evaluator_.size(); ++iExpression)
	{
		//cout << "Compute FEM metric for expresion " << iExpression << endl;
		fem_evaluator_[iExpression].ReTetrahedralizeCushion();

		fem_values[iExpression] = fem_evaluator_[iExpression].Evaluate_AccurateForce_WithGradient(fem_gradients[iExpression]);
	}
	for (int iExpression = 0; iExpression < fem_evaluator_.size(); ++iExpression)
	{
		value += fem_values[iExpression];
		gradient_fem += fem_gradients[iExpression];
	}
	value /= double(fem_evaluator_.size());
	gradient_fem /= double(fem_evaluator_.size());
	//cout << "value =" << value << endl;

	++function_evaluation_times_;
	if (function_evaluation_times_ % 100 == 0)
	{
		//cout << endl << "ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo" << endl <<
		//	function_evaluation_times_ << "iterations " << endl <<
		//	"ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo" << endl;
		cout << " " << function_evaluation_times_ << " objective function evaluations " << endl;
	}

	_gradient = gradient_fem + 
		weighted_height_field_gradient + 
		weighted_width_regularization_gradient + 
		weighted_length_variance_gradient + 
		weighted_shape_regularization_gradient +
		weighted_size_regularization_gradient;

	return value +
		weighted_height_field+
		weighted_width_regularization +
		weighted_length_variance +
		weighted_shape_regularization +
		weighted_size_regularization;

}


void CrossSectionOptimizeInfo::Initialization(const CushionSurface* _initial_cushion, int _sample_v, int _half1_degree, int _half2_degree)
{
	sample_v_ = _sample_v;
	half1_degree_ = _half1_degree;
	half2_degree_ = _half2_degree;
	if (half2_degree_ < 2)
	{
		cout << "[WARNING CrossSectionOptimizeInfo::CrossSectionOptimizeInfo] half2 degree < 2; doesn't support such case" << endl;
	}

	initial_x_ = nullptr;
	SetInitialCushion(_initial_cushion, sample_v_);

	//evaluator.setData(_head, _head_tree, &cushion_);
	regularization_term_weight_ = 0.;
	sharp_angle_penalty_weight_ = 0.0;
	regularity_weight_ = 0.;
	global_length_variance_weight_ = 0.0;
	curvature_barrier_weight_ = 0.0;
	weight_height_field_ = 0.0;
	weight_width_ = 0.0;
	weight_shape_ = 0.0;
	weight_size_ = 0.0;
	weight_convexity_ = 0.0;

	box_half1_end_radius_ = 0.;
	box_point_radius_ = 0.;
	max_length_ = 0.;
	box_angle_radius_ = 0.;
	min_length_ = 0.;

	function_evaluation_times_ = 0;

	cushion_thickness_ = 3.0;
}

void CrossSectionOptimizeInfo::FindHeightFieldBarrierFunctionParameterX0()
{
	for (auto& iCrossSection : key_cross_sections_)
	{
		iCrossSection.UpdateMoreInfo();
	}
	height_field_barrier_parameter_for_term3_.resize(key_cross_sections_.size());
	for (int iCrossSection = 0; iCrossSection < key_cross_sections_.size(); ++iCrossSection)
	{
		const CrossSection& cross_section = key_cross_sections_[iCrossSection];

		double theta = cross_section.getTangentAngle();

		// term: half1 point1.x > point2.x
		// B(point0.x - l*cos(theta) - point2.x)
		double length = cross_section.GetHalf1Vector01Length();
		double positive_value = cross_section.getConnectPointX() - length * cos(theta) - cross_section.getHalf1Point2X();
		if (positive_value < 0.0)
		{
			cout << "ERROR !!!!!!!!!!!CrossSectionOptimizeInfo::FindHeightFieldBarrierFunctionParameterX0 initial value doesn't meet height field constraint!!" << endl;
		}
		height_field_barrier_parameter_for_term3_[iCrossSection] = positive_value;
	}
}

double CrossSectionOptimizeInfo::Regularization(const double* _x) const
{
	double sum = 0.0;

	// L2 norm (with weight on elements)
	for (int iCross = 0; iCross < num_key_cross_section_; ++iCross)
	{		
	}
	return sum;
}

double CrossSectionOptimizeInfo::CurvatureBarrier(const std::vector<CrossSection>& _key_cross_sections, bool& _is_violated) const
{
	double sum = 0.0;
	double max_curvature_required = 1 / cushion_thickness_;
	int violation_counter = 0;
	for (auto& iCrossSection : _key_cross_sections)
	{
		CrossSectionCurvature compute_curvature;
		double half1_start_max_curvature = 0.0;
		double half1_mid_max_curvature = 0.0;
		double half1_end_max_curvature = 0.0;
		double half2_start_max_curvature = 0.0;
		double half2_mid_max_curvature = 0.0;
		double half2_end_max_curvature = 0.0;
		compute_curvature.ComputeMaxCurvatureCandidates(&iCrossSection,
			half1_start_max_curvature, half1_mid_max_curvature, half1_end_max_curvature,
			half2_start_max_curvature, half2_mid_max_curvature, half2_end_max_curvature);

		double temp1 = 1.0 / half1_start_max_curvature - cushion_thickness_;
		if (temp1 < 0.0)
		{
			violation_counter += 1;
		}
		else
		{			
			sum -= log(temp1);
		}
		double temp2 = 1.0 / half1_mid_max_curvature - cushion_thickness_;
		if (temp2 < 0.0)
		{
			violation_counter += 1;

		}
		else
		{
			sum -= log(temp2);
		}
		double temp3 = 1.0 / half1_end_max_curvature - cushion_thickness_;
		if (temp3 < 0.0)
		{
			violation_counter += 1;

		}
		else
		{
			sum -= log(temp3);
		}
		double temp4 = 1.0 / half2_start_max_curvature - cushion_thickness_;
		if (temp4 < 0.0)
		{
			violation_counter += 1;

		}
		else
		{
			sum -= log(temp4);
		}
		double temp5 = 1.0 / half2_mid_max_curvature - cushion_thickness_;
		if (temp5 < 0.0)
		{
			violation_counter += 1;

		}
		else
		{
			sum -= log(temp5);
		}
		double temp6 = 1.0 / half2_end_max_curvature - cushion_thickness_;
		if (temp6 < 0.0)
		{
			violation_counter += 1;

		}
		else
		{
			sum -= log(temp6);
		}

	}

	sum += double(violation_counter) * 1e10;
	if (violation_counter > 0)
	{
		_is_violated = true;
	}
	else
	{
		_is_violated = false;
	}
	return sum;
}

double CrossSectionOptimizeInfo::Convexity(Eigen::VectorXd& _gradient) const
{
	// add penalty at a control point when:
	// angle over 180  (concave penalty)

	double penalty_sum = 0.;

	_gradient.setZero(dim_x_);
	int iKey = 0;
	for (auto& iCrossSection : key_cross_sections_)
	{
		if (half1_degree_ > 1)
		{
			double cross_prduct1 = iCrossSection.CrossProductHalf1Vector10AndVector12();
			if (cross_prduct1 < 0.0)
			{
				penalty_sum += cross_prduct1 * cross_prduct1;

				// get gradient
				Eigen::Vector2d vector10 = -(*iCrossSection.GetHalf1Vector01());
				const Eigen::Vector2d& vector12 = *iCrossSection.GetHalf1Vector12();
				const Eigen::Vector2d& tangent = *iCrossSection.GetTangent();
				double length = iCrossSection.GetHalf1Vector01Length();

				int base_index = iKey * dim_each_cross_section_;
				// l1
				_gradient(base_index) = tangent.x() * vector12.y() + vector10.x() * tangent.y() - (tangent.y() * vector12.x() + vector10.y() * tangent.x());
				// theta
				_gradient(base_index + 1) = length * ( -tangent.y() * vector12.y() + length * tangent.x()* tangent.x() - (tangent.x() * vector12.x() - vector10.y() * tangent.y()));
				// point2.x
				int start_index_half1 = base_index + 3;
				_gradient(start_index_half1) = -vector10.y();
				//point2.y
				_gradient(start_index_half1 + 1) = vector10.x();
				// point0.x
				int point0_x = base_index + dim_each_cross_section_ - 2;
				_gradient(point0_x) = vector10.y();
				// point0.y
				int point0_y = base_index + dim_each_cross_section_ - 1;
				_gradient(point0_y) = -vector10.x();
			}
			//if (half1_degree_ > 2)
		}

		if (half2_degree_ > 1)
		{
			double cross_prduct2 = iCrossSection.CrossProductHalf2Vector12AndVector10();
			if (cross_prduct2 < 0.0)
			{
				penalty_sum += cross_prduct2 * cross_prduct2;

				// get gradient
				Eigen::Vector2d vector10 = -(*iCrossSection.GetHalf2Vector01());
				const Eigen::Vector2d& vector12 = *iCrossSection.GetHalf2Vector12();
				const Eigen::Vector2d& tangent = *iCrossSection.GetTangent();
				double length = iCrossSection.GetHalf2Vector01Length();

				int base_index = iKey * dim_each_cross_section_;
				// theta
				_gradient(base_index + 1) = length * (tangent.y() * vector10.y() + vector12.x() * (-tangent.x()) - (-tangent.x() * vector10.x() + vector12.y() * tangent.y()));
				// l2
				_gradient(base_index + 2) = -tangent.x() * vector10.y() + vector12.x() * (-tangent.y()) - (-tangent.y() * vector10.x() + vector12.y() * (-tangent.x()));
				// point2.x
				int start_index_half2 = base_index + 3 + 2 * (half1_degree_ - 1);
				_gradient(start_index_half2) = vector10.y();
				//point2.y
				_gradient(start_index_half2 + 1) = -vector10.x();
				// point0.x
				int point0_x = base_index + dim_each_cross_section_ - 2;
				_gradient(point0_x) = -vector10.y();
				// point0.y
				int point0_y = base_index + dim_each_cross_section_ - 1;
				_gradient(point0_y) = vector10.x();
			}
			//if (half2_degree_ > 2)
		}
		iKey++;
	}

	return penalty_sum;

}

bool CrossSectionOptimizeInfo::CurvatureBarrier(/*const Eigen::VectorXd& _x*/)
{
	// check feasibility of curvature
	//ChangeCushionFromX(_x.data(), cushion_, sample_v_);

	std::vector<CrossSection> all_cross_sections;
	cushion_.FindCrossSectionForAllFrames(all_cross_sections);
	for (auto& iCrossSection : all_cross_sections)
	{
		iCrossSection.UpdateMoreInfo();
	}
	assert(dim_each_cross_section_ == 7);
	double sum = 0.0;
	double max_curvature_required = 1 / cushion_thickness_;
	bool is_feasible = true;


	for (auto& iCrossSection : all_cross_sections)
	{
		CrossSectionCurvature compute_curvature;
		double half1_max_curvature = 0.0;
		double half2_max_curvature = 0.0;
		compute_curvature.ComputeMaxCurvature(&iCrossSection, half1_max_curvature, half2_max_curvature);

		double temp1 = max_curvature_required - half1_max_curvature;
		if (temp1 < 0.0)
		{
			is_feasible = false;
			break;
		}

		double temp2 = max_curvature_required - half2_max_curvature;
		if (temp2 < 0.0)
		{
			is_feasible = false;
			break;
		}
	}
	return is_feasible;
}

bool CrossSectionOptimizeInfo::HeightFieldBarrier(/*const Eigen::VectorXd& _x*/)
{
	// NOTE assume key_cross_sections_ already updated
	// ChangeCushionFromX(_x.data(), cushion_, sample_v_);
	// 
	for (auto& iCrossSection : key_cross_sections_)
	{
		iCrossSection.UpdateMoreInfo();
	}
	// check feasibility of height field
	for (int iCrossSection = 0; iCrossSection < key_cross_sections_.size(); ++iCrossSection)
	{
		const CrossSection& cross_section = key_cross_sections_[iCrossSection];

		double theta = cross_section.getTangentAngle();

		// term1:  theta > - M_PI_2
		if (theta <= -M_PI_2)
		{
			return false;

		}

		// term2:  theta < M_PI_2
		if (theta >= M_PI_2)
		{
			return false;
		}

		// term3: half1 point1.x > point2.x
		// B(point0.x - l*cos(theta) - point2.x)
		double length = cross_section.GetHalf1Vector01Length();
		if ((cross_section.getConnectPointX() - length * cos(theta) - cross_section.getHalf1Point2X()) <= 0.0)
		{
			return false;
		}

	}
	return true;
}

bool CrossSectionOptimizeInfo::IsFeasible2(const Eigen::VectorXd& _x)
{
	ChangeCushionFromX(_x.data(), cushion_, sample_v_);
	if (weight_height_field_ > 0.0)
	{
		if (!HeightFieldBarrier())
		{
			//cout << "height filed constraint violated" << endl;
			return false;
		}
	}

	if (curvature_barrier_weight_ > 0.0)
	{
		if (!CurvatureBarrier(/*_x*/))
		{
			cout << "curvature > 1/cushion_thickness_, constraint violated" << endl;
			return false;
		}
	}

	fem_evaluator_[0].ReTetrahedralizeCushion();
	if (!fem_evaluator_[0].IsInversionFree())
	{
		//cout << "tetrahdedral inversion constraint violated" << endl;
		return false;

	}	

	return true;
}

void CrossSectionOptimizeInfo::DecreaseCurvatureUntilFeasible()
{
	// get key cross-section and frame u values
	std::vector<double> key_cross_section_u;
	cushion_.getKeyCrossSectionsU(key_cross_section_u);
	std::vector<double> frames_u(cushion_.getDiscreteCushion()->size());
	for (int iFrame = 0; iFrame < cushion_.getDiscreteCushion()->size(); ++iFrame)
	{
		frames_u[iFrame] = (cushion_.getDiscreteCushion()->at(iFrame).u_);
	}
	int count2 = 0;

	do
	{
		// 1 find problematic frame v values
		std::vector<CrossSection> all_cross_sections;
		cushion_.FindCrossSectionForAllFrames(all_cross_sections);
		for (auto& iCrossSection : all_cross_sections)
		{
			iCrossSection.UpdateMoreInfo();
		}
		assert(dim_each_cross_section_ == 7);
		double max_curvature_required = 1 / cushion_thickness_;
		// pair (i,1) => frame i half1 curveture bad
		// pair (i,2) => frame i half2 curveture bad
		// pair (i,3) => frame i half1 and 2 both curveture bad
		std::vector<std::pair<int, int>> problem_frame;

		int count = 0;
		for (auto& iCrossSection : all_cross_sections)
		{
			CrossSectionCurvature compute_curvature;
			double half1_max_curvature = 0.0;
			double half2_max_curvature = 0.0;
			compute_curvature.ComputeMaxCurvature(&iCrossSection, half1_max_curvature, half2_max_curvature);
			bool is_half1_curvature_bad = false;
			bool is_half2_curvature_bad = false;
			double temp1 = max_curvature_required - half1_max_curvature;
			if (temp1 < 0.0)
			{
				is_half1_curvature_bad = true;
				//cout << "temp1=" << temp1 << endl;
			}

			double temp2 = max_curvature_required - half2_max_curvature;
			if (temp2 < 0.0)
			{
				is_half2_curvature_bad = true;
				//cout << "temp2=" << temp2 << endl;

			}

			if (is_half1_curvature_bad && is_half2_curvature_bad)
			{
				problem_frame.push_back(std::make_pair(count, 3));
				//cout << "frame=" << count<<" u="<< frames_u[count] << ", is_half1_curvature_bad && is_half2_curvature_bad" << endl;
			}
			else if (is_half1_curvature_bad && !is_half2_curvature_bad)
			{
				problem_frame.push_back(std::make_pair(count, 1));
				//cout << "frame=" << count << " u=" << frames_u[count] << ", is_half1_curvature_bad" << endl;

			}
			else if (!is_half1_curvature_bad && is_half2_curvature_bad)
			{
				problem_frame.push_back(std::make_pair(count, 2));
				//cout << "frame=" << count << " u=" << frames_u[count] << ", is_half2_curvature_bad" << endl;

			}
			count++;
		}
	
		// 2 find the involved key cross sections
		std::vector<int> involved_key_cross_sections(num_key_cross_section_, 0);
		for (int iFrame = 0; iFrame < problem_frame.size(); ++iFrame)
		{
			if (frames_u[problem_frame[iFrame].first] > key_cross_section_u.back())
			{
				// this problem frame beteween (last, first]
				SetInvolvedKeyCrossSection(involved_key_cross_sections, 0, problem_frame[iFrame].second);
				SetInvolvedKeyCrossSection(involved_key_cross_sections, num_key_cross_section_ - 1, problem_frame[iFrame].second);
				continue;
			}
			for (int iKey = 0; iKey < num_key_cross_section_; ++iKey)
			{
				if (key_cross_section_u[iKey] > frames_u[problem_frame[iFrame].first])
				{
					// this problem frame beteween (iKey-1, iKey]
					int previous_key = (iKey - 1 == -1) ? (num_key_cross_section_ - 1) : (iKey - 1);
					if (previous_key != num_key_cross_section_ - 1)
					{
						if (frames_u[problem_frame[iFrame].first] <= key_cross_section_u[previous_key])
						{
							cout << "ERROR CrossSectionOptimizeInfo::DecreaseCurvatureUntilFeasible unexpected" << endl;
						}
					}
					SetInvolvedKeyCrossSection(involved_key_cross_sections, iKey, problem_frame[iFrame].second);
					SetInvolvedKeyCrossSection(involved_key_cross_sections, previous_key, problem_frame[iFrame].second);
					break;
				}
				else if (key_cross_section_u[iKey] == frames_u[problem_frame[iFrame].first])
				{
					// this problem frame on iKey
					SetInvolvedKeyCrossSection(involved_key_cross_sections, iKey, problem_frame[iFrame].second);
					break;
				}
			}
		}

		double rate = 0.1;

		for (int iKey = 0; iKey < num_key_cross_section_; ++iKey)
		{
			if (involved_key_cross_sections[iKey] == 1 || involved_key_cross_sections[iKey] == 3)
			{
				DecreaseHaf1Curvature(iKey, rate);
				//cout << "DecreaseHaf1Curvature " << iKey << " u="<< key_cross_section_u[iKey] << endl;
			}
			if (involved_key_cross_sections[iKey] == 2 || involved_key_cross_sections[iKey] == 3)
			{
				DecreaseHaf2Curvature(iKey, rate);
				//cout << "DecreaseHaf2Curvature " << iKey << " u=" << key_cross_section_u[iKey] << endl;

			}
		}

		for (auto& iCrossSection : key_cross_sections_)
		{
			iCrossSection.UpdateMoreInfo();
		}
		GetXFromKeyCrossSections(key_cross_sections_, initial_x_);
		//for (int iX = 0; iX < dim_x_; ++iX)
		//{
		//	x[iX] = initial_x_[iX];
		//}
		count2++;

		ChangeCushionFromX(initial_x_, cushion_, sample_v_);


	} while (!CurvatureBarrier(/*x*/) /*&& count2 < 100*/);
}

double CrossSectionOptimizeInfo::CurvatureBarrier_WithGradient(const std::vector<CrossSection>& _cross_sections, Eigen::VectorXd& _gradient) const
{
	//assert(dim_each_cross_section_ == 5);
	double sum = 0.0;
	double max_curvature_required = 1 / cushion_thickness_;
	int violation_counter = 0;
	int current_cross_section_id = 0;
	_gradient.setZero(dim_x_);
	for (auto& iCrossSection : _cross_sections)
	{
		CrossSectionCurvature compute_curvature;
		double half1_start_max_curvature = 0.0;
		double half1_mid_max_curvature = 0.0;
		double half1_end_max_curvature = 0.0;
		double half2_start_max_curvature = 0.0;
		double half2_mid_max_curvature = 0.0;
		double half2_end_max_curvature = 0.0;
		Eigen::VectorXd gradient_half1_start_max_curvature; Eigen::VectorXd gradient_half1_mid_max_curvature; Eigen::VectorXd gradient_half1_end_max_curvature;
		Eigen::VectorXd gradient_half2_start_max_curvature; Eigen::VectorXd gradient_half2_mid_max_curvature; Eigen::VectorXd gradient_half2_end_max_curvature;
		compute_curvature.ComputeMaxCurvatureCandidates_WithGradient(&iCrossSection,
			half1_start_max_curvature, half1_mid_max_curvature, half1_end_max_curvature,
			half2_start_max_curvature, half2_mid_max_curvature, half2_end_max_curvature, 
			gradient_half1_start_max_curvature, gradient_half1_mid_max_curvature, gradient_half1_end_max_curvature,
			gradient_half2_start_max_curvature, gradient_half2_mid_max_curvature, gradient_half2_end_max_curvature);

		double temp1 = max_curvature_required - half1_start_max_curvature;
		if (temp1 < 0.0)
		{
			violation_counter += 1;
		}
		else
		{
			sum -= log(temp1);
			gradient_half1_start_max_curvature /= temp1;
		}
		double temp2 = max_curvature_required - half1_mid_max_curvature;
		if (temp2 < 0.0)
		{
			violation_counter += 1;
		}
		else
		{
			sum -= log(temp2);
			gradient_half1_mid_max_curvature /= temp2;
		}
		double temp3 = max_curvature_required - half1_end_max_curvature ;
		if (temp3 < 0.0)
		{
			violation_counter += 1;
		}
		else
		{
			sum -= log(temp3);
			gradient_half1_end_max_curvature /= temp3;
		}
		double temp4 = max_curvature_required - half2_start_max_curvature ;
		if (temp4 < 0.0)
		{
			violation_counter += 1;
		}
		else
		{
			sum -= log(temp4);
			gradient_half2_start_max_curvature /= temp4;
		}
		double temp5 = max_curvature_required - half2_mid_max_curvature ;
		if (temp5 < 0.0)
		{
			violation_counter += 1;
		}
		else
		{
			sum -= log(temp5);
			gradient_half2_mid_max_curvature /= temp5;
		}
		double temp6 = max_curvature_required - half2_end_max_curvature ;
		if (temp6 < 0.0)
		{
			violation_counter += 1;
		}
		else
		{
			sum -= log(temp6);
			gradient_half2_end_max_curvature /= temp6;
		}

		// arrange gredients
		Eigen::VectorXd curvature_barrier_gradient_sum =
			gradient_half1_start_max_curvature + gradient_half1_mid_max_curvature + gradient_half1_end_max_curvature +
			gradient_half2_start_max_curvature + gradient_half2_mid_max_curvature + gradient_half2_end_max_curvature;
		const CrossSectionSurface& cross_section_surface = cushion_.getCrossSectionSurface();
		const std::vector<Eigen::VectorXd>& control_curve_derivative_key_control_points = *(cross_section_surface.getControlCurveDerivativeKeyControlPoints());
		// each row is one element(e.g. theta) gradient w.r.t. key elements(e.g. key theta)
		// size: (row,col)= (7, num_key_cross_section_)
		Eigen::MatrixXd gradient_current_section;
		gradient_current_section = curvature_barrier_gradient_sum * (control_curve_derivative_key_control_points[current_cross_section_id]).transpose();
		// 
		// put them in the final gradient
		//already include point0
		for (int iKey = 0; iKey < num_key_cross_section_; ++iKey)
		{
			int base_index = iKey * dim_each_cross_section_;

			_gradient.block(base_index, 0, dim_each_cross_section_, 1) += gradient_current_section.block(0, iKey, dim_each_cross_section_, 1);
		}

		current_cross_section_id += 1;
	}

	//sum += double(violation_counter) * 1e5;
	if (violation_counter > 0)
	{
		cout << "ERROR!!!!!!!!!!CrossSectionOptimizeInfo::CurvatureBarrier_WithGradient curvature of radius < cushion_thickness_; violation_counter="<< violation_counter << endl;
	}

	return sum;

}

double CrossSectionOptimizeInfo::HeightFieldBarrier_WithGradient(const std::vector<CrossSection>& _cross_sections, Eigen::VectorXd& _gradient) const
{
	double sum = 0.0;
	_gradient.setZero(dim_x_);
	// note, the cross-sections here are considered as key cross-sections!!! (we don't need the interpolation matrix for derivative)
	if (_cross_sections.size() != num_key_cross_section_)
	{
		cout << "ERROR!!!!!CrossSectionOptimizeInfo::HeightFieldBarrier_WithGradient _cross_sections.size() != num_key_cross_section_" << endl;
	}
	for (int iCrossSection = 0; iCrossSection < _cross_sections.size(); ++iCrossSection)
	{
		const CrossSection& cross_section = _cross_sections[iCrossSection];

		double theta = cross_section.getTangentAngle();

		double theta_term1_derivative = 0.0;
		//sum += BarrierFunction_WithDerivative(theta + M_PI_2, theta_term1_derivative);
		// NOTE we set the x0 for this barrier function as M_PI_2 / 2.0 here
		sum += fem_evaluator_[0].IpcBarrierFunction_WithDerivative(theta + M_PI_2, M_PI_2 / 2.0, theta_term1_derivative);

		double theta_term2_derivative = 0.0;
		//sum += BarrierFunction_WithDerivative(M_PI_2 - theta, theta_term2_derivative);
		sum += fem_evaluator_[0].IpcBarrierFunction_WithDerivative(M_PI_2 - theta, M_PI_2 / 2.0, theta_term2_derivative);

		// term: half1 point1.x > point2.x
		// B(point0.x - l*cos(theta) - point2.x)
		double length = cross_section.GetHalf1Vector01Length();
		double term3_derivative = 0.0;
		//sum += BarrierFunction_WithDerivative(cross_section.getConnectPointX() - length * cos(theta) - cross_section.getHalf1Point2X(), term3_derivative);
		sum += fem_evaluator_[0].IpcBarrierFunction_WithDerivative(
			cross_section.getConnectPointX() - length * cos(theta) - cross_section.getHalf1Point2X(), 
			height_field_barrier_parameter_for_term3_[iCrossSection],
			term3_derivative);

		// derivative
		int base_index = iCrossSection * dim_each_cross_section_;
		// term1
		_gradient(base_index + 1) += theta_term1_derivative;
		// term2
		_gradient(base_index + 1) += (-theta_term2_derivative);
		// term3: point0.x - l*cos(theta) - point2.x
		// 
		// point0.x
		int start_index_point0 = base_index + dim_each_cross_section_ - 2;
		_gradient(start_index_point0) += term3_derivative;
		// theta
		_gradient(base_index + 1) -= length * (-sin(theta))* term3_derivative;
		// length half1
		_gradient(base_index) -= cos(theta)* term3_derivative;
		// point2.x
		int start_index_half1 = base_index + 3;
		_gradient(start_index_half1) -= term3_derivative;



	}

	return sum;
}

double CrossSectionOptimizeInfo::BarrierFunction_WithDerivative(double _x, double& _derivative) const
{
	double value = 1e10;
	if (_x > 0.0)
	{
		value = -log(_x);
		_derivative = -1.0 / _x;
	}
	else
	{
		cout << "WARNING!!!!!!!!!!! CrossSectionOptimizeInfo::BarrierFunction_WithDerivative x<=0, exceed barrier" << endl;
	}
	return value;
}

double CrossSectionOptimizeInfo::WidthRegularization(const std::vector<CrossSection>& _cross_sections, Eigen::VectorXd& _gradient) const
{
	// (sum point0.y)^2
	double value = 0.0;
	double sum = 0.0;
	if (_cross_sections.size() != num_key_cross_section_)
	{
		cout << "ERROR!!!!!!!!!! CrossSectionOptimizeInfo::WidthRegularization _cross_sections.size() != num_key_cross_section_" << endl;
	}
	for (int iCrossSection = 0; iCrossSection < num_key_cross_section_; ++iCrossSection)
	{
		const CrossSection& cross_section = _cross_sections[iCrossSection];
		sum += cross_section.GetPoint0()->y();
	}
	value = sum* sum;

	//gradient 
	// each partial gradient on point0.y is 2*(sum point0.y)
	_gradient.setZero(dim_x_);
	double partial_gradient = 2.0 * sum;
	for (int iCrossSection = 0; iCrossSection < num_key_cross_section_; ++iCrossSection)
	{
		_gradient(iCrossSection * dim_each_cross_section_ + dim_each_cross_section_ - 1) = partial_gradient;
	}


	return value;
}

double CrossSectionOptimizeInfo::LengthVariance(const std::vector<CrossSection>& _key_cross_sections, Eigen::VectorXd& _gradient) const
{
	if (_key_cross_sections.size() != num_key_cross_section_)
	{
		cout << "ERROR!!!!!!!!!! CrossSectionOptimizeInfo::LengthVariance _key_cross_sections.size() != num_key_cross_section_" << endl;
	}

	std::vector<double> half1_end_x(num_key_cross_section_);
	for (int iCrossSection = 0; iCrossSection < num_key_cross_section_; ++iCrossSection)
	{
		const CrossSection& cross_section = _key_cross_sections[iCrossSection];
		half1_end_x[iCrossSection] = cross_section.getHalf1Point2X();
	}

	Variance compute_variance;
	Eigen::VectorXd local_gradient;
	double length_variance = compute_variance.computeVariance_WithGradient(half1_end_x, local_gradient);

	_gradient.setZero(dim_x_);
	for (int iCrossSection = 0; iCrossSection < num_key_cross_section_; ++iCrossSection)
	{
		int base_index = iCrossSection * dim_each_cross_section_;
		int start_index_half1 = base_index + 3;
		_gradient(start_index_half1) = local_gradient[iCrossSection];
	}

	return length_variance;
}

double CrossSectionOptimizeInfo::ShapeRegularization(const double* _x, Eigen::VectorXd& _gradient) const
{
	std::vector<std::vector<double>> difference(dim_each_cross_section_, std::vector<double>(num_key_cross_section_));// difference[j][i] is for j-th parameter on i-th key cross-section. value = x_i^j - x_{i+1}^j
	for (int iParameter = 0; iParameter < dim_each_cross_section_; ++iParameter)
	{
		for (int iKey = 0; iKey < num_key_cross_section_; ++iKey)
		{
			int variable_id = iKey * dim_each_cross_section_ + iParameter;
			int variable_id_next_key = ((iKey + 1) % num_key_cross_section_) * dim_each_cross_section_ + iParameter;
			difference[iParameter][iKey] = _x[variable_id] - _x[variable_id_next_key];
		}
	}
	double shape_regularization = 0.0;
	for (int iParameter = 0; iParameter < dim_each_cross_section_; ++iParameter)
	{
		for (int iKey = 0; iKey < num_key_cross_section_; ++iKey)
		{
			shape_regularization += difference[iParameter][iKey] * difference[iParameter][iKey];
		}
	}
	// gradient
	_gradient.setZero(dim_x_);
	for (int iParameter = 0; iParameter < dim_each_cross_section_; ++iParameter)
	{
		for (int iKey = 0; iKey < num_key_cross_section_; ++iKey)
		{
			int variable_id = iKey * dim_each_cross_section_ + iParameter;
			int previous_key = (iKey - 1) == -1 ? num_key_cross_section_ - 1 : (iKey - 1);
			_gradient(variable_id) = 2.0 * (difference[iParameter][iKey] - difference[iParameter][previous_key]);
		}
	}
	return shape_regularization;
}

double CrossSectionOptimizeInfo::SizeRegularization(Eigen::VectorXd& _gradient) const
{
	double value = 0.0;
	for (int iKey = 0; iKey < key_cross_sections_.size(); ++iKey)
	{
		value += key_cross_sections_[iKey].getConnectPointX() * key_cross_sections_[iKey].getConnectPointX();
	}
	// gradient
	_gradient.setZero(dim_x_);
	for (int iKey = 0; iKey < key_cross_sections_.size(); ++iKey)
	{
		int point0_x_id = iKey * dim_each_cross_section_ + dim_each_cross_section_ - 2;
		_gradient(point0_x_id) = 2.0 * key_cross_sections_[iKey].getConnectPointX();
	}
	return value;
}

void CrossSectionOptimizeInfo::SetInvolvedKeyCrossSection(std::vector<int>& _involved_key_cross_sections, int _key, int _frame_state)
{
	if (_frame_state == 3)
	{
		_involved_key_cross_sections[_key] = 3;
	}
	else if (_frame_state == 2)
	{
		if (_involved_key_cross_sections[_key] == 1)
		{
			_involved_key_cross_sections[_key] = 3;
		}
		else if (_involved_key_cross_sections[_key] == 0)
		{
			_involved_key_cross_sections[_key] = 2;
		}
	}
	else if (_frame_state == 1)
	{
		if (_involved_key_cross_sections[_key] == 2)
		{
			_involved_key_cross_sections[_key] = 3;
		}
		else if (_involved_key_cross_sections[_key] == 0)
		{
			_involved_key_cross_sections[_key] = 1;
		}
	}
}

void CrossSectionOptimizeInfo::DecreaseHaf1Curvature(int _key, double _rate)
{
	key_cross_sections_[_key].UpdateMoreInfo();
	Eigen::Vector2d vector10 = -(*key_cross_sections_[_key].GetHalf1Vector01());
	Eigen::Vector2d vector12 = *key_cross_sections_[_key].GetHalf1Vector12();
	Eigen::Vector2d point1 = *key_cross_sections_[_key].GetHalf1Point1();
	Eigen::Vector2d point2_new;// to be found
	if (!key_cross_sections_[_key].IsConvex_Half1Point1())
	{
		// just straighten vector12
		point2_new = point1 - vector12.norm() * vector10.normalized();
	}
	else
	{
		// rotate vector12 toward paralell with vector01
		//cout << "key_cross_sections_" << _key;
		double angle = acos(vector10.normalized().dot(vector12.normalized()));
		//cout << " angle "<<angle << " ";
		double supplement_angle = M_PI - angle;
		double angle_moved = _rate * supplement_angle + angle + FindVectorAngle(vector10);
		//cout << "angle_moved " << angle_moved;

		Eigen::Vector2d vector12_new = vector12.norm() * Eigen::Vector2d(cos(angle_moved), sin(angle_moved));
		//cout << " vector12 " << vector12.transpose();
		//cout << " vector12_new " << vector12_new.transpose() << endl;
		point2_new = vector12_new + point1;
	}
	key_cross_sections_[_key].setHalf1EndPoint(point2_new);
}

void CrossSectionOptimizeInfo::DecreaseHaf2Curvature(int _key, double _rate)
{
	// move half2_point1 to the midpoint of point0 and half2_point2
	key_cross_sections_[_key].UpdateMoreInfo();

	// 1. find midpoint
	const Eigen::Vector2d* point2 = key_cross_sections_[_key].GetHalf2Point2();
	const Eigen::Vector2d* point0 = key_cross_sections_[_key].GetPoint0();
	Eigen::Vector2d midpoint = (*point2 + *point0) / 2.0;

	// 2. move
	if (_rate < 0.0 || _rate >1.0)
	{
		cout << "ERROR!!!!!!!!!!!!!! CrossSectionOptimizeInfo::DecreaseHaf2Curvature move _rate < 0.0 || _rate >1.0" << endl;
	}
	Eigen::Vector2d moved_point1 = (1.0 - _rate) * *key_cross_sections_[_key].GetHalf2Point1() + _rate * midpoint;

	// 3. find the changed tangent angle
	Eigen::Vector2d new_vector01 = moved_point1 - *point0;
	Eigen::Vector2d old_vector01 = *key_cross_sections_[_key].GetHalf2Vector01();
	double changed_angle = acos(new_vector01.normalized().dot(old_vector01.normalized()));

	// 4. determine + or - this changed angle
	Eigen::Vector2d vector02 = (*key_cross_sections_[_key].GetHalf2Point2() - *key_cross_sections_[_key].GetPoint0());
	Eigen::Matrix4d cross;
	cross.block(0, 0, 2, 1) = old_vector01;
	cross.block(0, 1, 2, 1) = vector02;
	double new_tangent_angle = key_cross_sections_[_key].getTangentAngle();
	if (cross.determinant() < 0.0)
	{
		// -
		new_tangent_angle -= changed_angle;
	}
	else
	{
		// +
		new_tangent_angle += changed_angle;
	}

	// 5. find new half2_point1_length
	double new_half2_point1_length = new_vector01.norm();

	// find the new half1 point2
	Eigen::Vector2d half1_point1_new = -key_cross_sections_[_key].GetHalf1Vector01Length() * Eigen::Vector2d(cos(new_tangent_angle), sin(new_tangent_angle))
		+ *key_cross_sections_[_key].GetPoint0();
	Eigen::Vector2d half1_point2_new = half1_point1_new + *key_cross_sections_[_key].GetHalf1Vector12();
	key_cross_sections_[_key].setHalf1EndPoint(half1_point2_new);
	key_cross_sections_[_key].setTangentAngle(new_tangent_angle);
	key_cross_sections_[_key].setHalf2Point1Length(new_half2_point1_length);
}

double CrossSectionOptimizeInfo::FindVectorAngle(const Eigen::Vector2d& _vector) const
{
	double angle = acos(_vector.normalized()[0]);
	if (_vector[1] < 0.0)
	{
		angle = - angle;
	}
	return angle;
}

double CrossSectionOptimizeInfo::SharpAnglePenalty(const std::vector<CrossSection>& _key_cross_sections) const
{
	// add penalty at a control point when:
	// angle over 180  (concave penalty)
	// angle less than 90 (curve non-increasing penalty)

	double penalty_sum = 0.;
	//int counter = 1;
	bool is_non_increasing_on = false;

	for (auto& iCrossSection : _key_cross_sections)
	{

		//cout << "key" << counter << endl;
		//++counter;

		if (half1_degree_ > 1)
		{
			bool is_convex_half1_point1 = iCrossSection.IsConvex_Half1Point1();
			double angle11 = iCrossSection.AngleAtHalf1Point1();// the vector angle in [0,pi]
			//cout <<"is_convex_half1_point1  "<< is_convex_half1_point1 << ", angle11=" << angle11/M_PI*180 << endl;
			if (!is_convex_half1_point1) // concave penalty
			{
				double extra_angle_over180 = M_PI - angle11;
				penalty_sum += PenaltyOnExtraAnlge(extra_angle_over180);
			}
			else if (is_non_increasing_on)
			{
				if (angle11 < M_PI_2)
				{
					double angle_less90 = M_PI_2 - angle11;
					penalty_sum += PenaltyOnExtraAnlge(angle_less90);
				}
			}

			if (half1_degree_ > 2)
			{
				bool is_convex_half1_point2 = iCrossSection.IsConvex_Half1Point2();
				double angle12 = iCrossSection.AngleAtHalf1Point2();;
				//cout << "is_convex_half1_point2  " << is_convex_half1_point2 << ", angle12=" << angle12 / M_PI * 180 << endl;

				if (!is_convex_half1_point2)
				{
					double extra_angle_over180 = M_PI - angle12;
					penalty_sum += PenaltyOnExtraAnlge(extra_angle_over180);
				}
				// notice this is different from previous
				// for the last point of half1, we add penalty to over 90 angle between (vector01 and vector 23)
				// this is trying to keep the control point in an increasing manner (not turning its head and travel back)
				double angle_between_01and23 = iCrossSection.AngleBetweenHalf1Vector01And23();
				//cout << "angle_between_01and23=  " << angle_between_01and23 / M_PI * 180 << endl;
				if (angle_between_01and23 > M_PI_2)
				{
					double extra_angle_over90 = angle_between_01and23 - M_PI_2;
					penalty_sum += 100.0 * PenaltyOnExtraAnlge(extra_angle_over90);
				}
			}
		}


		if (half2_degree_ > 1)
		{
			bool is_convex_half2_point1 = iCrossSection.IsConvex_Half2Point1();
			double angle21 = iCrossSection.AngleAtHalf2Point1();
			//cout << "is_convex_half2_point1  " << is_convex_half2_point1 << ", angle21=" << angle21 / M_PI * 180 << endl;

			if (!is_convex_half2_point1)
			{
				double extra_angle_over180 = M_PI - angle21;
				penalty_sum += PenaltyOnExtraAnlge(extra_angle_over180);
			}
			else if(is_non_increasing_on)
			{
				if (angle21 < M_PI_2)
				{
					double angle_less90 = M_PI_2 - angle21;
					penalty_sum += PenaltyOnExtraAnlge(angle_less90);
				}
			}

			if (half2_degree_ > 2)
			{
				bool is_convex_half2_point2 = iCrossSection.IsConvex_Half2Point2();
				double angle22 = iCrossSection.AngleAtHalf2Point2();
				//cout << "is_convex_half2_point2  " << is_convex_half2_point2 << ", angle22=" << angle22 / M_PI * 180 << endl;

				if (!is_convex_half2_point2)
				{
					double extra_angle_over180 = M_PI - angle22;
					penalty_sum += PenaltyOnExtraAnlge(extra_angle_over180);
				}
				else
				{
					if (angle22 < M_PI_2)
					{
						double angle_less90 = M_PI_2 - angle22;
						penalty_sum += PenaltyOnExtraAnlge(angle_less90);
					}
				}
			}
		}
	}

	return penalty_sum;
}

double CrossSectionOptimizeInfo::PenaltyOnExtraAnlge(double _extra_angle) const
{
	return _extra_angle * _extra_angle;
}

double CrossSectionOptimizeInfo::Regularity(const std::vector<CrossSection>& _key_cross_sections) const
{
	double penalty_sum = 0.;

	for (auto& iCrossSection : _key_cross_sections)
	{

		if (half1_degree_ > 1)
		{
			std::vector<double> lengths;
			lengths.reserve(2);
			lengths.push_back(iCrossSection.GetHalf1Vector01Length());
			lengths.push_back(iCrossSection.GetHalf1Vector12Length());

			if (half1_degree_ > 2)
			{
				lengths.push_back(iCrossSection.GetHalf1Vector23Length());
			}
			Variance compute_variance;
			penalty_sum += compute_variance.computeVariance(lengths);
		}

		if (half2_degree_ > 1)
		{
			std::vector<double> lengths;
			lengths.reserve(2);
			lengths.push_back(iCrossSection.GetHalf2Vector01Length());
			lengths.push_back(iCrossSection.GetHalf2Vector12Length());

			if (half2_degree_ > 2)
			{
				lengths.push_back(iCrossSection.GetHalf2Vector23Length());
			}
			Variance compute_variance;
			penalty_sum += compute_variance.computeVariance(lengths);
		}
	}

	return penalty_sum;

}

double CrossSectionOptimizeInfo::GlobalLengthVarianceRegularity(const std::vector<CrossSection>& _key_cross_sections) const
{
	std::vector<double> lengths;
	lengths.reserve(_key_cross_sections.size());
	for (auto& iCrossSection : _key_cross_sections)
	{
		lengths.push_back(iCrossSection.GetHalf1EndNorm());
	}
	Variance compute_variance;
	double half1_end_norm_variance = compute_variance.computeVariance(lengths);

	return half1_end_norm_variance;
}