#pragma once
#include "../MaskInterface/CushionSurface.h"
#include "../FEM/TetrahedralCushion.h"
#include "../FEM/FemManager.h"
#include <Discregrid/All>

// TODO provide deep copy function for the private pointer member (FemManager*)
// 
class LinearTetrahedralFemEvaluation
{
public:
	LinearTetrahedralFemEvaluation();
	~LinearTetrahedralFemEvaluation();

	// set cushion and tetrahedralization information
	void SetCushion(const CushionSurface* _cushion);
	void SetFemParameters(double _E, double _nu, const tmd::TriangleMeshDistance* _head_sdf, double _collision_weight,
		int _max_iteration, double _epsilon_gradient, double _epsilon_fucntion, double _epsilon_step, int _sdf_function_index = 0, const Discregrid::CubicLagrangeDiscreteGrid* _real_head_sdf = nullptr, const tmd::TriangleMeshDistance* _simplified_head_sdf = nullptr);

	void SetSigmoidParameters(double _scale, double _epsilon_force);
	void SetObjectiveFunctionWeight(double _weight_average_pressure, double _weight_force_distribution,	double _weight_area_distribution);
	void SetTetrahdedralInversionWeight(double _weight_tertrahedra_inversion);
	void TetrahedralizeCushion(int _num_layer, double _height = 3.0);
	// if cushion surface is changed, just call this again to regenerate tetrahedral cushion
	void ReTetrahedralizeCushion();// use the num_layer and height set previously

	double Evaluate(bool _use_saved_simulation_result = false, std::string _filename = std::string());
	double Evaluate_AccurateForce_WithGradient(Eigen::VectorXd& _gradient);
	double Evaluate_AccurateForce();// copy from Evaluate_AccurateForce_WithGradient..........

	bool IsSimulationGood() const;
	void SaveSimulationResultXml(std::string _filepath, std::string _filename) const;

	// tetrahedral inversion barrier
	// pre-request: call EvaluationPreparation
	// gradient w.r.t. all original nodes
	double TetrahedralInversionBarrierWithGradient(Eigen::VectorXd& _gradient);
	double IpcBarrierFunction_WithDerivative(double _x, double _x0, double& _derivative) const;// -(x-x0)^2*log(x/x0)
	void FindBarrierFunctionX0();
	bool IsInversionFree();

public:
	// for visualization
	void ForceDistributionColor(std::vector<double>& _force_color) const;
	void AreaDistributionColor(std::vector<double>& _area_color) const;
	void SigmoidColor(std::vector<double>& _sigmoid_color) const;
	void GetDeformedCushion(TetrahedralMesh& _tetrahedral_cushion, Mesh& _cushion_surface) const;// input surface mesh alreay have the correct topology
	void GetTetrahedralCushionForces(const TetrahedralMesh& _deformed_mesh, Mesh& _force_mesh, std::vector<double>& _node_force_magnitude, double _scale = 1.0) const;// input surface mesh alreay have the correct topology
	void GetTetrahedralCushionPressure(const std::vector<double>& _node_force_magnitude, std::vector<double>& _node_pressure_magnitude) const;// input all node force, output all node pressure, but only the cushion surface nodes have values
	void GetFreeSurfaceNodesForce(std::vector<double>& _cushion_surface_nodes_force_magnitude) const;
	void GetTetrahedralCushion(TetrahedralMesh* _tetrahedral_cushion) const;// the undeformed cushion

private:
	// mind the correct order to call these functions, some on previous ones  
	void FindCushionSurfaceNodes();
	void ComputeSigmoidForce();
	void ComputeAreasOfAdjacentTriangles();
	void ComputeAreaOfSurfaceNodes(const std::vector<int>& _cushion_surface_nodes, const std::vector<std::vector<int>>* _faces_of_all_nodes);
	double TriangleArea(const Eigen::Vector3d& _node1, const Eigen::Vector3d& _node2, const Eigen::Vector3d& _node3) const;


	// 0.5 + 0.5 * (tanh(scale * (x - epsilon_force)))
	// https://en.wikipedia.org/wiki/Hyperbolic_functions
	double sigmoid(double _x, double _epsilon) const;

	bool EvaluationPreparation(bool _use_saved_simulation_result = false, std::string _filename = std::string(), double _inexact_gradient_epsilon = 0.0);


	//=================================== squared force version (experimented evaluation functions (not used anymore))=======================================
	void ComputeSquaredVector(std::vector<double>& _squared, const std::vector<double>& _original) const;
	void ComputeSigmoidSquaredForce();
	void ComputeMetricAlongU(std::vector<double>& _metric_along_u, const std::vector<double>& _metric_value, const std::vector<double>& _sigmoid) const;

	// adding gradients

	// row i of gradient is w.r.t. node i 
	void ComputeTriangleAreaGradient(double _area, const Eigen::Vector3d& _node0, const Eigen::Vector3d& _node1, const Eigen::Vector3d& _node2, Eigen::Matrix3d& _gradient) const;
	////note compute areas_all_adjacent_triangles_ before calling this
	// and FindCushionSurfaceNodes();
	void ComputeAreasOfAdjacentTrianglesGradient();
	// pre-request: call 	FindCushionSurfaceNodes(); ComputeAreasOfAdjacentTrianglesGradient();
	void ComputeAreaOfSurfaceNodesGradient(const std::vector<std::vector<int>>* _faces_of_all_nodes);
	// pre-request: let FemManager update its force
	void FindGradient_SquaredForce_DeformedAllNodes();
	void ComputeSigmoidSquaredForceGradient();

	// alpha, beta, gamma definition in comments of the class members
	// pre-request: call ComputeSigmoidSquaredForceGradient(), ComputeAreaOfSurfaceNodesGradient, FindGradient_SquaredForce_DeformedAllNodes
	void ComputeGradientBeta_DeformedAllNodes();
	void ComputeGradientBeta_OriginalAllNodes(std::vector<Eigen::SparseMatrix<double>>& _beta_gradient_original_all_nodes);
	// pre-request: call ComputeSigmoidSquaredForceGradient(), FindGradient_SquaredForce_DeformedAllNodes
	void ComputeGradientGamma_DeformedAllNodes();
	// pre-request: call ComputeGradientGamma_DeformedAllNodes()
	// its just negative of gamma_gradient_deformed_all_nodes
	void ComputeGradientGamma_OriginalAllNodes(std::vector<Eigen::SparseMatrix<double>>& _gamma_gradient_original_all_nodes);


	// For the following gradient computation functions, the variable that the gradient is w.r.t. depends on the input gradient (beta and gamma)
	//
	// given gradient of beta and gamma, find the gradient and the gradient of sum beta
	// pre-request: update area_of_surface_nodes and sigmoid_squared_force
	// this is actual the gradient of Average Squared Pressure
	void ComputeAveragePressureGradient(Eigen::VectorXd& _gradient, Eigen::VectorXd& _gradient_sum_beta,
		const std::vector<Eigen::SparseMatrix<double>>& _gamma_gradient,
		const std::vector<Eigen::SparseMatrix<double>>& _beta_gradient,
		double _sum_beta, const std::vector<double>& _alpha, double _sum_alpha) const;
	void ComputeAreaVarianceGradient(Eigen::VectorXd& _gradient, 
		const std::vector<Eigen::SparseMatrix<double>>& _beta_gradient, 
		const Eigen::VectorXd& _gradient_sum_beta, double _sum_beta) const;
	void ComputeForceVarianceGradient(Eigen::VectorXd& _gradient,
		const std::vector<Eigen::SparseMatrix<double>>& _gradient_squared_force, double _sum_force) const;
	// D(|F_i|^2*s_i)
	// note what is the pre-requests
	// this formula is very close to gradient of gamma
	void ComputeGradientSquaredForceMultiplySigmoid(std::vector<Eigen::SparseMatrix<double>>& _gradients,
		const std::vector<Eigen::SparseMatrix<double>>& _sqared_force_gradient) const;

	// utility functions
	// assume FEM is done
	void Evaluate_SquaredForce(double& _sum_beta, double& _sum_force, double& _average_squared_pressure, double& _force_variance,double& _area_variance);

	//=================================== accurate force version=======================================
	void ComputeSigmoidAccurateForce();
	void ComputeSigmoidAccurateForceGradient();
	void ComputeGradientBeta_DeformedAllNodes_type3();
	void ComputeGradientBeta_OriginalAllNodes_type3(std::vector<Eigen::SparseMatrix<double>>& _beta_gradient_original_all_nodes);

	void ComputeAveragePressureGradient_type3(Eigen::VectorXd& _gradient, Eigen::VectorXd& _gradient_sum_beta, Eigen::VectorXd& _gradient_sum_force,
		const std::vector<Eigen::SparseMatrix<double>>& _beta_gradient, const std::vector<Eigen::SparseMatrix<double>>& _squared_force_gradient,
		double _sum_beta, double _average_pressure) const;
	void ComputeForceVarianceGradient_type3(Eigen::VectorXd& _gradient,
		const std::vector<Eigen::SparseMatrix<double>>& _gradient_squared_force, 
		const Eigen::VectorXd& _gradient_sum_force, double _sum_force) const;

private:
	const CushionSurface* cushion_ = nullptr;
	// for FEM simulation
	FemManager* simulation_ = nullptr;
	bool is_first_time_evaluation_ = true;
	double E_ = 0.0;
	double nu_ = 0.0;
	const tmd::TriangleMeshDistance* head_sdf_ = nullptr;
	int sdf_function_index_ = 0;//for CubicLagrangeDiscreteGrid function query
	const Discregrid::CubicLagrangeDiscreteGrid* real_head_sdf_ = nullptr;
	double collision_weight_ = 0.0;
	int max_iteration_ = 0;
	double epsilon_gradient_ = 0.0;
	double epsilon_fucntion_ = 0.0;
	double epsilon_step_ = 0.0;
	int num_layer_ = 0;
	double height_ = 0.0;

	double scale_ = 0.0;// used in sigmoid function
	double epsilon_force_ = 0.0;// used in sigmoid function, force threshold
	double epsilon_force_squared_ = 0.0;// used in sigmoid function, force threshold
	double weight_average_pressure_ = 0.0;
	double weight_force_distribution_ = 0.0;
	double weight_area_distribution_ = 0.0;
	double weight_tertrahedra_inversion_ = 0.0;
	std::vector<double> barrier_function_x0_each_tetrahedra_;// store the parameter for the IPC barrier function

	TetrahedralCushion tetrahedral_cushion_;
	TetrahedralMesh tetrahedral_cushion_mesh_;// read from a saved result, undeformed
	TetrahedralMesh deformed_tetrahedral_cushion_mesh_;

	// nodes are in v major order; id of node (u,v) is u * (num_v - 1) + v (use num_v - 1 since we don't include the connector nodes)
	std::vector<int> cushion_surface_nodes_;// global id of the free nodes; global id of node (u,v) is u * num_v + v; size = num_u * (num_v - 1)
	// following vectors have size same with cushion_surface_nodes_ (if no comment)
	std::vector<double> cushion_surface_nodes_force_magnitude_;
	std::vector<double> cushion_surface_nodes_force_square_magnitude_;
	std::vector<double> sigmoid_force_;
	std::vector<double> sigmoid_squared_force_;
	std::vector<double> sigmoid_accurate_force_;
	std::vector<double> areas_all_adjacent_triangles_;// size = all surface triangles, but only those adjacent to the cushion surface nodes has correct area
	std::vector<double> area_of_surface_nodes_;

	std::vector<double> area_along_u_;// area sum along v; size = sample number of u
	std::vector<double> force_along_u_;
	std::vector<double> squared_force_along_u_;

	// adding gradient
	// 
	// size = all surface triangles; gradient correponding to areas_all_adjacent_triangles_
	// each row in the 3x3 matrix is a gradient to the vertex(local index based on tetrahedral_cushion_)
	std::vector<Eigen::Matrix3d> areas_all_adjacent_triangles_gradient_;

	// for all the gradient members below the stored SparseMatrix object has just one column
	// 
	// gradient w.r.t. deformed coordinates: x (on all nodes, not just surface nodes; but those gradient on the inner nodes should be 0)
	std::vector<Eigen::SparseMatrix<double>> area_of_surface_nodes_gradient_;
	// D|F|^2/Dx
	// the gradient of cushion_surface_nodes_force_square_magnitude_
	// note D|F|^2/Dx0 = -D|F|^2/Dx
	std::vector<Eigen::SparseMatrix<double>> squared_force_gradient_deformed_all_nodes_;
	std::vector<double> sigmoid_squared_force_gradient_;// only the gradient of sigmoid to its variable "x"; the squared force here means we assign squared force value to x
	std::vector<double> sigmoid_accurate_force_gradient_;// only the gradient of sigmoid to its variable "x"; the squared force here means we assign squared force value to x

	// define beta_i = A_i*s_i =  area_of_surface_nodes[i] * sigmoid_squared_force[i]
	// Dbeta/Dx
	std::vector<Eigen::SparseMatrix<double>> beta_gradient_deformed_all_nodes_;

	// define gamma_i = |F_i|^2 * (s_i)^2 = cushion_surface_nodes_force_square_magnitude_[i] * sigmoid_squared_force[i]^2
	// Dgamma_i/Dx
	std::vector<Eigen::SparseMatrix<double>> gamma_gradient_deformed_all_nodes_;

	// define alpha_i = |F_i|^2 * s_i / A_i = gamma_i/beta_i
	// this will be updated in Evaluate_SquaredForce when evaluating average squared pressure
	std::vector<double> alpha_;
	double sum_alpha_ = 0.0;

};