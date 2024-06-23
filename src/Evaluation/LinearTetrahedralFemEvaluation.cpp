#include "LinearTetrahedralFemEvaluation.h"
#include "Variance.h"
#include <windows.h>

#include<Eigen/SparseCholesky>	
#include <finitediff.hpp>

#include <iostream>
using std::cout;
using std::endl;
LinearTetrahedralFemEvaluation::LinearTetrahedralFemEvaluation()
{
}

LinearTetrahedralFemEvaluation::~LinearTetrahedralFemEvaluation()
{
	if (simulation_ != nullptr)
	{
		delete simulation_;
		simulation_ = nullptr;
	}
}

void LinearTetrahedralFemEvaluation::SetCushion(const CushionSurface* _cushion)
{
	cushion_ = _cushion;
}

void LinearTetrahedralFemEvaluation::SetFemParameters(double _E, double _nu, const tmd::TriangleMeshDistance* _head_sdf, double _collision_weight,
	int _max_iteration, double _epsilon_gradient, double _epsilon_fucntion, double _epsilon_step, int _sdf_function_index, const Discregrid::CubicLagrangeDiscreteGrid* _real_head_sdf, const tmd::TriangleMeshDistance* _simplified_head_sdf)
{
	E_ = _E;
	nu_ = _nu;
	head_sdf_ = _head_sdf;
	collision_weight_ = _collision_weight;
	max_iteration_ = _max_iteration;
	epsilon_gradient_ = _epsilon_gradient;
	epsilon_fucntion_ = _epsilon_fucntion;
	epsilon_step_ = _epsilon_step;
	real_head_sdf_ = _real_head_sdf;
	sdf_function_index_ = _sdf_function_index;
}

void LinearTetrahedralFemEvaluation::SetSigmoidParameters(double _scale, double _epsilon_force)
{
	scale_ = _scale;
	epsilon_force_ = _epsilon_force;
	epsilon_force_squared_ = epsilon_force_ * epsilon_force_;
}

void LinearTetrahedralFemEvaluation::SetObjectiveFunctionWeight(double _weight_average_pressure, double _weight_force_distribution, double _weight_area_distribution)
{
	weight_average_pressure_ = _weight_average_pressure;
	weight_force_distribution_ = _weight_force_distribution;
	weight_area_distribution_ = _weight_area_distribution;
}

void LinearTetrahedralFemEvaluation::SetTetrahdedralInversionWeight(double _weight_tertrahedra_inversion)
{
	weight_tertrahedra_inversion_ = _weight_tertrahedra_inversion;
}

double LinearTetrahedralFemEvaluation::Evaluate(bool _use_saved_simulation_result, std::string _filename)
{
	double value = 0.0;
	value = Evaluate_AccurateForce();
	return value;
}


double LinearTetrahedralFemEvaluation::Evaluate_AccurateForce_WithGradient(Eigen::VectorXd& _gradient)
{
	bool is_FEM_good = EvaluationPreparation();
	if (!is_FEM_good)
	{
		return 1e+10;
	}
	LARGE_INTEGER t11, t12, tc1;
	QueryPerformanceFrequency(&tc1);
	QueryPerformanceCounter(&t11);

	// prepare information for evaluation
	simulation_->FindResultNodalForce();
	FindCushionSurfaceNodes();
	simulation_->FindForceManitudeAtGivenNodes(cushion_surface_nodes_, cushion_surface_nodes_force_magnitude_);
	ComputeSquaredVector(cushion_surface_nodes_force_square_magnitude_, cushion_surface_nodes_force_magnitude_);
	ComputeSigmoidAccurateForce();
	ComputeAreasOfAdjacentTriangles();

	// compute area distribution (also get the sum of node contact areas)
	ComputeAreaOfSurfaceNodes(cushion_surface_nodes_, tetrahedral_cushion_.GetFacesOfNodes());
	ComputeMetricAlongU(area_along_u_, area_of_surface_nodes_, sigmoid_accurate_force_);
	Variance compute_variance;
	double sum_beta = 0.0;
	double area_variance = compute_variance.computeVariance(area_along_u_, sum_beta);

	// compute force distribution (also get the sum of forces)
	std::vector<double> ones(sigmoid_accurate_force_.size(), 1.0);
	ComputeMetricAlongU(squared_force_along_u_, cushion_surface_nodes_force_square_magnitude_, ones);// no need add sigmoid
	double sum_force = 0.0;
	double force_variance = compute_variance.computeVariance(squared_force_along_u_, sum_force);

	// compute average pressure
	double average_squared_pressure = sum_force / sum_beta / sum_beta;
	double total_energy = weight_average_pressure_ * average_squared_pressure + weight_force_distribution_ * force_variance + weight_area_distribution_ * area_variance;
	//cout << "average_squared_pressure=" << average_squared_pressure << " ";
	//cout << "force_variance=" << force_variance << " ";
	//cout << "area_variance=" << area_variance << endl;
	//cout << "total_fem_metric =" << total_energy << endl;

	// =================================== gradient ===================================
	// prepare gradient components
	simulation_->UpdateDerivativeInfo();
	ComputeAreasOfAdjacentTrianglesGradient();
	ComputeAreaOfSurfaceNodesGradient(tetrahedral_cushion_.GetFacesOfNodes());
	FindGradient_SquaredForce_DeformedAllNodes();
	ComputeSigmoidAccurateForceGradient();
	ComputeGradientBeta_DeformedAllNodes_type3();


	// find derivative w.r.t. deformed nodes x
	// DE_average_pressure/Dx
	Eigen::VectorXd average_pressure_gradient_deformed_all_nodes;
	Eigen::VectorXd sum_beta_gradient_deformed_all_nodes;
	Eigen::VectorXd sum_force_gradient_deformed_all_nodes;
	ComputeAveragePressureGradient_type3(average_pressure_gradient_deformed_all_nodes, 
		sum_beta_gradient_deformed_all_nodes, sum_force_gradient_deformed_all_nodes,
		beta_gradient_deformed_all_nodes_, squared_force_gradient_deformed_all_nodes_, sum_beta, average_squared_pressure);

	// DE_area_variance/Dx
	Eigen::VectorXd area_variance_gradient_deformed_all_nodes;

	ComputeAreaVarianceGradient(area_variance_gradient_deformed_all_nodes, beta_gradient_deformed_all_nodes_, sum_beta_gradient_deformed_all_nodes, sum_beta);

	// DE_force_variance/Dx
	Eigen::VectorXd force_variance_gradient_deformed_all_nodes;
	ComputeForceVarianceGradient_type3(force_variance_gradient_deformed_all_nodes, squared_force_gradient_deformed_all_nodes_, sum_force_gradient_deformed_all_nodes, sum_force);

	// DE_all/Dx
	Eigen::VectorXd weighted_average_pressure_gradient_deformed_all_nodes = weight_average_pressure_ * average_pressure_gradient_deformed_all_nodes;
	Eigen::VectorXd weighted_area_variance_gradient_deformed_all_nodes = weight_area_distribution_ * area_variance_gradient_deformed_all_nodes;
	Eigen::VectorXd weighted_force_variance_gradient_deformed_all_nodes = weight_force_distribution_ * force_variance_gradient_deformed_all_nodes;
	Eigen::VectorXd total_enegy_gradient_deformed_all_nodes =
		weighted_average_pressure_gradient_deformed_all_nodes +
		weighted_area_variance_gradient_deformed_all_nodes +
		weighted_force_variance_gradient_deformed_all_nodes;

	// find derivative w.r.t. original nodes x0
	std::vector<Eigen::SparseMatrix<double>> beta_gradient_original_all_nodes;
	ComputeGradientBeta_OriginalAllNodes_type3(beta_gradient_original_all_nodes);
	// DE_average_pressure/Dx0
	std::vector<Eigen::SparseMatrix<double>> squared_force_gradient_original_all_nodes = squared_force_gradient_deformed_all_nodes_;
	for (auto& i : squared_force_gradient_original_all_nodes)
	{
		i = -i;
	}
	Eigen::VectorXd average_pressure_gradient_original_all_nodes;
	Eigen::VectorXd sum_beta_gradient_original_all_nodes;
	Eigen::VectorXd sum_force_gradient_original_all_nodes;

	ComputeAveragePressureGradient_type3(average_pressure_gradient_original_all_nodes,
		sum_beta_gradient_original_all_nodes, sum_force_gradient_original_all_nodes,
		beta_gradient_original_all_nodes, squared_force_gradient_original_all_nodes, sum_beta, average_squared_pressure);
	// DE_area_variance/Dx0
	Eigen::VectorXd area_variance_gradient_original_all_nodes;

	ComputeAreaVarianceGradient(area_variance_gradient_original_all_nodes, beta_gradient_original_all_nodes, sum_beta_gradient_original_all_nodes, sum_beta);
	// DE_force_variance/Dx0 = - DE_force_variance/Dx

	Eigen::VectorXd force_variance_gradient_original_all_nodes = -force_variance_gradient_deformed_all_nodes;

	// DE_all/Dx0
	Eigen::VectorXd weighted_average_pressure_gradient_original_all_nodes = weight_average_pressure_ * average_pressure_gradient_original_all_nodes;
	Eigen::VectorXd weighted_area_variance_gradient_original_all_nodes = weight_area_distribution_ * area_variance_gradient_original_all_nodes;
	Eigen::VectorXd weighted_force_variance_gradient_original_all_nodes = weight_force_distribution_ * force_variance_gradient_original_all_nodes;
	Eigen::VectorXd total_enegy_gradient_original_all_nodes =
		weighted_average_pressure_gradient_original_all_nodes +
		weighted_area_variance_gradient_original_all_nodes +
		weighted_force_variance_gradient_original_all_nodes;


	// denote w = (DE_all/Dx * Dx/Dx0)^T (DE_all/Dx is regarded as a row vector, w is a column vector)
	Eigen::VectorXd w;

	Eigen::VectorXd w1;//for debug
	Eigen::VectorXd w2;
	Eigen::VectorXd w3;
	{
		// compute w this way:
		// P * Dx/Dx0 = M1^(-1) * M2 (P is a row permutation matrix)
		// 
		// M1^T * z = P * (DE_all/Dx)^T (note M1 is positive-definite)
		// w = M2^T * z
		Eigen::SparseMatrix<double> M1;
		Eigen::SparseMatrix<double> M2;
		simulation_->Derivative_Deformed_Original_Term1(M1);
		simulation_->Derivative_Deformed_Original_Term2(M2);
		const std::vector<int>& P = simulation_->GetDerivativeVariableRowOrderingMap();

		// construct P * (DE_all/Dx)^T
		Eigen::VectorXd permuted_total_enegy_gradient_deformed_all_nodes;
		permuted_total_enegy_gradient_deformed_all_nodes.resize(total_enegy_gradient_deformed_all_nodes.rows());
		for (int i = 0; i < total_enegy_gradient_deformed_all_nodes.rows(); ++i)
		{
			permuted_total_enegy_gradient_deformed_all_nodes(P[i]) = total_enegy_gradient_deformed_all_nodes(i);
		}

		//for debug
		Eigen::VectorXd permuted_average_pressure_gradient_deformed_all_nodes;
		permuted_average_pressure_gradient_deformed_all_nodes.resize(weighted_average_pressure_gradient_deformed_all_nodes.rows());
		for (int i = 0; i < weighted_average_pressure_gradient_deformed_all_nodes.rows(); ++i)
		{
			permuted_average_pressure_gradient_deformed_all_nodes(P[i]) = weighted_average_pressure_gradient_deformed_all_nodes(i);
		}
		Eigen::VectorXd permuted_area_variance_gradient_deformed_all_nodes;
		permuted_area_variance_gradient_deformed_all_nodes.resize(weighted_area_variance_gradient_deformed_all_nodes.rows());
		for (int i = 0; i < weighted_area_variance_gradient_deformed_all_nodes.rows(); ++i)
		{
			permuted_area_variance_gradient_deformed_all_nodes(P[i]) = weighted_area_variance_gradient_deformed_all_nodes(i);
		}
		Eigen::VectorXd permuted_force_variance_gradient_deformed_all_nodes;
		permuted_force_variance_gradient_deformed_all_nodes.resize(weighted_force_variance_gradient_deformed_all_nodes.rows());
		for (int i = 0; i < weighted_force_variance_gradient_deformed_all_nodes.rows(); ++i)
		{
			permuted_force_variance_gradient_deformed_all_nodes(P[i]) = weighted_force_variance_gradient_deformed_all_nodes(i);
		}


		// solve M1 * z = P * (DE_all/Dx)^T
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver_M1;
		solver_M1.compute(M1);
		if (solver_M1.info() != Eigen::Success) {
			// decomposition failed
			cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ERROR LinearTetrahedralFemEvaluation::Evaluate_AccurateForce_WithGradient M1 decomposition failed" << endl;
			/*return 1e+10;*/
		}
		if (solver_M1.info() == Eigen::NumericalIssue)
		{
			cout << "WARNING!!!!!!!! LinearTetrahedralFemEvaluation::Evaluate_AccurateForce_WithGradient M1 may not positive-definite" << endl;
		}
		Eigen::VectorXd z = solver_M1.solve(permuted_total_enegy_gradient_deformed_all_nodes);
		if (solver_M1.info() != Eigen::Success) {
			// solving failed
			cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ERROR LinearTetrahedralFemEvaluation::Evaluate_AccurateForce_WithGradient M1 solving failed" << endl;
			return 1e+10;
		}

		////for debug
		//Eigen::VectorXd z1 = solver_M1.solve(permuted_average_pressure_gradient_deformed_all_nodes);
		//Eigen::VectorXd z2 = solver_M1.solve(permuted_area_variance_gradient_deformed_all_nodes);
		//Eigen::VectorXd z3 = solver_M1.solve(permuted_force_variance_gradient_deformed_all_nodes);
		//w1 = M2.transpose() * z1;
		//w2 = M2.transpose() * z2;
		//w3 = M2.transpose() * z3;


		// w = M2^T * z
		w = M2.transpose() * z;
	}

	// Dx0/Dp
	Eigen::MatrixXd Dx0_Dp;
	tetrahedral_cushion_.GetDerivative_Nodes_CrossSectionKeyParameters(Dx0_Dp);

	// note w and DE_all/Dx0 is stored as column Eigen::Vector
	// DE_all/Dp = ((w + DE_all/Dx0).transpose * Dx0/Dp).transpose
	//before that, we add gradient of tetrahedral inversion term energy
	Eigen::VectorXd gradient_tetradedral_intersion;
	double tetradedral_intersion = TetrahedralInversionBarrierWithGradient(gradient_tetradedral_intersion);
	_gradient = ((w + total_enegy_gradient_original_all_nodes + weight_tertrahedra_inversion_ * gradient_tetradedral_intersion).transpose() * Dx0_Dp).transpose();


	//// debug
	//Eigen::VectorXd gradient_tertrahedra_inversion = (weight_tertrahedra_inversion_ * gradient_tetradedral_intersion).transpose() * Dx0_Dp;
	//Eigen::VectorXd gradient_fem_metric = (w + total_enegy_gradient_original_all_nodes).transpose() * Dx0_Dp;
	//Eigen::VectorXd gradient_fem_metric1 = (w1 + weighted_average_pressure_gradient_original_all_nodes).transpose() * Dx0_Dp;
	//Eigen::VectorXd gradient_fem_metric2 = (w2 + weighted_area_variance_gradient_original_all_nodes).transpose() * Dx0_Dp;
	//Eigen::VectorXd gradient_fem_metric3 = (w3 + weighted_force_variance_gradient_original_all_nodes).transpose() * Dx0_Dp;
	////cout << "_gradient" << _gradient.transpose() << endl;
	////cout << "gradient_fem_metric" << gradient_fem_metric.transpose() << endl;
	////cout << "gradient_tertrahedra_inversion" << gradient_tertrahedra_inversion.transpose() << endl;

	////cout << "gradient_fem_metric norm=" << gradient_fem_metric.norm()
	////	<< "\n gradient_average_pressure norm=" << gradient_fem_metric1.norm() /*<< "\ngradient_w1 norm=" << w1.norm() << "\nweighted_average_pressure_gradient_original_all_nodes norm="<< weighted_average_pressure_gradient_original_all_nodes.norm()*/
	////	<< "\n gradient_area_variance norm=" << gradient_fem_metric2.norm() /*<< "\ngradient_w2 norm=" << w2.norm() << "\nweighted_area_variance_gradient_original_all_nodes norm=" << weighted_area_variance_gradient_original_all_nodes.norm()*/
	////	<< "\n gradient_force_variance norm=" << gradient_fem_metric3.norm() /*<< "\ngradient_w3 norm=" << w3.norm() << "\nweighted_force_variance_gradient_original_all_nodes norm=" << weighted_force_variance_gradient_original_all_nodes.norm()*/
	////	<< "\ngradient_tertrahedra_inversion norm=" << gradient_tertrahedra_inversion.norm() << endl;
	////cout << "weighted tetradedral_intersion=" << weight_tertrahedra_inversion_ * tetradedral_intersion << endl;

	//QueryPerformanceCounter(&t12);
	////cout << "gradient computation time = " << (double)(t12.QuadPart - t11.QuadPart) / (double)tc1.QuadPart << "s" << endl;

	return total_energy + weight_tertrahedra_inversion_ * tetradedral_intersion;
}

double LinearTetrahedralFemEvaluation::Evaluate_AccurateForce()
{
	bool is_FEM_good = EvaluationPreparation();
	if (!is_FEM_good)
	{
		return 1e+10;
	}
	LARGE_INTEGER t11, t12, tc1;
	QueryPerformanceFrequency(&tc1);
	QueryPerformanceCounter(&t11);

	// prepare information for evaluation
	simulation_->FindResultNodalForce();
	FindCushionSurfaceNodes();
	simulation_->FindForceManitudeAtGivenNodes(cushion_surface_nodes_, cushion_surface_nodes_force_magnitude_);
	ComputeSquaredVector(cushion_surface_nodes_force_square_magnitude_, cushion_surface_nodes_force_magnitude_);
	ComputeSigmoidAccurateForce();
	ComputeAreasOfAdjacentTriangles();

	// compute area distribution (also get the sum of node contact areas)
	ComputeAreaOfSurfaceNodes(cushion_surface_nodes_, tetrahedral_cushion_.GetFacesOfNodes());
	ComputeMetricAlongU(area_along_u_, area_of_surface_nodes_, sigmoid_accurate_force_);
	Variance compute_variance;
	double sum_beta = 0.0;
	double area_variance = compute_variance.computeVariance(area_along_u_, sum_beta);

	// compute force distribution (also get the sum of forces)
	std::vector<double> ones(sigmoid_accurate_force_.size(), 1.0);
	ComputeMetricAlongU(squared_force_along_u_, cushion_surface_nodes_force_square_magnitude_, ones);// no need add sigmoid
	double sum_force = 0.0;
	double force_variance = compute_variance.computeVariance(squared_force_along_u_, sum_force);

	// compute average pressure
	double average_squared_pressure = sum_force / sum_beta / sum_beta;
	double total_energy = weight_average_pressure_ * average_squared_pressure + weight_force_distribution_ * force_variance + weight_area_distribution_ * area_variance;
	cout << "average_squared_pressure=" << average_squared_pressure << " ";
	cout << "force_variance=" << force_variance << " ";
	cout << "area_variance=" << area_variance << endl;
	cout << "total_fem_metric =" << total_energy << endl;
	return total_energy;
}

bool LinearTetrahedralFemEvaluation::IsSimulationGood() const
{
	bool is_it_good = false;
	if (simulation_ != nullptr)
	{
		is_it_good = simulation_->IsElementPositive();
	}
	return is_it_good;
}

void LinearTetrahedralFemEvaluation::SaveSimulationResultXml(std::string _filepath, std::string _filename) const
{
	if (simulation_ != nullptr)
	{
		simulation_->SaveResultXml(_filepath, _filename);
	}
}

double LinearTetrahedralFemEvaluation::TetrahedralInversionBarrierWithGradient(Eigen::VectorXd& _gradient)
{
	// use 6*volume determine the inversion

	// get element-wise values and gradient
	const TetrahedralMesh* tet_cushion = tetrahedral_cushion_.GetTetrahedralCushion();
	int num_element = tet_cushion->GetTopologyTetrahedra()->rows();
	std::vector<double> volumes6(num_element);//6 times of volumes
	std::vector<std::array<Eigen::Vector3d, 4>> gradient_elements(num_element);
	for (int iTetrahedra = 0; iTetrahedra < num_element; ++iTetrahedra)
	{
		// get vertex
		std::array<Eigen::Vector3d, 4> vertex;
		tet_cushion->GetTetrahedron(vertex, iTetrahedra);

		// vectors
		Eigen::Vector3d vector01 = vertex[1] - vertex[0];
		Eigen::Vector3d vector02 = vertex[2] - vertex[0];
		Eigen::Vector3d vector03 = vertex[3] - vertex[0];

		//6*vol
		volumes6[iTetrahedra] = (vector01.cross(vector02)).dot(vector03);

		Eigen::Vector3d vector13 = vertex[3] - vertex[1];
		Eigen::Vector3d vector12 = vertex[2] - vertex[1];
		//w.r.t. vertex0
		gradient_elements[iTetrahedra][0] = vector13.cross(vector12);
		//w.r.t. vertex1
		gradient_elements[iTetrahedra][1] = vector02.cross(vector03);
		//w.r.t. vertex2
		gradient_elements[iTetrahedra][2] = vector03.cross(vector01);
		//w.r.t. vertex3
		gradient_elements[iTetrahedra][3] = vector01.cross(vector02);
	}

	// values
	// sum of -log(6*vol)
	double sum = 0.0;
	std::vector<double> barrier_graident(num_element);
	for (int iTetrahedra = 0; iTetrahedra < num_element; ++iTetrahedra)
	{
		sum += IpcBarrierFunction_WithDerivative(volumes6[iTetrahedra], barrier_function_x0_each_tetrahedra_[iTetrahedra], barrier_graident[iTetrahedra]);
	}

	// arrange gradient of each tetrahdedron to gradient of the sum w.r.t. nodes
	_gradient.setZero(3 * tetrahedral_cushion_.GetNodeNumber());
	for (int iTetrahedra = 0; iTetrahedra < num_element; ++iTetrahedra)
	{
		for (int iVertex = 0; iVertex < 4; ++iVertex)
		{
			int node_id = (*tet_cushion->GetTopologyTetrahedra())(iTetrahedra, iVertex);
			//_gradient.block(3 * node_id, 0, 3, 1) += gradient_elements[iTetrahedra][iVertex] / (-volumes6[iTetrahedra]);
			_gradient.block(3 * node_id, 0, 3, 1) += gradient_elements[iTetrahedra][iVertex] * barrier_graident[iTetrahedra];
		}
	}

	// note, now gradient w.r.t. original nodes x0, will need to multiply it with Dx0/Dp to get the final form

	return sum;
}

double LinearTetrahedralFemEvaluation::IpcBarrierFunction_WithDerivative(double _x, double _x0, double& _derivative) const
{
	double value = 1e10;
	if (_x > 0.0)
	{
		if (_x >= _x0)
		{
			value = 0.0;
			_derivative = 0.0;
		}
		else
		{
			value = -(_x - _x0) * (_x - _x0) * log(_x / _x0);
			// (x0 - x)(2log(x/x0)-x0/x+1)
			_derivative = (_x0 - _x) * (2.0 * log(_x / _x0) - _x0 / _x + 1.0);
		}
	}
	else
	{
		cout << "WARNING!!!!!!!!!!! CrossSectionOptimizeInfo::BarrierFunction_WithDerivative x<=0, exceed barrier" << endl;
	}
	return value;
}

void LinearTetrahedralFemEvaluation::FindBarrierFunctionX0()
{
	// get volumes (copied from TetrahedralInversionBarrierWithGradient)
	const TetrahedralMesh* tet_cushion = tetrahedral_cushion_.GetTetrahedralCushion();
	int num_element = tet_cushion->GetTopologyTetrahedra()->rows();
	barrier_function_x0_each_tetrahedra_.resize(num_element);
	for (int iTetrahedra = 0; iTetrahedra < num_element; ++iTetrahedra)
	{
		// get vertex
		std::array<Eigen::Vector3d, 4> vertex;
		tet_cushion->GetTetrahedron(vertex, iTetrahedra);

		// vectors
		Eigen::Vector3d vector01 = vertex[1] - vertex[0];
		Eigen::Vector3d vector02 = vertex[2] - vertex[0];
		Eigen::Vector3d vector03 = vertex[3] - vertex[0];

		barrier_function_x0_each_tetrahedra_[iTetrahedra] = (vector01.cross(vector02)).dot(vector03)/1.0;
		//cout << "(vector01.cross(vector02)).dot(vector03)=" << (vector01.cross(vector02)).dot(vector03) << endl;
	}
}

bool LinearTetrahedralFemEvaluation::IsInversionFree()
{
	// use 6*volume determine the inversion

	const TetrahedralMesh* tet_cushion = tetrahedral_cushion_.GetTetrahedralCushion();
	int num_element = tet_cushion->GetTopologyTetrahedra()->rows();
	for (int iTetrahedra = 0; iTetrahedra < num_element; ++iTetrahedra)
	{
		// get vertex
		std::array<Eigen::Vector3d, 4> vertex;
		tet_cushion->GetTetrahedron(vertex, iTetrahedra);

		// vectors
		Eigen::Vector3d vector01 = vertex[1] - vertex[0];
		Eigen::Vector3d vector02 = vertex[2] - vertex[0];
		Eigen::Vector3d vector03 = vertex[3] - vertex[0];

		//6*vol
		double volume6 = (vector01.cross(vector02)).dot(vector03);

		if (volume6 < 1e-10)
		{
			return false;
		}
	}
	return true;
}

void LinearTetrahedralFemEvaluation::ForceDistributionColor(std::vector<double>& _force_color) const
{
	std::vector<double> force_along_u;

	force_along_u = squared_force_along_u_;
	
	_force_color.resize(tetrahedral_cushion_.GetUNumber() * tetrahedral_cushion_.GetVNumber());
	for (int iU = 0; iU < tetrahedral_cushion_.GetUNumber(); ++iU)
	{
		for (int iV = 0; iV < tetrahedral_cushion_.GetVNumber(); ++iV)
		{
			_force_color[iU * tetrahedral_cushion_.GetVNumber() + iV] = force_along_u[iU];
		}
	}
}

void LinearTetrahedralFemEvaluation::AreaDistributionColor(std::vector<double>& _area_color) const
{

	_area_color.resize(tetrahedral_cushion_.GetUNumber() * tetrahedral_cushion_.GetVNumber());
	for (int iU = 0; iU < tetrahedral_cushion_.GetUNumber(); ++iU)
	{
		for (int iV = 0; iV < tetrahedral_cushion_.GetVNumber(); ++iV)
		{
			_area_color[iU * tetrahedral_cushion_.GetVNumber() + iV] = area_along_u_[iU];
		}
	}
}

void LinearTetrahedralFemEvaluation::SigmoidColor(std::vector<double>& _sigmoid_color) const
{
	std::vector<double> sigmoid_force;
	sigmoid_force = sigmoid_accurate_force_;

	_sigmoid_color.resize(tetrahedral_cushion_.GetUNumber() * tetrahedral_cushion_.GetVNumber());
	for (int iU = 0; iU < tetrahedral_cushion_.GetUNumber(); ++iU)
	{
		for (int iV = 0; iV < tetrahedral_cushion_.GetVNumber() - 1; ++iV)
		{
			_sigmoid_color[iU * tetrahedral_cushion_.GetVNumber() + iV] = sigmoid_force[iU * (tetrahedral_cushion_.GetVNumber() - 1) + iV];
		}
		_sigmoid_color[iU * tetrahedral_cushion_.GetVNumber() + tetrahedral_cushion_.GetVNumber() - 1] = 0;
	}
}

void LinearTetrahedralFemEvaluation::GetDeformedCushion(TetrahedralMesh& _tetrahedral_cushion, Mesh& _cushion_surface) const
{
	if (simulation_ != nullptr)
	{

		simulation_->GetDeformedMesh(_tetrahedral_cushion);

		for (int iU = 0; iU < tetrahedral_cushion_.GetUNumber(); ++iU)
		{
			for (int iV = 0; iV < tetrahedral_cushion_.GetVNumber(); ++iV)
			{
				int vertex_id = iU * tetrahedral_cushion_.GetVNumber() + iV;
				Eigen::Vector3d vertex;
				_tetrahedral_cushion.GetVertex(vertex_id, vertex);
				_cushion_surface.set_point(_cushion_surface.vertex_handle(vertex_id), Mesh::Point(vertex(0), vertex(1), vertex(2)));
			}
		}
	}
}

void LinearTetrahedralFemEvaluation::GetTetrahedralCushionForces(const TetrahedralMesh& _deformed_mesh, Mesh& _force_mesh, std::vector<double>& _node_force_magnitude, double _scale) const
{
	if (simulation_ != nullptr)
	{
		simulation_->FindResultNodalForce();
		simulation_->FindResultSurfaceForceLineMeshAndAllForceMagnitude(_deformed_mesh, _force_mesh, _node_force_magnitude, _scale);
		// set the force of non-free nodes to zero
		std::vector<bool> is_vertex_non_free(_force_mesh.n_faces(), true);// notice face of this mesh correspond to the tetrahedral boundary nodes
		for (int iFree = 0; iFree < cushion_surface_nodes_.size(); ++iFree)
		{
			is_vertex_non_free[cushion_surface_nodes_[iFree]] = false;
		}
		for (int iVertex = 0; iVertex < _force_mesh.n_faces(); ++iVertex)
		{
			if (is_vertex_non_free[iVertex] == true)
			{
				_force_mesh.delete_face(_force_mesh.face_handle(iVertex), true);
			}
		}
		_force_mesh.garbage_collection();
		// set the force magnitude of non-free nodes to zero
		std::vector<int> node_fixed;
		simulation_->GetBoundaryConditionNodes(node_fixed);
		//cout << "node_fixed" << node_fixed.size() << endl;
		for (int iFix=0; iFix < node_fixed.size(); ++iFix)
		{
			_node_force_magnitude[node_fixed[iFix]] = 0.0;
		}
	}
}

void LinearTetrahedralFemEvaluation::GetTetrahedralCushionPressure(const std::vector<double>& _node_force_magnitude, std::vector<double>& _node_pressure_magnitude) const
{
	_node_pressure_magnitude.resize(_node_force_magnitude.size(), 0.0);
	for (int iNode = 0; iNode < cushion_surface_nodes_.size(); ++iNode)
	{
		_node_pressure_magnitude[cushion_surface_nodes_[iNode]] = 3.0 * _node_force_magnitude[cushion_surface_nodes_[iNode]] / area_of_surface_nodes_[iNode];
		//cout << "area_of_surface_nodes_[iNode]=" << area_of_surface_nodes_[iNode] << endl;
	}
}

void LinearTetrahedralFemEvaluation::GetFreeSurfaceNodesForce(std::vector<double>& _cushion_surface_nodes_force_magnitude) const
{
	_cushion_surface_nodes_force_magnitude = cushion_surface_nodes_force_magnitude_;
}

void LinearTetrahedralFemEvaluation::GetTetrahedralCushion(TetrahedralMesh* _tetrahedral_cushion) const
{
	if (tetrahedral_cushion_mesh_.GetVertex()->rows() != 0)
	{
		*_tetrahedral_cushion = tetrahedral_cushion_mesh_;
	}
	else
	{
		*_tetrahedral_cushion = *tetrahedral_cushion_.GetTetrahedralCushion();
	}
}



void LinearTetrahedralFemEvaluation::FindCushionSurfaceNodes()
{
	int num_v_free = (tetrahedral_cushion_.GetVNumber() - 1);
	int surface_node_num = tetrahedral_cushion_.GetUNumber() * num_v_free;
	if (cushion_surface_nodes_.size() != surface_node_num)
	{
		cushion_surface_nodes_.resize(surface_node_num);
		int count_node = 0;

		for (int iU = 0; iU < tetrahedral_cushion_.GetUNumber(); ++iU)
		{
			for (int iV = 0; iV < num_v_free; ++iV)
			{
				cushion_surface_nodes_[count_node] = iU * tetrahedral_cushion_.GetVNumber() + iV;
				++count_node;
			}
		}
	}
	
}

void LinearTetrahedralFemEvaluation::ComputeSigmoidForce()
{
	sigmoid_force_.resize(cushion_surface_nodes_force_magnitude_.size());
	for (int iNode = 0; iNode < cushion_surface_nodes_force_magnitude_.size(); ++iNode)
	{
		sigmoid_force_[iNode] = sigmoid(cushion_surface_nodes_force_magnitude_[iNode], epsilon_force_);
	}
}

void LinearTetrahedralFemEvaluation::ComputeAreasOfAdjacentTriangles()
{
	areas_all_adjacent_triangles_.clear();
	areas_all_adjacent_triangles_.resize(tetrahedral_cushion_.GetTriangleNumber(), -1);
	const std::vector<std::vector<int>>* faces_of_all_nodes = tetrahedral_cushion_.GetFacesOfNodes();
	for (int iNode = 0; iNode < cushion_surface_nodes_.size(); ++iNode)
	{
		const std::vector<int>& adjacent_faces = faces_of_all_nodes->at(cushion_surface_nodes_[iNode]);
		for (int iFace = 0; iFace < adjacent_faces.size(); ++iFace)
		{
			if (areas_all_adjacent_triangles_[adjacent_faces[iFace]] < 0)
			{
				int node1_id;
				int node2_id;
				int node3_id;

				tetrahedral_cushion_.GetNodesIdOfSurfaceFace(node1_id, node2_id, node3_id, adjacent_faces[iFace]);

				Eigen::Vector3d node1;
				Eigen::Vector3d node2;
				Eigen::Vector3d node3;
				deformed_tetrahedral_cushion_mesh_.GetVertex(node1_id, node1);
				deformed_tetrahedral_cushion_mesh_.GetVertex(node2_id, node2);
				deformed_tetrahedral_cushion_mesh_.GetVertex(node3_id, node3);
				areas_all_adjacent_triangles_[adjacent_faces[iFace]] = TriangleArea(node1, node2, node3);

			}
		}
	}
}

void LinearTetrahedralFemEvaluation::ComputeAreaOfSurfaceNodes(const std::vector<int>& _cushion_surface_nodes, const std::vector<std::vector<int>>* _faces_of_all_nodes)
{
	area_of_surface_nodes_.resize(_cushion_surface_nodes.size());
	for (int iNode = 0; iNode < _cushion_surface_nodes.size(); ++iNode)
	{
		double area = 0.0;
		const std::vector<int>& adjacent_faces = _faces_of_all_nodes->at(_cushion_surface_nodes[iNode]);
		for (int iFace = 0; iFace < adjacent_faces.size(); ++iFace)
		{
			area += areas_all_adjacent_triangles_[adjacent_faces[iFace]];
			assert(areas_all_adjacent_triangles_[adjacent_faces[iFace]] > 0);
		}
		area_of_surface_nodes_[iNode] = area;
	}
}

double LinearTetrahedralFemEvaluation::TriangleArea(const Eigen::Vector3d& _node1, const Eigen::Vector3d& _node2, const Eigen::Vector3d& _node3) const
{
	return 0.5 * ((_node2 - _node1).cross(_node3 - _node1)).norm();
}


double LinearTetrahedralFemEvaluation::sigmoid(double _x, double _epsilon) const
{
	return 0.5 + 0.5 * std::tanh(scale_ * (_x - _epsilon));
}

bool LinearTetrahedralFemEvaluation::EvaluationPreparation(bool _use_saved_simulation_result, std::string _filename, double _inexact_gradient_epsilon)
{
	bool is_FEM_good = true;// for indicating if last tetrahedral has no negative volume
	double gradient_epsilon = epsilon_gradient_;
	if (_inexact_gradient_epsilon > 0.0)
	{
		gradient_epsilon = _inexact_gradient_epsilon;
	}

	// simulation
	if (simulation_ != nullptr)
	{

		delete simulation_;
		simulation_ = nullptr;
	}
	simulation_ = new FemManager;
	if (!_use_saved_simulation_result)
	{
		simulation_->SetTetrahedralMesh(tetrahedral_cushion_.GetTetrahedralCushion());
		simulation_->SetMaterial(E_, nu_);
		std::vector<std::pair<int, Eigen::Vector3d>> displacement_boundary;
		tetrahedral_cushion_.GetFixConnectorBoundaryCondition(displacement_boundary);
		
		simulation_->SetBoundaryCondition(displacement_boundary);
		simulation_->SetCollisionInfo(head_sdf_, collision_weight_, sdf_function_index_, real_head_sdf_);
		if (deformed_tetrahedral_cushion_mesh_.GetVertex()->rows() == 0)
		{
			// initial displacement is set as zero
			is_FEM_good = simulation_->Compute_LBFGS_PenaltyCollision(Eigen::Vector3d::Zero(), max_iteration_, gradient_epsilon, epsilon_fucntion_, epsilon_step_);

			if (is_first_time_evaluation_)
			{
				is_first_time_evaluation_ = false;
			}
		}
		else
		{
			// initial displacement set as the last simulation result configuration
			is_FEM_good = simulation_->Compute_LBFGS_PenaltyCollision(&deformed_tetrahedral_cushion_mesh_, max_iteration_, gradient_epsilon, epsilon_fucntion_, epsilon_step_);
		}
	}
	else
	{
		tetrahedral_cushion_mesh_ = *tetrahedral_cushion_.GetTetrahedralCushion();
		simulation_->ReadResultXml(_filename, tetrahedral_cushion_mesh_);
	}

	if (is_FEM_good == false)
	{
		cout << "NOTICE from LinearTetrahedralFemEvaluation::EvaluationPreparation, exists an element Volume not Positive" << endl;
	}
	else 
	{
		// get the deformed cushion mesh of this simulation
		simulation_->GetDeformedMesh(deformed_tetrahedral_cushion_mesh_);
	}
	

	return is_FEM_good;
}

void LinearTetrahedralFemEvaluation::ComputeSquaredVector(std::vector<double>& _squared, const std::vector<double>& _original) const
{
	_squared.resize(_original.size());
	for (int iNode = 0; iNode < _original.size(); ++iNode)
	{
		_squared[iNode] = _original[iNode] * _original[iNode];
	}
}

void LinearTetrahedralFemEvaluation::ComputeSigmoidSquaredForce()
{
	sigmoid_squared_force_.resize(cushion_surface_nodes_force_square_magnitude_.size());
	for (int iNode = 0; iNode < cushion_surface_nodes_force_square_magnitude_.size(); ++iNode)
	{
		sigmoid_squared_force_[iNode] = sigmoid(cushion_surface_nodes_force_square_magnitude_[iNode], epsilon_force_squared_);
	}
}

void LinearTetrahedralFemEvaluation::ComputeMetricAlongU(std::vector<double>& _metric_along_u, const std::vector<double>& _metric_value, const std::vector<double>& _sigmoid) const
{	
	_metric_along_u.resize(tetrahedral_cushion_.GetUNumber());

	int num_v_free = tetrahedral_cushion_.GetVNumber() - 1;
	for (int iU = 0; iU < tetrahedral_cushion_.GetUNumber(); ++iU)
	{
		_metric_along_u[iU] = 0.0;
		for (int iV = 0; iV < num_v_free; ++iV)
		{
			int free_node_id = iU * num_v_free + iV;
			_metric_along_u[iU] += _metric_value[free_node_id] * _sigmoid[free_node_id];
		}
	}
	
}

void LinearTetrahedralFemEvaluation::ComputeTriangleAreaGradient(double _area, const Eigen::Vector3d& _node0, const Eigen::Vector3d& _node1, const Eigen::Vector3d& _node2, Eigen::Matrix3d& _gradient) const
{
	Eigen::Vector3d edge01 = _node1 - _node0;
	Eigen::Vector3d edge12 = _node2 - _node1;
	Eigen::Vector3d edge20 = _node0 - _node2;
	Eigen::Vector3d edge10 = -edge01;
	Eigen::Vector3d edge21 = -edge12;
	Eigen::Vector3d edge02 = -edge20;

	double temp = 4.0 * _area;
	_gradient.block(0, 0, 1, 3) = ((edge01.cross(edge02)).cross(edge12) / temp).transpose();
	_gradient.block(1, 0, 1, 3) = ((edge10.cross(edge12)).cross(edge02) / temp).transpose();
	_gradient.block(2, 0, 1, 3) = ((edge20.cross(edge21)).cross(edge01) / temp).transpose();
}

void LinearTetrahedralFemEvaluation::ComputeAreasOfAdjacentTrianglesGradient()
{
	// all area intialized to -1; adjacent ones will be computed, others remains -1
	areas_all_adjacent_triangles_gradient_.clear();
	areas_all_adjacent_triangles_gradient_.resize(tetrahedral_cushion_.GetTriangleNumber());
	const std::vector<std::vector<int>>* faces_of_all_nodes = tetrahedral_cushion_.GetFacesOfNodes();
	std::vector<bool> is_face_visited(tetrahedral_cushion_.GetTriangleNumber(), false);// indicate whether the gradient is computed for a triangle
	for (int iNode = 0; iNode < cushion_surface_nodes_.size(); ++iNode)
	{
		const std::vector<int>& adjacent_faces = faces_of_all_nodes->at(cushion_surface_nodes_[iNode]);
		for (int iFace = 0; iFace < adjacent_faces.size(); ++iFace)
		{
			int current_face_id = adjacent_faces[iFace];
			if (!is_face_visited[current_face_id])
			{
				int node1_id;
				int node2_id;
				int node3_id;

				tetrahedral_cushion_.GetNodesIdOfSurfaceFace(node1_id, node2_id, node3_id, adjacent_faces[iFace]);

				Eigen::Vector3d node1;
				Eigen::Vector3d node2;
				Eigen::Vector3d node3;
				deformed_tetrahedral_cushion_mesh_.GetVertex(node1_id, node1);
				deformed_tetrahedral_cushion_mesh_.GetVertex(node2_id, node2);
				deformed_tetrahedral_cushion_mesh_.GetVertex(node3_id, node3);

				//compute gradient
				ComputeTriangleAreaGradient(areas_all_adjacent_triangles_[current_face_id], node1, node2, node3,
					areas_all_adjacent_triangles_gradient_[current_face_id]);
				//cout << "face" << current_face_id << endl << areas_all_adjacent_triangles_gradient_[current_face_id] << endl;

				is_face_visited[current_face_id] = true;
			}

		}
	}
}

void LinearTetrahedralFemEvaluation::ComputeAreaOfSurfaceNodesGradient(const std::vector<std::vector<int>>* _faces_of_all_nodes)
{
	area_of_surface_nodes_gradient_.clear();
	area_of_surface_nodes_gradient_.resize(cushion_surface_nodes_.size());
	int num_all_variable = 3 * tetrahedral_cushion_.GetNodeNumber();// graident w.r.t. all nodes
	for (int iNode = 0; iNode < cushion_surface_nodes_.size(); ++iNode)
	{
		std::vector<Eigen::Triplet<double>> triplets;

		const std::vector<int>& adjacent_faces = _faces_of_all_nodes->at(cushion_surface_nodes_[iNode]);
		for (int iFace = 0; iFace < adjacent_faces.size(); ++iFace)
		{
			int current_face_id = adjacent_faces[iFace];

			// find the node id of current triangle
			int nodes_id[3] = { 0,0,0 };
			tetrahedral_cushion_.GetNodesIdOfSurfaceFace(nodes_id[0], nodes_id[1], nodes_id[2], current_face_id);

			for (int iNode = 0; iNode < 3; ++iNode)
			{
				for (int iCoordinate = 0; iCoordinate < 3; ++iCoordinate)
				{
					int row = 3 * nodes_id[iNode] + iCoordinate;
					int col = 0;
					double value = areas_all_adjacent_triangles_gradient_[current_face_id](iNode, iCoordinate);
					triplets.emplace_back(Eigen::Triplet<double>(row, col, value));
				}
			}
		}

		area_of_surface_nodes_gradient_[iNode].resize(num_all_variable, 1);
		area_of_surface_nodes_gradient_[iNode].setFromTriplets(triplets.begin(), triplets.end());
	}
}

void LinearTetrahedralFemEvaluation::FindGradient_SquaredForce_DeformedAllNodes()
{
	simulation_->Derivative_ElasticNodalForceSquaredNorm_DeformedAllNodes(cushion_surface_nodes_, squared_force_gradient_deformed_all_nodes_);
}

void LinearTetrahedralFemEvaluation::ComputeSigmoidSquaredForceGradient()
{
	sigmoid_squared_force_gradient_.resize(sigmoid_squared_force_.size());
	for (int i = 0; i < sigmoid_squared_force_.size(); ++i)
	{
		sigmoid_squared_force_gradient_[i] = scale_ * 0.5 * (1.0 - pow(2.0 * sigmoid_squared_force_[i] - 1.0, 2));
	}
}

void LinearTetrahedralFemEvaluation::ComputeGradientBeta_DeformedAllNodes()
{
	// Dbeta_i/Dx = DA_i * s_i + A_i * Ds_i * D|F_i|^2
	beta_gradient_deformed_all_nodes_.resize(cushion_surface_nodes_.size());
	for (int iNode = 0; iNode < cushion_surface_nodes_.size(); ++iNode)
	{
		beta_gradient_deformed_all_nodes_[iNode] =
			sigmoid_squared_force_[iNode] * area_of_surface_nodes_gradient_[iNode] // DA_i * s_i
			+ area_of_surface_nodes_[iNode] * sigmoid_squared_force_gradient_[iNode] * squared_force_gradient_deformed_all_nodes_[iNode];  //A_i * Ds_i * D | F_i | ^ 2
	}
}

void LinearTetrahedralFemEvaluation::ComputeGradientBeta_OriginalAllNodes(std::vector<Eigen::SparseMatrix<double>>& _beta_gradient_original_all_nodes)
{
	// Dbeta_i/Dx0 = DA_i/Dx0 * s_i + A_i * Ds_i * D|F_i|^2/Dx0
	// = -A_i * Ds_i * D|F_i|^2/Dx
	// (since DA_i/Dx0=0, D|F_i|^2/Dx0 = -D|F_i|^2/Dx)
	_beta_gradient_original_all_nodes.resize(cushion_surface_nodes_.size());
	for (int iNode = 0; iNode < cushion_surface_nodes_.size(); ++iNode)
	{
		_beta_gradient_original_all_nodes[iNode] =
			-area_of_surface_nodes_[iNode] * sigmoid_squared_force_gradient_[iNode] * squared_force_gradient_deformed_all_nodes_[iNode];  //-A_i * Ds_i * D|F_i|^2/Dx
	}
}

void LinearTetrahedralFemEvaluation::ComputeGradientGamma_DeformedAllNodes()
{
	// Dgamma_i/Dx = s_i*D|F_i|^2 * [s_i + 2*|F_i|^2*Ds_i]
	gamma_gradient_deformed_all_nodes_.resize(cushion_surface_nodes_.size());
	for (int iNode = 0; iNode < cushion_surface_nodes_.size(); ++iNode)
	{
		gamma_gradient_deformed_all_nodes_[iNode] =
			sigmoid_squared_force_[iNode] * squared_force_gradient_deformed_all_nodes_[iNode] //s_i*D|F_i|^2
			* (sigmoid_squared_force_[iNode] + 2.0 * cushion_surface_nodes_force_square_magnitude_[iNode] * sigmoid_squared_force_gradient_[iNode]);// s_i + 2*|F_i|^2*Ds_i
	}
}

void LinearTetrahedralFemEvaluation::ComputeGradientGamma_OriginalAllNodes(std::vector<Eigen::SparseMatrix<double>>& _gamma_gradient_original_all_nodes)
{
	_gamma_gradient_original_all_nodes.resize(gamma_gradient_deformed_all_nodes_.size());
	for (int iNode = 0; iNode < gamma_gradient_deformed_all_nodes_.size(); ++iNode)
	{
		_gamma_gradient_original_all_nodes[iNode] = -gamma_gradient_deformed_all_nodes_[iNode];
	}
}

void LinearTetrahedralFemEvaluation::ComputeAveragePressureGradient(Eigen::VectorXd& _gradient, Eigen::VectorXd& _gradient_sum_beta,
	const std::vector<Eigen::SparseMatrix<double>>& _gamma_gradient, 
	const std::vector<Eigen::SparseMatrix<double>>& _beta_gradient, 
	double _sum_beta, const std::vector<double>& _alpha, double _sum_alpha) const
{
	if (_gamma_gradient.size() > 0)
	{
		int variable_num = _gamma_gradient[0].rows();
		assert(_gamma_gradient.size() == _beta_gradient.size());
		assert(_gamma_gradient[0].rows() == _beta_gradient[0].rows());// We check only first element, while in fact, all gradient should have the same row

		// gradient = term1 + term2

		// find term2 and _gradient_sum_beta
		_gradient_sum_beta.setZero(variable_num);
		for (int iNode = 0; iNode < _beta_gradient.size(); ++iNode)
		{
			_gradient_sum_beta += _beta_gradient[iNode].col(0);
		}
		Eigen::VectorXd term2 = - _sum_alpha * _gradient_sum_beta / (_sum_beta * _sum_beta);

		// find term1
		Eigen::VectorXd term1;
		term1.setZero(variable_num);
		for (int iNode = 0; iNode < _beta_gradient.size(); ++iNode)
		{
			term1 += (_gamma_gradient[iNode].col(0) - _alpha[iNode] * _beta_gradient[iNode].col(0)) / (area_of_surface_nodes_[iNode] * sigmoid_squared_force_[iNode]);
		}
		term1 /= _sum_beta;

		_gradient = term1 + term2;
	}
}

void LinearTetrahedralFemEvaluation::ComputeAreaVarianceGradient(Eigen::VectorXd& _gradient, 
	const std::vector<Eigen::SparseMatrix<double>>& _beta_gradient, const Eigen::VectorXd& _gradient_sum_beta, double _sum_beta) const
{

	if (_beta_gradient.size() > 0)
	{

		int num_u = area_along_u_.size();
		int num_v_free = tetrahedral_cushion_.GetVNumber() - 1;
		_gradient.setZero(_beta_gradient[0].rows());
		double average_beta = _sum_beta / double(num_u);
		Eigen::VectorXd negative_averaged_gradient_sum_beta = -_gradient_sum_beta / double(num_u);

		for (int iU = 0; iU < num_u; ++iU)
		{
			// the gradient at iU is: 2/num_u*(beta_iu - beta/num_u)*(Dbeta_iu - Dbeta/num_u)
			Eigen::VectorXd gradient_u = negative_averaged_gradient_sum_beta;
			for (int iV = 0; iV < num_v_free; ++iV)
			{
				int free_node_id = iU * num_v_free + iV;
				gradient_u += _beta_gradient[free_node_id].col(0);
			}
			gradient_u *= (area_along_u_[iU] - average_beta);
			_gradient += gradient_u;
		}

		_gradient *= 2.0 / double(num_u);
	}
}

void LinearTetrahedralFemEvaluation::ComputeForceVarianceGradient(Eigen::VectorXd& _gradient,
	const std::vector<Eigen::SparseMatrix<double>>& _gradient_squared_force, double _sum_force) const
{
	if (_gradient_squared_force.size() > 0)
	{
		//compute D(|F_i|^2*s_i)
		std::vector<Eigen::SparseMatrix<double>> gradients_squared_force_multiply_sigmoid;
		ComputeGradientSquaredForceMultiplySigmoid(gradients_squared_force_multiply_sigmoid, _gradient_squared_force);

		// sum D(|F_i|^2*s_i) on v (along cross-sections)
		// and then on all nodes
		int num_u = squared_force_along_u_.size();
		int num_v_free = tetrahedral_cushion_.GetVNumber() - 1;
		std::vector<Eigen::SparseMatrix<double>> sum_v_force_gradients(num_u);
		for (int iU = 0; iU < num_u; ++iU)
		{
			// intialize with first v, then trasverse v starting with index 1
			sum_v_force_gradients[iU] = gradients_squared_force_multiply_sigmoid[iU * num_v_free];
			for (int iV = 1; iV < num_v_free; ++iV)
			{
				int free_node_id = iU * num_v_free + iV;
				sum_v_force_gradients[iU] += gradients_squared_force_multiply_sigmoid[free_node_id];
			}
		}
		Eigen::VectorXd sum_force_gradients;
		sum_force_gradients.setZero(_gradient_squared_force[0].rows());
		for (int iU = 0; iU < num_u; ++iU)
		{
			sum_force_gradients += sum_v_force_gradients[iU].col(0);
		}
		Eigen::VectorXd negative_averaged_sum_force_gradients = -sum_force_gradients / double(num_u);


		// compute gradient
		_gradient.setZero(_gradient_squared_force[0].rows());
		double average_force = _sum_force / double(num_u);
		for (int iU = 0; iU < num_u; ++iU)
		{
			// the gradient at iU is: 
			// 2/num_u*(force_iu - sum_force/num_u)*(sum_v_force_gradients[iU] - sum_force_gradients/num_u)
			Eigen::VectorXd gradient_u = negative_averaged_sum_force_gradients;
			gradient_u += sum_v_force_gradients[iU];
			gradient_u *= (squared_force_along_u_[iU] - average_force);
			_gradient += gradient_u;
		}
		_gradient *= 2.0 / double(num_u);
	}
}

void LinearTetrahedralFemEvaluation::ComputeGradientSquaredForceMultiplySigmoid(std::vector<Eigen::SparseMatrix<double>>& _gradients,
	const std::vector<Eigen::SparseMatrix<double>>& _sqared_force_gradient) const
{
	// D|F_i|^2 * [s_i + |F_i|^2*Ds_i]
	_gradients.resize(cushion_surface_nodes_.size());
	for (int iNode = 0; iNode < cushion_surface_nodes_.size(); ++iNode)
	{
		_gradients[iNode] =
			_sqared_force_gradient[iNode] //D|F_i|^2
			* (sigmoid_squared_force_[iNode] + cushion_surface_nodes_force_square_magnitude_[iNode] * sigmoid_squared_force_gradient_[iNode]);// s_i + 2*|F_i|^2*Ds_i
	}
}

void LinearTetrahedralFemEvaluation::Evaluate_SquaredForce(double& _sum_beta, double& _sum_force, double& _average_squared_pressure, double& _force_variance, double& _area_variance)
{
	// prepare information for evaluation
	simulation_->FindResultNodalForce();
	FindCushionSurfaceNodes();
	simulation_->FindForceManitudeAtGivenNodes(cushion_surface_nodes_, cushion_surface_nodes_force_magnitude_);
	ComputeSquaredVector(cushion_surface_nodes_force_square_magnitude_, cushion_surface_nodes_force_magnitude_);
	ComputeSigmoidSquaredForce();
	ComputeAreasOfAdjacentTriangles();

	// compute area distribution (also get the sum of node contact areas)
	ComputeAreaOfSurfaceNodes(cushion_surface_nodes_, tetrahedral_cushion_.GetFacesOfNodes());
	ComputeMetricAlongU(area_along_u_, area_of_surface_nodes_, sigmoid_squared_force_);
	Variance compute_variance;
	_area_variance = compute_variance.computeVariance(area_along_u_, _sum_beta);

	// compute force distribution (also get the sum of forces)
	ComputeMetricAlongU(squared_force_along_u_, cushion_surface_nodes_force_square_magnitude_, sigmoid_squared_force_);
	_force_variance = compute_variance.computeVariance(squared_force_along_u_, _sum_force);

	// compute average pressure
	sum_alpha_ = 0.0;
	alpha_.resize(cushion_surface_nodes_force_square_magnitude_.size());
	for (int iNode = 0; iNode < cushion_surface_nodes_force_square_magnitude_.size(); ++iNode)
	{
		alpha_[iNode] = cushion_surface_nodes_force_square_magnitude_[iNode] / area_of_surface_nodes_[iNode] * sigmoid_squared_force_[iNode];
		sum_alpha_ += alpha_[iNode];
	}
	_average_squared_pressure = sum_alpha_ / _sum_beta;
}

void LinearTetrahedralFemEvaluation::ComputeSigmoidAccurateForce()
{
	sigmoid_accurate_force_.resize(cushion_surface_nodes_force_square_magnitude_.size());
	for (int iNode = 0; iNode < cushion_surface_nodes_force_square_magnitude_.size(); ++iNode)
	{
		sigmoid_accurate_force_[iNode] = std::tanh(scale_ * cushion_surface_nodes_force_square_magnitude_[iNode]);
	}
}

void LinearTetrahedralFemEvaluation::ComputeSigmoidAccurateForceGradient()
{
	sigmoid_accurate_force_gradient_.resize(cushion_surface_nodes_force_square_magnitude_.size());
	for (int iNode = 0; iNode < cushion_surface_nodes_force_square_magnitude_.size(); ++iNode)
	{
		sigmoid_accurate_force_gradient_[iNode] = scale_ * (1.0 - pow(sigmoid_accurate_force_[iNode], 2.0));
	}
}

void LinearTetrahedralFemEvaluation::ComputeGradientBeta_DeformedAllNodes_type3()
{
	// Dbeta_i/Dx = DA_i * s_i + A_i * Ds_i * D|F_i|^2
	beta_gradient_deformed_all_nodes_.resize(cushion_surface_nodes_.size());
	for (int iNode = 0; iNode < cushion_surface_nodes_.size(); ++iNode)
	{
		beta_gradient_deformed_all_nodes_[iNode] =
			sigmoid_accurate_force_[iNode] * area_of_surface_nodes_gradient_[iNode] // DA_i * s_i
			+ area_of_surface_nodes_[iNode] * sigmoid_accurate_force_gradient_[iNode] * squared_force_gradient_deformed_all_nodes_[iNode];  //A_i * Ds_i * D | F_i | ^ 2
	}

}

void LinearTetrahedralFemEvaluation::ComputeGradientBeta_OriginalAllNodes_type3(std::vector<Eigen::SparseMatrix<double>>& _beta_gradient_original_all_nodes)
{
	// Dbeta_i/Dx0 = DA_i/Dx0 * s_i + A_i * Ds_i * D|F_i|^2/Dx0
	// = -A_i * Ds_i * D|F_i|^2/Dx
	// (since DA_i/Dx0=0, D|F_i|^2/Dx0 = -D|F_i|^2/Dx)
	_beta_gradient_original_all_nodes.resize(cushion_surface_nodes_.size());
	for (int iNode = 0; iNode < cushion_surface_nodes_.size(); ++iNode)
	{
		_beta_gradient_original_all_nodes[iNode] =
			-area_of_surface_nodes_[iNode] * sigmoid_accurate_force_gradient_[iNode] * squared_force_gradient_deformed_all_nodes_[iNode];  //-A_i * Ds_i * D|F_i|^2/Dx
	}

}

void LinearTetrahedralFemEvaluation::ComputeAveragePressureGradient_type3(Eigen::VectorXd& _gradient, 
	Eigen::VectorXd& _gradient_sum_beta, Eigen::VectorXd& _gradient_sum_force, 
	const std::vector<Eigen::SparseMatrix<double>>& _beta_gradient, const std::vector<Eigen::SparseMatrix<double>>& _squared_force_gradient, 
	double _sum_beta, double _average_pressure) const
{
	if (_beta_gradient.size() > 0)
	{
		int variable_num = _beta_gradient[0].rows();

		// gradient = term1 + term2

		// find term2 and _gradient_sum_beta
		_gradient_sum_beta.setZero(variable_num);
		for (int iNode = 0; iNode < _beta_gradient.size(); ++iNode)
		{
			_gradient_sum_beta += _beta_gradient[iNode].col(0);
		}
		Eigen::VectorXd term2 = -_average_pressure * 2.0 * _gradient_sum_beta / _sum_beta;

		// find term1 and _gradient_sum_force
		_gradient_sum_force.setZero(variable_num);
		for (int iNode = 0; iNode < cushion_surface_nodes_.size(); ++iNode)
		{
			_gradient_sum_force += _squared_force_gradient[iNode].col(0);
		}
		Eigen::VectorXd term1 = _gradient_sum_force / _sum_beta / _sum_beta;

		_gradient = term1 + term2;
	}
}

void LinearTetrahedralFemEvaluation::ComputeForceVarianceGradient_type3(Eigen::VectorXd& _gradient, 
	const std::vector<Eigen::SparseMatrix<double>>& _gradient_squared_force, const Eigen::VectorXd& _gradient_sum_force, double _sum_force) const
{
	// compute gradient
	int num_u = squared_force_along_u_.size();
	int num_v_free = tetrahedral_cushion_.GetVNumber() - 1;
	_gradient.setZero(_gradient_squared_force[0].rows());
	double average_force = _sum_force / double(num_u);
	Eigen::VectorXd negative_averaged_sum_force_gradients = -_gradient_sum_force / double(num_u);
	for (int iU = 0; iU < num_u; ++iU)
	{
		// the gradient at iU is: 
		// 2/num_u*(force_iu - sum_force/num_u)*(sum_v_force_gradients[iU] - sum_force_gradients/num_u)
		Eigen::VectorXd gradient_u = negative_averaged_sum_force_gradients;
		for (int iV = 0; iV < num_v_free; ++iV)
		{
			int free_node_id = iU * num_v_free + iV;
			gradient_u += _gradient_squared_force[free_node_id];
		}
		gradient_u *= (squared_force_along_u_[iU] - average_force);
		_gradient += gradient_u;
	}
	_gradient *= 2.0 / double(num_u);
}


void LinearTetrahedralFemEvaluation::TetrahedralizeCushion(int _num_layer, double _height)
{
	tetrahedral_cushion_.Tetrahedralize(cushion_, _num_layer, _height);	
	num_layer_ = _num_layer;
	height_ = _height;
}

void LinearTetrahedralFemEvaluation::ReTetrahedralizeCushion()
{
	tetrahedral_cushion_.Tetrahedralize(cushion_, num_layer_, height_);
}
