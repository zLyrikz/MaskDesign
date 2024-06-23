#include "FemManager.h"
#include "../Utility/EigenOperations.h"
#include<Eigen/SparseCholesky>	

#include <iostream>
#include <fstream>
#include <windows.h>


using std::cout;
using std::endl;

FemManager::FemManager()
{
	tetrahedral_mesh_ = nullptr;
}

FemManager::~FemManager()
{
}

void FemManager::SetTetrahedralMesh(const TetrahedralMesh* _mesh)
{
	tetrahedral_mesh_ = _mesh;
	elements_.SetBaseMesh(tetrahedral_mesh_);
	total_node_num_ = tetrahedral_mesh_->GetVertex()->rows();
}

void FemManager::SetMaterial(double _E, double _nu)
{
	elements_.SetMaterial(_E, _nu);
	E_ = _E;
	nu_ = _nu;
}

void FemManager::SetBoundaryCondition(const std::vector<std::pair<int, Eigen::Vector3d>>& _displacement_boundary, const std::vector<std::pair<int, Eigen::Vector3d>>& _force_boundary)
{
	// 4 jobs in this function:
	// 1
	elements_.SetBoundaryCondition(_displacement_boundary, _force_boundary);

	// 2
	displacement_boundary_ = _displacement_boundary;
	force_boundary_ = _force_boundary;

	// 3 
	// get unknown_displace_node_id_
	// find the unknown variable index
	std::vector<bool> is_unknown(total_node_num_, true);
	for (auto& iNode : _displacement_boundary)
	{
		is_unknown[iNode.first] = false;
	}
	// arrange the index in a vector
	unknown_displace_node_id_.clear();
	unknown_displace_node_id_.reserve(total_node_num_ - _displacement_boundary.size());
	for (int iNode = 0; iNode < total_node_num_; ++iNode)
	{
		if (is_unknown[iNode])
		{
			unknown_displace_node_id_.push_back(iNode);
		}
	}
	assert(unknown_displace_node_id_.size() == total_node_num_ - _displacement_boundary.size());

	// 4
	unknown_node_num_ = unknown_displace_node_id_.size();
}

void FemManager::AddUniformBodyForce(const Eigen::Vector3d& _body_force)
{
	elements_.SetUniformBodyForce(_body_force);
}

void FemManager::Compute()
{
	cout << "FEM start ... " << endl;
	LARGE_INTEGER t11, t12, tc1;
	QueryPerformanceFrequency(&tc1);
	QueryPerformanceCounter(&t11);

	// main procedure
	elements_.ConstructEquation();
	elements_.Solve();
	elements_.ComputeBoundaryForceLoad();


	QueryPerformanceCounter(&t12);

	cout << "FEM compute done; time spent " << (double)(t12.QuadPart - t11.QuadPart) / ((double)tc1.QuadPart) << "s" << endl;

	// get result
	all_displacement_ = *elements_.GetDisplacement();
}

void FemManager::SetCollisionInfo(const tmd::TriangleMeshDistance* _sdf_other_objects, double _collision_weight, int _sdf_function_index, const Discregrid::CubicLagrangeDiscreteGrid* _real_sdf)
{
	sdf_other_objects_ = _sdf_other_objects;
	collision_weight_ = _collision_weight;
	sdf_function_index_ = _sdf_function_index;
	real_sdf_ = _real_sdf;
}

void FemManager::Compute_CG_PenaltyCollision(const Eigen::Vector3d& _uniform_initial_displacement, int _max_iteration, double _epsilon_gradient, double _epsilon_fucntion, double _epsilon_step)
{
	//QuasiStaticFem solver;
	SolverSetup(_uniform_initial_displacement, solver_);

	// optimize
	cout << "CG optimization quasi-static FEM start ... " << endl;
	LARGE_INTEGER t11, t12, tc1;
	QueryPerformanceFrequency(&tc1);
	QueryPerformanceCounter(&t11);

	solver_.CgOptimize(_max_iteration, _epsilon_gradient, _epsilon_fucntion, _epsilon_step);

	QueryPerformanceCounter(&t12);
	cout << "FEM compute done; time spent " << (double)(t12.QuadPart - t11.QuadPart) / ((double)tc1.QuadPart) << "s" << endl;

	// result 
	Eigen::VectorXd displacement_variable_final;
	solver_.GetFinalDisplacement(displacement_variable_final);
	DisplacementVariable2AllDisplacement(displacement_variable_final, all_displacement_);
}

bool FemManager::Compute_LBFGS_PenaltyCollision(const Eigen::Vector3d& _uniform_initial_displacement, int _max_iteration, double _epsilon_gradient, double _epsilon_fucntion, double _epsilon_step)
{
	Eigen::VectorXd displacement_variable_final;
	bool is_good = Compute_LBFGS_PenaltyCollision(_uniform_initial_displacement, displacement_variable_final, _max_iteration, _epsilon_gradient, _epsilon_fucntion, _epsilon_step);
	if (is_good)
	{
		DisplacementVariable2AllDisplacement(displacement_variable_final, all_displacement_);
	}
	return is_good;
}

bool FemManager::Compute_LBFGS_PenaltyCollision(const TetrahedralMesh* _initial_displacement, int _max_iteration, double _epsilon_gradient, double _epsilon_fucntion, double _epsilon_step)
{
	//QuasiStaticFem solver;
	SolverSetup(_initial_displacement, solver_);

	// this can be used for checking if element orientation correct
	if (!elements_.IsVolumePositive())
	{
		cout << "ERROR FemManager::SolverSetup exists an element Volume not Positive" << endl;
		return false;
	}

	// optimize
	//cout << "L-BFGS optimization quasi-static FEM start ... " << endl;
	LARGE_INTEGER t11, t12, tc1;
	QueryPerformanceFrequency(&tc1);
	QueryPerformanceCounter(&t11);

	solver_.LbfgsOptimize(_max_iteration, _epsilon_gradient, _epsilon_fucntion, _epsilon_step, false);// false to not print iteration times and termination type

	QueryPerformanceCounter(&t12);
	//cout << "FEM compute done; time spent " << (double)(t12.QuadPart - t11.QuadPart) / ((double)tc1.QuadPart) << "s" << endl;

	// result 
	Eigen::VectorXd displacement_variable_final;
	solver_.GetFinalDisplacement(displacement_variable_final);
	DisplacementVariable2AllDisplacement(displacement_variable_final, all_displacement_);

	return true;
}

void FemManager::ComputePenaltyMethod(const Eigen::Vector3d& _uniform_initial_displacement, 
	int _max_outer_iteration, double _max_collision, double _penalty_increase, 
	int _max_inner_iteration, double _epsilon_gradient, double _epsilon_fucntion, double _epsilon_step, int _M,
	std::vector<Eigen::VectorXd>& _all_displacement_each_step)
{
	//QuasiStaticFem solver;
	SolverSetup(_uniform_initial_displacement, solver_);

	//// this can be used for checking if element orientation correct
	//if (!elements_.IsVolumePositive())
	//{
	//	cout << "ERROR FemManager::SolverSetup exists an element Volume not Positive" << endl;
	//	return;
	//}

	// optimize
	cout << "penalty method optimization (L-BFGS as inner solver) quasi-static FEM start ... " << endl;
	LARGE_INTEGER t11, t12, tc1;
	QueryPerformanceFrequency(&tc1);
	QueryPerformanceCounter(&t11);

	std::vector<Eigen::VectorXd> _unknown_displacement_each_step;
	solver_.LbfgsPenaltyOptimize(_max_outer_iteration, _max_collision, _penalty_increase,
		_max_inner_iteration, _epsilon_gradient, _epsilon_fucntion, _epsilon_step, _M, 
		_unknown_displacement_each_step);

	QueryPerformanceCounter(&t12);
	cout << "FEM compute done; time spent " << (double)(t12.QuadPart - t11.QuadPart) / ((double)tc1.QuadPart) << "s" << endl;

	// result 
	Eigen::VectorXd displacement_variable_final;
	solver_.GetFinalDisplacement(displacement_variable_final);
	DisplacementVariable2AllDisplacement(displacement_variable_final, all_displacement_);

	// result in each outer iteration
	for (auto& iDisplace : _unknown_displacement_each_step)
	{
		_all_displacement_each_step.emplace_back(Eigen::VectorXd());
		DisplacementVariable2AllDisplacement(iDisplace, _all_displacement_each_step.back());
	}
}

void FemManager::ReadResultDisplacement(std::string _file)
{
	std::ifstream file(_file);
	if (!file.is_open()) {
		std::cerr << "WARNING FemManager::ReadResultDisplacement can't open file!" << std::endl;
		return;
	}

	// get row number
	int num_rows = 0;
	std::string line;
	while (std::getline(file, line)) {
		++num_rows;
	}

	// return to begining
	file.clear();
	file.seekg(0, std::ios::beg);

	// read
	all_displacement_.resize(num_rows);
	for (int i = 0; i < num_rows; ++i) {
		std::getline(file, line);
		std::istringstream line_stream(line);
		double value;
		line_stream >> value;
		all_displacement_(i) = value;
	}
	file.close();

}

void FemManager::ReadBoundaryCondition(std::string _boundary_file, std::vector<std::pair<int, Eigen::Vector3d>>& _boundary_condition) const
{
	std::ifstream file(_boundary_file);
	if (!file.is_open()) {
		std::cerr << "FemManager::ReadBoundaryCondition cannot open file" << std::endl;
		return;
	}

	// read each row
	int intValue;
	double x, y, z;

	while (file >> intValue >> x >> y >> z) {
		Eigen::Vector3d vector3dValue(x, y, z);
		_boundary_condition.emplace_back(intValue, vector3dValue);
	}
	file.close();

}

bool FemManager::IsElementPositive() const
{
	if (!elements_.IsVolumePositive())
	{
		return false;
	}
	return true;
}

void FemManager::GetDeformedMesh(TetrahedralMesh& _deformed_mesh) const
{
	//_deformed_mesh = *tetrahedral_mesh_;
	//_deformed_mesh.TranslateIndividually(all_displacement_);

	//// print result info
	//ElementVolumeEnergyInfo_AllDisplace(all_displacement_);

	GetDeformedMesh(_deformed_mesh, all_displacement_);
}

void FemManager::GetDeformedMesh(TetrahedralMesh& _deformed_mesh, const Eigen::VectorXd& _all_displacement) const
{
	_deformed_mesh = *tetrahedral_mesh_;
	_deformed_mesh.TranslateIndividually(_all_displacement);

	// print result info
	//ElementVolumeEnergyInfo_AllDisplace(_all_displacement);
}

void FemManager::GetDeformedMesh_UnkownDisplacement(TetrahedralMesh& _deformed_mesh, const Eigen::VectorXd& _unknown_displacement) const
{
	Eigen::VectorXd all_displacement;
	DisplacementVariable2AllDisplacement(_unknown_displacement, all_displacement);
	GetDeformedMesh(_deformed_mesh, all_displacement);
}

void FemManager::GetBoundaryForceLineMesh(const TetrahedralMesh& _deformed_mesh, Mesh& _boundary_force, double _scale) const
{
	const std::vector<Eigen::Vector3d>* boundary_force = elements_.GetBoundaryForceLoad();
	const std::vector<std::pair<int, Eigen::Vector3d>>* displacement_boundary = elements_.GetDisplacementBoundary();
	GetBoundaryForceLineMesh(_deformed_mesh, _boundary_force, *boundary_force, *displacement_boundary, _scale);
	//const std::vector<Eigen::Vector3d>* boundary_force = elements_.GetBoundaryForceLoad();
	//const std::vector<std::pair<int, Eigen::Vector3d>>* displacement_boundary = elements_.GetDisplacementBoundary();
	//assert(boundary_force->size() == displacement_boundary->size());


	//for (int iBoundary = 0; iBoundary < displacement_boundary->size(); ++iBoundary)
	//{
	//	Eigen::Vector3d start_eigen;
	//	int vertex_id = displacement_boundary->at(iBoundary).first;
	//	_deformed_mesh.GetVertex(vertex_id, start_eigen);
	//	Eigen::Vector3d end_eigen = start_eigen + _scale * boundary_force->at(iBoundary);

	//	
	//	Mesh::Point start(start_eigen(0), start_eigen(1), start_eigen(2));
	//	Mesh::Point end(end_eigen(0), end_eigen(1), end_eigen(2));
	//	//Mesh::Point end = OpenMesh::vector_cast<Mesh::Point, Eigen::Vector3d>(end_eigen);

	//	//cout << "force magnitude = " << boundary_force->at(iBoundary).norm() << endl;

	//	Mesh::VertexHandle start_handle = _boundary_force.add_vertex(start);
	//	Mesh::VertexHandle end_handle = _boundary_force.add_vertex(end);
	//	Mesh::VertexHandle end2_handle = _boundary_force.add_vertex(end);

	//	_boundary_force.add_face(start_handle, end_handle, end2_handle);// a triangle face degenerated to a line segment
	//}
}

void FemManager::GetForceBoundaryForceLineMesh(const TetrahedralMesh& _deformed_mesh, Mesh& _boundary_force, double _scale) const
{
	std::vector<Eigen::Vector3d> boundary_force;
	elements_.ComputeForceBoundaryForceLoad(boundary_force);
	const std::vector<std::pair<int, Eigen::Vector3d>>* force_boundary = elements_.GetForceBoundary();
	GetBoundaryForceLineMesh(_deformed_mesh, _boundary_force, boundary_force, *force_boundary, _scale);
}

void FemManager::SaveResultXml(std::string _filepath, std::string _filename) const
{
	// save displacement
	std::string displacement_file = _filepath + _filename + "_displace.txt";
	OutputDisplacement(displacement_file);

	// save tetrahedral mesh
	std::string mesh_file = _filepath + _filename + "_mesh.mesh";
	OutputTetrahedralMesh(mesh_file);

	// save boundary condition
	std::string displace_boundary_file = _filepath + _filename + "_displaceBoundary.txt";;
	std::string force_boundary_file = _filepath + _filename + "_forceBoundary.txt";
	OutputBoundaryCondition(displace_boundary_file, force_boundary_file);

	// arrange info in xml
	std::string xml_file = _filepath + _filename + ".xml";
	// TODO implement a XmlParser
	//XmlParser::saveSiumationResult(xml_file, displacement_file, mesh_file, displace_boundary_file, force_boundary_file, nu_, E_);
}

void FemManager::OutputDisplacement(std::string _file) const
{
	std::ofstream outfile(_file);
	for (const auto& value : all_displacement_) {
		outfile << value << "\n";
	}
	outfile.close();

	cout << "FemManager::OutputDisplacement success" << endl;
}

void FemManager::OutputTetrahedralMesh(std::string _file) const
{
	tetrahedral_mesh_->WriteMesh(_file);
	cout << "FemManager::OutputTetrahedralMesh success" << endl;

}

void FemManager::OutputBoundaryCondition(std::string _displace_boundary_file, std::string _force_boundary_file) const
{
	std::ofstream outfile(_displace_boundary_file);
	for (const auto& value : displacement_boundary_)
	{
		outfile << value.first << " " << value.second(0) << " " << value.second(1) << " " << value.second(2) << endl;
	}
	outfile.close();

	outfile.open(_force_boundary_file);
	for (const auto& value : force_boundary_)
	{
		outfile << value.first << " " << value.second(0) << " " << value.second(1) << " " << value.second(2) << endl;
	}
	outfile.close();

	cout << "FemManager::OutputBoundaryCondition success" << endl;
}

void FemManager::ReadResultXml(std::string _file, TetrahedralMesh& _tet_mesh)
{
	// read
	std::string displacement_file;
	std::string mesh_file;
	std::string displace_boundary_file;
	std::string force_boundary_file;
	//// TODO implement a XmlParser
	//XmlParser::readSimulationResult(_file, displacement_file, mesh_file, displace_boundary_file, force_boundary_file, nu_, E_);

	// set 
	_tet_mesh.SetMesh(mesh_file);

	SetTetrahedralMesh(&_tet_mesh);
	SetMaterial(E_, nu_);

	ReadBoundaryCondition(displace_boundary_file, displacement_boundary_);
	ReadBoundaryCondition(force_boundary_file, force_boundary_);
	SetBoundaryCondition(displacement_boundary_, force_boundary_);

	elements_.ConstructEquation();
	ReadResultDisplacement(displacement_file);
}

void FemManager::FindResultNodalForce(const Eigen::VectorXd& _all_displacement, Eigen::VectorXd& _force) const
{
	const Eigen::SparseMatrix<double>* stiffness_matrix = elements_.GetStiffnessMatrix();
	_force = (*stiffness_matrix) * _all_displacement;
}

void FemManager::FindResultNodalForce()
{
	FindResultNodalForce(all_displacement_, nodal_force_);
}

void FemManager::FindForceManitudeAtGivenNodes(const std::vector<int>& _nodes, std::vector<double>& _node_force_magnitude) const
{
	_node_force_magnitude.resize(_nodes.size());
	for (int iNode = 0; iNode < _nodes.size(); ++iNode)
	{
		_node_force_magnitude[iNode] = sqrt(nodal_force_(3 * _nodes[iNode]) * nodal_force_(3 * _nodes[iNode]) +
			nodal_force_(3 * _nodes[iNode] + 1) * nodal_force_(3 * _nodes[iNode] + 1) +
			nodal_force_(3 * _nodes[iNode] + 2) * nodal_force_(3 * _nodes[iNode] + 2));
	}
}

void FemManager::FindResultSurfaceForceLineMeshAndAllForceMagnitude(const TetrahedralMesh& _deformed_mesh, Mesh& _force_mesh, std::vector<double>& _node_force_magnitude, double _scale) const
{
	const std::vector<int>* boundary_node = _deformed_mesh.GetBoundaryVertex();
	for (int iNode = 0; iNode < boundary_node->size(); ++iNode)
	{
		Eigen::Vector3d start_eigen;
		int vertex_id = boundary_node->at(iNode);
		_deformed_mesh.GetVertex(vertex_id, start_eigen);
		Eigen::Vector3d node_force(nodal_force_(3 * vertex_id), nodal_force_(3 * vertex_id + 1), nodal_force_(3 * vertex_id + 2));
		Eigen::Vector3d end_eigen = start_eigen + _scale * node_force;


		Mesh::Point start(start_eigen(0), start_eigen(1), start_eigen(2));
		Mesh::Point end(end_eigen(0), end_eigen(1), end_eigen(2));
		//Mesh::Point end = OpenMesh::vector_cast<Mesh::Point, Eigen::Vector3d>(end_eigen);

		//cout << "force magnitude = " << boundary_force->at(iBoundary).norm() << endl;

		Mesh::VertexHandle start_handle = _force_mesh.add_vertex(start);
		Mesh::VertexHandle end_handle = _force_mesh.add_vertex(end);
		Mesh::VertexHandle end2_handle = _force_mesh.add_vertex(end);

		_force_mesh.add_face(start_handle, end_handle, end2_handle);// a triangle face degenerated to a line segment
	}

	_node_force_magnitude.resize(total_node_num_);
	for (int iNode = 0; iNode < total_node_num_; ++iNode)
	{
		_node_force_magnitude[iNode] = sqrt(nodal_force_(3 * iNode) * nodal_force_(3 * iNode) +
			nodal_force_(3 * iNode + 1) * nodal_force_(3 * iNode + 1) +
			nodal_force_(3 * iNode + 2) * nodal_force_(3 * iNode + 2));
	}

}

void FemManager::FindNormalSurfaceForceLineMeshAndMagnitude(const std::vector<Mesh::Point>* _vertex_normal, const Mesh& _force_mesh, Mesh& _normal_force_mesh, std::vector<double>& _force_magnitude, double _scale) const
{
	int num_boundary_vertex = _vertex_normal->size();
	assert(_force_mesh.n_vertices() == 3 * num_boundary_vertex);
	_force_magnitude.resize(total_node_num_, 0);
	const std::vector<int>* boundary_vertex_id = tetrahedral_mesh_->GetBoundaryVertex();

	for (int iVertex = 0; iVertex < num_boundary_vertex; ++iVertex)
	{
		const Mesh::Point& force_start = _force_mesh.point(_force_mesh.vertex_handle(3 * iVertex));
		const Mesh::Point& force_end = _force_mesh.point(_force_mesh.vertex_handle(3 * iVertex + 1));
		Mesh::Point force_vector = (force_end - force_start) / _scale;
		Mesh::Point normal_force = force_vector.dot(_vertex_normal->at(iVertex)) * _vertex_normal->at(iVertex);
		//cout << "vertex normal" << _vertex_normal->at(iVertex) << " " << _vertex_normal->at(iVertex).norm() << " ";
		//cout << "force_vector" << force_vector << " " << force_vector.norm() << " ";
		//cout << "normal_force" << normal_force << " " << normal_force.norm() << "          ";

		_force_magnitude[boundary_vertex_id->at(iVertex)] = normal_force.norm();

		Mesh::Point normal_force_end =  force_start + _scale * normal_force;
		Mesh::VertexHandle start_handle = _normal_force_mesh.add_vertex(force_start);
		Mesh::VertexHandle end_handle = _normal_force_mesh.add_vertex(normal_force_end);
		Mesh::VertexHandle end2_handle = _normal_force_mesh.add_vertex(normal_force_end);
		_normal_force_mesh.add_face(start_handle, end_handle, end2_handle);// a triangle face degenerated to a line segment
	}
}

void FemManager::GetBoundaryConditionNodes(std::vector<int>& _nodes_id) const
{
	_nodes_id.clear();
	_nodes_id.resize(displacement_boundary_.size());
	for (int iNode = 0; iNode < _nodes_id.size(); ++iNode)
	{
		_nodes_id[iNode] = displacement_boundary_[iNode].first;
	}
}

void FemManager::FindDisplacementVariable(const Eigen::Vector3d& _uniform_initial_displacement, Eigen::VectorXd& _initial_displacement) const
{
	_initial_displacement.resize(3 * unknown_node_num_);
	for (int iNode = 0; iNode < unknown_node_num_; ++iNode)
	{
		for (int iCoordinate = 0; iCoordinate < 3; ++iCoordinate)
		{
			_initial_displacement[3 * iNode + iCoordinate] = _uniform_initial_displacement[iCoordinate];
		}
	}
}

void FemManager::DisplacementVariable2AllDisplacement(const Eigen::VectorXd& _unknown_displacement, Eigen::VectorXd& _all_displacement) const
{
	_all_displacement.resize(3 * total_node_num_);
	// the unknown
	for (int iNode = 0; iNode < unknown_node_num_; ++iNode)
	{
		for (int iCoordinate = 0; iCoordinate < 3; ++iCoordinate)
		{
			_all_displacement[3 * unknown_displace_node_id_[iNode] + iCoordinate] = _unknown_displacement[3 * iNode + iCoordinate];
		}
	}
	// the known
	for (const std::pair<int, Eigen::Vector3d>& iBoundary : displacement_boundary_)
	{
		for (int iCoordinate = 0; iCoordinate < 3; ++iCoordinate)
		{
			int id = 3 * iBoundary.first + iCoordinate;

			_all_displacement[id] = iBoundary.second(iCoordinate);
		}
	}
}

void FemManager::FindInitialNodesVariable(Eigen::VectorXd& _initial_nodes) const
{
	_initial_nodes.resize(3 * unknown_node_num_);
	for (int iNode = 0; iNode < unknown_node_num_; ++iNode)
	{
		for (int iCoordinate = 0; iCoordinate < 3; ++iCoordinate)
		{
			_initial_nodes[3 * iNode + iCoordinate] = tetrahedral_mesh_->GetVertex()->row(unknown_displace_node_id_[iNode])[iCoordinate];
		}
	}
}

double FemManager::ElasticEnergyFunction(const Eigen::VectorXd& _unknown_displacement, const Eigen::SparseMatrix<double>* _stiffness_matrix) const
{
	Eigen::VectorXd all_displacement;
	DisplacementVariable2AllDisplacement(_unknown_displacement, all_displacement);
	return 0.5 * all_displacement.transpose() * (*_stiffness_matrix) * all_displacement;
}

Eigen::VectorXd FemManager::ElasticEnergyGradient(const Eigen::VectorXd& _unknown_displacement, const Eigen::SparseMatrix<double>* _stiffness_matrix) const
{
	Eigen::VectorXd all_displacement;
	DisplacementVariable2AllDisplacement(_unknown_displacement, all_displacement);

	Eigen::VectorXd gradient;
	gradient.resize(3 * unknown_node_num_);

	Eigen::VectorXd full_gradient = (*_stiffness_matrix) * all_displacement;// this is the gradient if include the fixed displacement
	for (int iNode = 0; iNode < unknown_node_num_; ++iNode)
	{
		for (int iCoordinate = 0; iCoordinate < 3; ++iCoordinate)
		{
			gradient[3 * iNode + iCoordinate] = full_gradient(3 * unknown_displace_node_id_[iNode] + iCoordinate);
		}
	}
	return gradient;
}

double FemManager::ElasticEnergyValueAndGradient(const Eigen::VectorXd& _unknown_displacement, Eigen::VectorXd& _gradient, const Eigen::SparseMatrix<double>* _stiffness_matrix) const
{
	
	Eigen::VectorXd all_displacement;
	DisplacementVariable2AllDisplacement(_unknown_displacement, all_displacement);

	// gradient
	_gradient.resize(3 * unknown_node_num_);

	Eigen::VectorXd full_gradient = (*_stiffness_matrix) * all_displacement;// this is the gradient if include the fixed displacement
	for (int iNode = 0; iNode < unknown_node_num_; ++iNode)
	{
		for (int iCoordinate = 0; iCoordinate < 3; ++iCoordinate)
		{
			_gradient[3 * iNode + iCoordinate] = full_gradient(3 * unknown_displace_node_id_[iNode] + iCoordinate);
		}
	}

	// value
	//double energy = 0.5 * all_displacement.transpose() * (*_stiffness_matrix) * all_displacement;
	double energy = 0.5 * all_displacement.transpose() * full_gradient;

	return energy;
}

void FemManager::ElementVolumeEnergyInfo(const Eigen::VectorXd& _unknown_displacement) const
{
	Eigen::VectorXd all_displacement;
	DisplacementVariable2AllDisplacement(_unknown_displacement, all_displacement);
	ElementVolumeEnergyInfo_AllDisplace(all_displacement);
	//std::vector<double> volumes;
	//elements_.computeVolumesGivenDisplacement(all_displacement, volumes);
	//std::vector<double> energys;
	//elements_.computeElementEnergyGivenDisplacement(all_displacement, energys);

	//// volume info
	//int num_element = volumes.size();
	//int num_negative_volume = 0;
	//double average_negative_volume = 0.0;
	//for (int iElement = 0; iElement < num_element; ++iElement)
	//{
	//	if (volumes[iElement] < 0)
	//	{
	//		num_negative_volume += 1;
	//		average_negative_volume += volumes[iElement];
	//	}
	//}
	//average_negative_volume /= double(num_negative_volume);

	////energy info
	//double total_energy = 0.0;
	//double total_negative_volume_element_energy = 0.0;
	//double average_negative_volume_element_energy = 0.0;
	//for (int iElement = 0; iElement < num_element; ++iElement)
	//{
	//	total_energy += energys[iElement];
	//	if (volumes[iElement] < 0)
	//	{
	//		total_negative_volume_element_energy += energys[iElement];
	//	}
	//}
	//average_negative_volume_element_energy = total_negative_volume_element_energy / double(num_negative_volume);

	//cout << "total element num = " << num_element << "//// num_negative_volume =" << num_negative_volume << ",  average_negative_volume="<< average_negative_volume 
	//	<< "//// total_energy = "<< total_energy<< ", total_negative_volume_element_energy="<< total_negative_volume_element_energy<< ", average_negative_volume_element_energy=" << average_negative_volume_element_energy <<endl;
}

void FemManager::ElementVolumeEnergyInfo_AllDisplace(const Eigen::VectorXd& _all_displacement) const
{
	std::vector<double> volumes;
	elements_.computeVolumesGivenDisplacement(_all_displacement, volumes);
	std::vector<double> energys;
	elements_.computeElementEnergyGivenDisplacement(_all_displacement, energys);

	// volume info
	int num_element = volumes.size();
	int num_negative_volume = 0;
	double average_negative_volume = 0.0;
	for (int iElement = 0; iElement < num_element; ++iElement)
	{
		if (volumes[iElement] < 0)
		{
			num_negative_volume += 1;
			average_negative_volume += volumes[iElement];
		}
	}
	average_negative_volume /= double(num_negative_volume);

	//energy info
	double total_energy = 0.0;
	double total_negative_volume_element_energy = 0.0;
	double average_negative_volume_element_energy = 0.0;
	for (int iElement = 0; iElement < num_element; ++iElement)
	{
		total_energy += energys[iElement];
		if (volumes[iElement] < 0)
		{
			total_negative_volume_element_energy += energys[iElement];
		}
	}
	average_negative_volume_element_energy = total_negative_volume_element_energy / double(num_negative_volume);

	if (num_negative_volume > 0)
	{

		cout << "WARNING FROM FemManager, FEM result have inverted element\n total element num = " << num_element << "//// num_negative_volume =" << num_negative_volume << ",  average_negative_volume=" << average_negative_volume
			<< "//// total_energy = " << total_energy << ", total_negative_volume_element_energy=" << total_negative_volume_element_energy << ", average_negative_volume_element_energy=" << average_negative_volume_element_energy << endl;
	}
}

bool FemManager::Compute_LBFGS_PenaltyCollision(const Eigen::Vector3d& _uniform_initial_displacement, Eigen::VectorXd& _displacement_variable_final, int _max_iteration, double _epsilon_gradient, double _epsilon_fucntion, double _epsilon_step)
{
	//QuasiStaticFem solver;
	SolverSetup(_uniform_initial_displacement, solver_);

	// this can be used for checking if element orientation correct
	if (!elements_.IsVolumePositive())
	{
		cout << "ERROR FemManager::SolverSetup exists an element Volume not Positive" << endl;
		return false;
	}

	// optimize
	//cout << "L-BFGS optimization quasi-static FEM start ... " << endl;
	LARGE_INTEGER t11, t12, tc1;
	QueryPerformanceFrequency(&tc1);
	QueryPerformanceCounter(&t11);

	solver_.LbfgsOptimize(_max_iteration, _epsilon_gradient, _epsilon_fucntion, _epsilon_step, false);

	QueryPerformanceCounter(&t12);
	//cout << "FEM compute done; time spent " << (double)(t12.QuadPart - t11.QuadPart) / ((double)tc1.QuadPart) << "s" << endl;

	// result 
	solver_.GetFinalDisplacement(_displacement_variable_final);
	return true;
}

void FemManager::UpdateDerivativeInfo()
{
	FindUnknownVariablesGlobalId();
	FindStiffnessSubmatrixOnUnknownNodes();
	// find the row major stiffness matrix
	stiffness_matrix_row_major_ = *elements_.GetStiffnessMatrix();
	FindDerivativeOrderingMap();
}

const Eigen::SparseMatrix<double, Eigen::RowMajor>* FemManager::GetRowMajorStiffnessMatrix() const
{
	return &stiffness_matrix_row_major_;
}

void FemManager::Derivative_Deformed_Original_Term1(Eigen::SparseMatrix<double>& _matrix) const
{
	std::vector<Eigen::Triplet<double>> triplets;

	// the fixed nodes block (identity)
	int num_fix_dof = 3 * (total_node_num_ - unknown_node_num_);// degree of freedom of fixed nodes 
	for (int iFix = 0; iFix < num_fix_dof; ++iFix)
	{
		triplets.emplace_back(Eigen::Triplet<double>(iFix, iFix, 1.0));
	}

	// the free nodes block (-DF_free_Dx_free)
	Eigen::SparseMatrix<double> DF_free_Dx_free;
	Derivative_ForceOnFreeNodes_DeformedFreeNodes(DF_free_Dx_free);
	// Iterating over the nonzero coefficients of DF_free_Dx_free
	for (int iOuter = 0; iOuter < DF_free_Dx_free.outerSize(); ++iOuter)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(DF_free_Dx_free, iOuter); it; ++it)
		{
			// row and col += num_fix_dof because the elements are placed below the first block
			// value = -1.0 * it.value();
			triplets.emplace_back(Eigen::Triplet<double>(it.row() + num_fix_dof, it.col() + num_fix_dof, -it.value()));
		}
	}
	_matrix.resize(3 * total_node_num_, 3 * total_node_num_);
	_matrix.setFromTriplets(triplets.begin(), triplets.end());

}

void FemManager::Derivative_Deformed_Original_Term2(Eigen::SparseMatrix<double>& _matrix) const
{
	std::vector<Eigen::Triplet<double>> triplets;

	// top rows: fix nodes Dx_fix/Dx0
	int num_fix_node = displacement_boundary_.size();
	for (int iFix = 0; iFix < num_fix_node; ++iFix)
	{
		for (int iCoordinate = 0; iCoordinate < 3; ++iCoordinate)
		{
			triplets.emplace_back(Eigen::Triplet<double>(3 * iFix + iCoordinate, 3 * displacement_boundary_[iFix].first + iCoordinate, 1.0));
		}
	}

	// bottom rows: free nodes (DF_free/Dx0)
	Eigen::SparseMatrix<double> DF_free_Dx0;
	Derivative_ForceOnFreeNodes_OriginalAllNodes(DF_free_Dx0);
	// Iterating over the nonzero coefficients of DF_free/Dx0
	int num_fix_dof = 3 * num_fix_node;// degree of freedom of fixed nodes 
	for (int iOuter = 0; iOuter < DF_free_Dx0.outerSize(); ++iOuter)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(DF_free_Dx0, iOuter); it; ++it)
		{
			// row += num_fix_dof because the elements are placed at the bottom
			// col index is the same
			triplets.emplace_back(Eigen::Triplet<double>(it.row() + num_fix_dof, it.col(), it.value()));
		}
	}
	_matrix.resize(3 * total_node_num_, 3 * total_node_num_);
	_matrix.setFromTriplets(triplets.begin(), triplets.end());
}

const std::vector<int>& FemManager::GetDerivativeVariableRowOrderingMap() const
{
	return derivative_ordering_map_;
}

void FemManager::Derivative_ElasticNodalForceSquaredNorm_DeformedAllNodes(const std::vector<int>& _nodes, std::vector<Eigen::SparseMatrix<double>>& _derivatives) const
{
	// D|F|^2/Dx = 2*F^T*DF/Dx = 2* (-F^T) * (-DF/Dx)
	// -DF/Dx takes the rows from K
	// 
	_derivatives.resize(_nodes.size());
	for (int iNode = 0; iNode < _nodes.size(); ++iNode)
	{
		// at current node, find derivative: 2* (-F^T) * (-DF/Dx)
		std::vector<Eigen::Triplet<double>> triplets;

		// iterate over 3 elements of current nodal force F
		for (int iCoordinate = 0; iCoordinate < 3; ++iCoordinate)
		{
			// current F element global index
			// (-DF/Dx) at this element i.e. this row of stiffness_matrix K
			int row = 3 * _nodes[iNode] + iCoordinate;

			// mutipliy (-DF/Dx) with (-F^T) 
			// then assign value to triplets
			for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(stiffness_matrix_row_major_, row); it; ++it)
			{
				// it.index() i.e. inner index, here it is equal to it.col()
				// note the constructed sparse matrix is column major, so put values on corresponding rows
				triplets.emplace_back(Eigen::Triplet(it.index(), 0, nodal_force_(row) * it.value()));

			}
			

		}

		_derivatives[iNode].resize(3 * total_node_num_, 1);
		_derivatives[iNode].setFromTriplets(triplets.begin(), triplets.end());
		_derivatives[iNode] *= 2.0;
	}
}

void FemManager::GetBoundaryForceLineMesh(const TetrahedralMesh& _deformed_mesh, Mesh& _boundary_force_mesh,
	const std::vector<Eigen::Vector3d>& _boundary_force, const std::vector<std::pair<int, Eigen::Vector3d>>& _boundary_id, double _scale) const
{
	assert(_boundary_force.size() == _boundary_id.size());
	for (int iBoundary = 0; iBoundary < _boundary_id.size(); ++iBoundary)
	{
		Eigen::Vector3d start_eigen;
		int vertex_id = _boundary_id.at(iBoundary).first;
		_deformed_mesh.GetVertex(vertex_id, start_eigen);
		Eigen::Vector3d end_eigen = start_eigen + _scale * _boundary_force.at(iBoundary);


		Mesh::Point start(start_eigen(0), start_eigen(1), start_eigen(2));
		Mesh::Point end(end_eigen(0), end_eigen(1), end_eigen(2));
		//Mesh::Point end = OpenMesh::vector_cast<Mesh::Point, Eigen::Vector3d>(end_eigen);

		//cout << "force magnitude = " << boundary_force->at(iBoundary).norm() << endl;

		Mesh::VertexHandle start_handle = _boundary_force_mesh.add_vertex(start);
		Mesh::VertexHandle end_handle = _boundary_force_mesh.add_vertex(end);
		Mesh::VertexHandle end2_handle = _boundary_force_mesh.add_vertex(end);

		_boundary_force_mesh.add_face(start_handle, end_handle, end2_handle);// a triangle face degenerated to a line segment
	}
}

void FemManager::SolverSetup(const Eigen::Vector3d& _uniform_initial_displacement, QuasiStaticFem& _solver)
{	
	Eigen::VectorXd initial_displacement;
	FindDisplacementVariable(_uniform_initial_displacement, initial_displacement);
	SolverSetup(initial_displacement, _solver);


	//Eigen::VectorXd initial_nodes;
	//FindInitialNodesVariable(initial_nodes);


	//// initial setting
	//_solver.SetInitialConfiguration(initial_nodes);
	//_solver.SetInitialDisplacement(initial_displacement);
	//_solver.SetCollisionSdf(sdf_other_objects_);
	//_solver.SetCollisionWeight(collision_weight_);

	//// set energy function and gradient
	//// getting the stiffness matrix
	//elements_.FindStiffnessMatrix();

	//const Eigen::SparseMatrix<double>* stiffness_matrix = elements_.GetStiffnessMatrix();
	//// TODO pre-separate the stiffness matrix, remove the constant term where only known displacement involved
	//// TODO test function evaluation efficiency! (sparse matrix multiplication)

	//*(_solver.GetElasticEnergyValueAndGradient()) = std::bind(&FemManager::ElasticEnergyValueAndGradient, this, std::placeholders::_1, std::placeholders::_2, stiffness_matrix);

	//// for debug and get more info in the optimization process
	////*(_solver.GetInfoFunction()) = std::bind(&FemManager::ElementVolumeEnergyInfo, this, std::placeholders::_1);
}

void FemManager::SolverSetup(const TetrahedralMesh* _initial_displacement, QuasiStaticFem& _solver)
{
	Eigen::VectorXd initial_displacement;

	initial_displacement.resize(3 * unknown_node_num_);
	for (int iNode = 0; iNode < unknown_node_num_; ++iNode)
	{
		for (int iCoordinate = 0; iCoordinate < 3; ++iCoordinate)
		{
			initial_displacement[3 * iNode + iCoordinate] = 
				_initial_displacement->GetVertex()->row(unknown_displace_node_id_[iNode])[iCoordinate] - 
				tetrahedral_mesh_->GetVertex()->row(unknown_displace_node_id_[iNode])[iCoordinate];

			//cout << initial_displacement[3 * iNode + iCoordinate] << " ";
		}
	}

	SolverSetup(initial_displacement, _solver);

}

void FemManager::SolverSetup(const Eigen::VectorXd& _initial_displacement, QuasiStaticFem& _solver)
{
	Eigen::VectorXd initial_nodes;
	FindInitialNodesVariable(initial_nodes);

	// initial setting
	_solver.SetInitialConfiguration(initial_nodes);
	_solver.SetInitialDisplacement(_initial_displacement);
	_solver.SetCollisionSdf(sdf_other_objects_, sdf_function_index_, real_sdf_);
	_solver.SetCollisionWeight(collision_weight_);

	// set energy function and gradient
	// getting the stiffness matrix
	elements_.FindStiffnessMatrix();

	const Eigen::SparseMatrix<double>* stiffness_matrix = elements_.GetStiffnessMatrix();
	// TODO pre-separate the stiffness matrix, remove the constant term where only known displacement involved
	// TODO test function evaluation efficiency! (sparse matrix multiplication)

	*(_solver.GetElasticEnergyValueAndGradient()) = std::bind(&FemManager::ElasticEnergyValueAndGradient, this, std::placeholders::_1, std::placeholders::_2, stiffness_matrix);

	// for debug and get more info in the optimization process
	//*(_solver.GetInfoFunction()) = std::bind(&FemManager::ElementVolumeEnergyInfo, this, std::placeholders::_1);

}

void FemManager::FindUnknownVariablesGlobalId()
{
	unknown_variables_global_id_.resize(3 * unknown_node_num_);
	for (int iNode = 0; iNode < unknown_displace_node_id_.size(); ++iNode)
	{
		unknown_variables_global_id_[3 * iNode] = 3 * unknown_displace_node_id_[iNode];
		unknown_variables_global_id_[3 * iNode + 1] = 3 * unknown_displace_node_id_[iNode] + 1;
		unknown_variables_global_id_[3 * iNode + 2] = 3 * unknown_displace_node_id_[iNode] + 2;
	}
}

void FemManager::FindStiffnessSubmatrixOnUnknownNodes()
{
	const Eigen::SparseMatrix<double>* stiffness_matrix = elements_.GetStiffnessMatrix();// K
	stiffness_submatrix_unknown_nodes_ = Eigen::SparseMatrix<double>();
	EigenOperations::ExtractSubmatrix_InputColumnMajor(*stiffness_matrix,
		stiffness_submatrix_unknown_nodes_, unknown_variables_global_id_, unknown_variables_global_id_);
}

void FemManager::Derivative_ForceOnFreeNodes_DeformedFreeNodes(Eigen::SparseMatrix<double>& _jacobian) const
{
	// DF/Dx = -D^2E/Dx^2
	// E = elastic energy + collision energy

	// D^2E_elastic = K_Nf (submatrix of stiffness matrix K, at the free nodes)
	//FindStiffnessSubmatrixOnUnknownNodes();// should call this function beforehand

	// D^2E collision
	Eigen::SparseMatrix<double> hessian_collision;
	solver_.FinalCollisionHessian(hessian_collision);

	_jacobian = -1.0 * (stiffness_submatrix_unknown_nodes_ + hessian_collision);


}

void FemManager::Derivative_ForceOnFreeNodes_OriginalAllNodes(Eigen::SparseMatrix<double>& _jacobian) const
{
	const Eigen::SparseMatrix<double>* stiffness_matrix = elements_.GetStiffnessMatrix();// K

	std::vector<int> all_column_index(stiffness_matrix->cols());
	std::iota(all_column_index.begin(), all_column_index.end(), 0);// all_column_index[i] = i
	EigenOperations::ExtractSubmatrix_InputColumnMajor(*stiffness_matrix,
		_jacobian, unknown_variables_global_id_, all_column_index);
}

void FemManager::FindDerivativeOrderingMap()
{
	derivative_ordering_map_.resize(3 * total_node_num_);
	for (int iFix = 0; iFix < displacement_boundary_.size(); ++iFix)
	{
		for (int iCoordinate = 0; iCoordinate < 3; ++iCoordinate)
		{
			int current_variable_id = 3 * displacement_boundary_[iFix].first + iCoordinate;

			// map the current node variable index to front rows(its index in the boundary condition)
			derivative_ordering_map_[current_variable_id] = 3 * iFix + iCoordinate;
		}
	}

	int num_fixed_variables = 3 * displacement_boundary_.size();
	for (int iFree = 0; iFree < unknown_variables_global_id_.size(); ++iFree)
	{
		// map the free nodes to the back of the fix nodes
		derivative_ordering_map_[unknown_variables_global_id_[iFree]] = iFree + num_fixed_variables;
	}
}


Eigen::VectorXd FemManager::SimulationFunction(const Eigen::VectorXd& _x0)
{
	// 
	static TetrahedralMesh tetrahedral_mesh = *tetrahedral_mesh_;
	tetrahedral_mesh.ChangeVertex(_x0);


	SetTetrahedralMesh(&tetrahedral_mesh);


	static TetrahedralMesh deformed_tetrahedral_cushion_mesh;
	bool is_FEM_good = true;
	if (deformed_tetrahedral_cushion_mesh.GetVertex()->rows() == 0)
	{
		// initial displacement is set as zero
		is_FEM_good = Compute_LBFGS_PenaltyCollision(Eigen::Vector3d::Zero(), 5000, 1e-10/*1*/, 0, 0);
	}
	else
	{
		// initial displacement set as the last simulation result configuration
		is_FEM_good = Compute_LBFGS_PenaltyCollision(&deformed_tetrahedral_cushion_mesh, 5000, 1e-10/*1*/, 0, 0);
	}

	if (is_FEM_good == false)
	{
		cout << "NOTICE from LinearTetrahedralFemEvaluation::EvaluationPreparation, exists an element Volume not Positive" << endl;
	}
	else
	{
		// get the deformed cushion mesh of this simulation
		GetDeformedMesh(deformed_tetrahedral_cushion_mesh);
	}

	return _x0 + all_displacement_;
}

Eigen::Vector3d FemManager::SimulationFunction_OneNode(const Eigen::VectorXd& _x0, int _node_id)
{
	Eigen::VectorXd x_all = SimulationFunction(_x0);
	Eigen::Vector3d x_i(x_all(3 * _node_id + 0), x_all(3 * _node_id + 1), x_all(3 * _node_id + 2));

	return x_i;
}

Eigen::VectorXd FemManager::TotalForceOnFreeNodesFunction(const Eigen::VectorXd& _deformed_free_nodes)
{
	Eigen::VectorXd unknown_displacement(_deformed_free_nodes.rows());
	for (int iFree = 0; iFree < unknown_node_num_; ++iFree)
	{
		unknown_displacement(3 * iFree + 0) = _deformed_free_nodes(3 * iFree + 0) - (*tetrahedral_mesh_->GetVertex())(unknown_displace_node_id_[iFree], 0);
		unknown_displacement(3 * iFree + 1) = _deformed_free_nodes(3 * iFree + 1) - (*tetrahedral_mesh_->GetVertex())(unknown_displace_node_id_[iFree], 1);
		unknown_displacement(3 * iFree + 2) = _deformed_free_nodes(3 * iFree + 2) - (*tetrahedral_mesh_->GetVertex())(unknown_displace_node_id_[iFree], 2);
	}
	Eigen::VectorXd all_displacement;
	DisplacementVariable2AllDisplacement(unknown_displacement, all_displacement);

	const Eigen::SparseMatrix<double>* stiffness_matrix = elements_.GetStiffnessMatrix();
	Eigen::VectorXd elastic_force = -(*stiffness_matrix) * all_displacement;
	Eigen::VectorXd elastic_force_free_nodes(_deformed_free_nodes.rows());
	for (int i = 0; i < _deformed_free_nodes.rows(); ++i)
	{
		elastic_force_free_nodes[i] = elastic_force[unknown_variables_global_id_[i]];
	}
	double collision_energy = 0.0;
	Eigen::VectorXd collision_gradient;//gradient w.r.t. displacement <=> w.r.t. deformed nodes
	solver_.CollisionValueAndGradient(unknown_displacement, collision_energy, collision_gradient);

	return elastic_force_free_nodes - collision_weight_ * collision_gradient;
}



Eigen::VectorXd FemManager::TotalForceOnFreeNodesFunction_OnOriginalNodes(const Eigen::VectorXd& _original_nodes)
{
	Eigen::VectorXd deformed_nodes(3 * total_node_num_);
	for (int iNode = 0; iNode < total_node_num_; ++iNode)
	{
		for (int iCoordinate = 0; iCoordinate < 3; ++iCoordinate)
		{
			deformed_nodes[3 * iNode + iCoordinate] = all_displacement_[3 * iNode + iCoordinate] + tetrahedral_mesh_->GetVertex()->row(iNode)[iCoordinate];
		}
	}
	Eigen::VectorXd all_displacement(3 * total_node_num_);//current ones
	for (int iNode = 0; iNode < total_node_num_; ++iNode)
	{
		for (int iCoordinate = 0; iCoordinate < 3; ++iCoordinate)
		{
			all_displacement[3 * iNode + iCoordinate] = deformed_nodes[3 * iNode + iCoordinate] - _original_nodes[3 * iNode + iCoordinate];
		}
	}

	const Eigen::SparseMatrix<double>* stiffness_matrix = elements_.GetStiffnessMatrix();
	Eigen::VectorXd elastic_force = -(*stiffness_matrix) * all_displacement;
	Eigen::VectorXd elastic_force_free_nodes(3*unknown_node_num_);
	for (int i = 0; i < 3 * unknown_node_num_; ++i)
	{
		elastic_force_free_nodes[i] = elastic_force[unknown_variables_global_id_[i]];
	}

	return elastic_force_free_nodes;
}
