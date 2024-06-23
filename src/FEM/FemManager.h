#pragma once
#include "QuasiStaticFem.h"
#include "TetrahedralElements.h"

#include "../MeshProcessing/TetrahedralMesh.h"
#include "../MeshViewer/MeshDefinition.h"
#include <Discregrid/All>
#include <string>

// this should be called the [Tetrahedral Mesh] FEM manager
// for solving static FEM on tetrahedral mesh with linear element
// note: mind the unit! normally meters are used
// expanded to support quasi-static FEM with contact
class FemManager
{
public:
	FemManager();
	~FemManager();

	// basic settings
	// call "set" functions in the same order 
	void SetTetrahedralMesh(const TetrahedralMesh* _mesh);
	void SetMaterial(double _E, double _nu);
	// input pairs: node index and displacement 
	// note if _force_boundary misses some vertex, then they are regarded as 0 force applied
	void SetBoundaryCondition(const std::vector<std::pair<int, Eigen::Vector3d>>& _displacement_boundary,
		const std::vector<std::pair<int, Eigen::Vector3d>>& _force_boundary = std::vector<std::pair<int, Eigen::Vector3d>>());
	void AddUniformBodyForce(const Eigen::Vector3d& _body_force);

	// this is the most basic one: linear tet mesh, single body, no contact, specify boundary condition and solve
	void Compute();

	// quasi-static FEM with contact
	// need some more settings
	void SetCollisionInfo(const tmd::TriangleMeshDistance* _sdf_other_objects, double _collision_weight, int _sdf_function_index = 0, const Discregrid::CubicLagrangeDiscreteGrid* _real_sdf = nullptr);
	// set a uniform initial displacement on the unknown variables
	void Compute_CG_PenaltyCollision(const Eigen::Vector3d& _uniform_initial_displacement, int _max_iteration = 5000, double _epsilon_gradient = 1e-5, double _epsilon_fucntion = 0, double _epsilon_step = 0.1);
	bool Compute_LBFGS_PenaltyCollision(const Eigen::Vector3d& _uniform_initial_displacement, int _max_iteration = 5000, double _epsilon_gradient = 1e-5, double _epsilon_fucntion = 0, double _epsilon_step = 0.1);
	// set initial displacement as a reference mesh configuration (exactly same topolgy, different vertex position)
	bool Compute_LBFGS_PenaltyCollision(const TetrahedralMesh* _initial_displacement, int _max_iteration = 5000, double _epsilon_gradient = 1e-5, double _epsilon_fucntion = 0, double _epsilon_step = 0.1);
	// the penalty term increases after each LBFGS iteration
	void ComputePenaltyMethod(const Eigen::Vector3d& _uniform_initial_displacement, 
		int _max_outer_iteration, double _max_collision, double _penalty_increase, 
		int _max_inner_iteration, double _epsilon_gradient = 1e-5, double _epsilon_fucntion = 0, double _epsilon_step = 0.1, int _M = 5,
		std::vector<Eigen::VectorXd>& _all_displacement_each_step = std::vector<Eigen::VectorXd>());

	void GetDeformedMesh(TetrahedralMesh& _deformed_mesh) const;
	void GetDeformedMesh(TetrahedralMesh& _deformed_mesh, const Eigen::VectorXd& _all_displacement) const;
	void GetDeformedMesh_UnkownDisplacement(TetrahedralMesh& _deformed_mesh, const Eigen::VectorXd& _unknown_displacement) const;

	// the static FEM result info
	// find the node load on the displacement boundary
	void GetBoundaryForceLineMesh(const TetrahedralMesh& _deformed_mesh, Mesh& _boundary_force, double _scale = 1.0) const;
	//find the node load on the force boundary
	void GetForceBoundaryForceLineMesh(const TetrahedralMesh& _deformed_mesh, Mesh& _boundary_force, double _scale = 1.0) const;

	// quasi-static FEM
	void FindResultNodalForce(const Eigen::VectorXd& _all_displacement, Eigen::VectorXd& _force) const;
	void FindResultNodalForce();
	void FindForceManitudeAtGivenNodes(const std::vector<int>& _nodes, std::vector<double>& _node_force_magnitude) const;
	void FindResultSurfaceForceLineMeshAndAllForceMagnitude(const TetrahedralMesh& _deformed_mesh, Mesh& _force_mesh, std::vector<double>& _node_force_magnitude, double _scale = 1.0) const;
	// the magnitude includes only boundary nodes
	// _force_mesh is the result from FindResultSurfaceForceLineMeshAndAllForceMagnitude()
	// notice the vertex index in the surface mesh corresponds to boundary_vertex member in TetrahedralMesh
	// mind the direction of _vertex_normal and force vector
	// note the scale is for getting the correct force from _force_mesh(since it may be scaled in the mesh)
	void FindNormalSurfaceForceLineMeshAndMagnitude(const std::vector<Mesh::Point>* _vertex_normal, const Mesh& _force_mesh, Mesh& _normal_force_mesh, std::vector<double>& _force_magnitude, double _scale = 1.0) const;
	void GetBoundaryConditionNodes(std::vector<int>& _nodes_id) const;

	// save & read quasi-satic simulation result
	// 
	// full file path = _filepath + _filename + ".xml"
	// will also save dislacement; original tetrahedral mesh .mesh file; boundary condition
	void SaveResultXml(std::string _filepath, std::string _filename) const;
	void OutputDisplacement(std::string _file) const;
	void OutputTetrahedralMesh(std::string _file) const;
	void OutputBoundaryCondition(std::string _displace_boundary_file, std::string _force_boundary_file) const;

	void ReadResultXml(std::string _file, TetrahedralMesh& _tet_mesh);// output an original tetrahedral mesh, other info remains in FemManager
	void ReadResultDisplacement(std::string _file);// read displacement from the previous calculation
	void ReadBoundaryCondition(std::string _boundary_file, std::vector<std::pair<int, Eigen::Vector3d>>& _boundary_condition) const;

	bool IsElementPositive() const;

	// provide derivative info
	// call this function before calling following derivative calculation functions
	void UpdateDerivativeInfo();
	const Eigen::SparseMatrix<double, Eigen::RowMajor>* GetRowMajorStiffnessMatrix() const;

	// denote x as deformed, x0 as original nodes(natural configuration)
	// expression with a nice index ordering of x (fixed nodes at top rows, free nodes at bottom rows):
	//          [I        0        ]^(-1)	   [Dx_fix/Dx0 ]
	// Dx/Dx0 = [                  ]       *   [           ]
	//          [0 -DF_free/Dx_free]           [DF_free/Dx0]
	//         =:        term1      ^(-1)  *       term2
	// size = 3*all_node_number
	// term1 is positive-definite
	void Derivative_Deformed_Original_Term1(Eigen::SparseMatrix<double>& _matrix) const;
	void Derivative_Deformed_Original_Term2(Eigen::SparseMatrix<double>& _matrix) const;
	// return the mapping vector that gives the "nice index ordering"
	// derivative_ordering_map_[i] is the row index of variable i(global index) in the term1 and term2 matrix
	// e.g. for vector v, construct a permuted vector permuted_v like this: permuted_v[derivative_ordering_map_[i]] = v[i]
	const std::vector<int>& GetDerivativeVariableRowOrderingMap() const;

	// for each node in a given set, find the D|F|^2/Dx, where F is elastic force at a node(= - DE_elastic/Dx = - K(x-x0))
	// the derivative is sparse, stored in SparseMatrix with only one column
	// require nodal_force_ already computed
	void Derivative_ElasticNodalForceSquaredNorm_DeformedAllNodes(const std::vector<int>& _nodes, std::vector<Eigen::SparseMatrix<double>>& _derivatives) const;

	// each gradient row= 3*4, w.r.t.  its four vertices
	//void GetElementVolumesAndGradient(std::vector<double>& _volumes, std::vector<Eigen::VectorXd>& _gradient);

	// only for derivative debug
	//void TestDxDx0();
	//void TestDx_freeDx0();
	Eigen::VectorXd SimulationFunction(const Eigen::VectorXd& _x0);
	Eigen::Vector3d SimulationFunction_OneNode(const Eigen::VectorXd& _x0, int _node_id);
	Eigen::VectorXd TotalForceOnFreeNodesFunction(const Eigen::VectorXd& _deformed_free_nodes);
	//void TestDerivative_ForceOnFreeNodes_DeformedFreeNodes();
	//void TestDerivative_ForceOnFreeNodes_OriginalNodes();
	Eigen::VectorXd TotalForceOnFreeNodesFunction_OnOriginalNodes(const Eigen::VectorXd& _original_nodes);


private:
	// for setting the quasi-static solver
	// arrange the unknown displacement variables into one vector
	void FindDisplacementVariable(const Eigen::Vector3d& _uniform_initial_displacement, Eigen::VectorXd& _initial_displacement) const;
	// combine the input unknown and those in the boundary condition displacement
	void DisplacementVariable2AllDisplacement(const Eigen::VectorXd& _unknown_displacement, Eigen::VectorXd& _all_displacement) const;
	void FindInitialNodesVariable(Eigen::VectorXd& _initial_nodes) const;

	// callback functions during optimization
	double ElasticEnergyFunction(const Eigen::VectorXd& _unknown_displacement, const Eigen::SparseMatrix<double>* _stiffness_matrix) const;
	Eigen::VectorXd ElasticEnergyGradient(const Eigen::VectorXd& _unknown_displacement,const Eigen::SparseMatrix<double>* _stiffness_matrix) const;
	double ElasticEnergyValueAndGradient(const Eigen::VectorXd& _unknown_displacement, Eigen::VectorXd& _gradient, const Eigen::SparseMatrix<double>* _stiffness_matrix) const;
	void ElementVolumeEnergyInfo(const Eigen::VectorXd& _unknown_displacement) const;
	void ElementVolumeEnergyInfo_AllDisplace(const Eigen::VectorXd& _all_displacement) const;

	bool Compute_LBFGS_PenaltyCollision(const Eigen::Vector3d& _uniform_initial_displacement, Eigen::VectorXd& _displacement_variable_final, int _max_iteration = 5000, double _epsilon_gradient = 1e-5, double _epsilon_fucntion = 0, double _epsilon_step = 0.1);




private:
	// utility functions
	void GetBoundaryForceLineMesh(const TetrahedralMesh& _deformed_mesh, Mesh& _boundary_force_mesh, 
		const std::vector<Eigen::Vector3d>& _boundary_force, const std::vector<std::pair<int, Eigen::Vector3d>>& _boundary_id, double _scale = 1.0) const;

	// quasi-static FEM with contact
	// CG and L-BFGS optimization setup
	void SolverSetup(const Eigen::Vector3d& _uniform_initial_displacement, QuasiStaticFem& _solver);
	void SolverSetup(const TetrahedralMesh* _initial_displacement, QuasiStaticFem& _solver);
	void SolverSetup(const Eigen::VectorXd& _initial_displacement, QuasiStaticFem& _solver);

	// for computing derivative
	void FindUnknownVariablesGlobalId();
	// pre-requst before calling the following 2 functions:
	// update FindUnknownVariablesGlobalId(); 
	void FindStiffnessSubmatrixOnUnknownNodes();
	void FindDerivativeOrderingMap();

	// jacobian matrix DF_free/Dx_free (it's a negative-definite matrix)
	// this is also known as the tangent stiffness matrix
	// row=col= 3 * free node number
	// pre-requst before calling this function:
	// update FindUnknownVariablesGlobalId();
	// update FindStiffnessSubmatrixOnUnknownNodes();
	// already performed quasi-static simulation (i.e. set solver_ ready)
	void Derivative_ForceOnFreeNodes_DeformedFreeNodes(Eigen::SparseMatrix<double>& _jacobian) const;
	// DF_free/Dx0, equal to a submatrix of K: rows of frees nodes on all columns
	// row= 3 * free node number
	// col= 3 * all node number
	// pre-requst before calling this function:
	// update FindUnknownVariablesGlobalId();
	void Derivative_ForceOnFreeNodes_OriginalAllNodes(Eigen::SparseMatrix<double>& _jacobian) const;



private:
	const TetrahedralMesh* tetrahedral_mesh_;
	TetrahedralElements elements_;
	double E_ = 0.0;
	double nu_ = 0.0;

	// for contact
	const tmd::TriangleMeshDistance* sdf_other_objects_ = nullptr;
	double collision_weight_ = 0.0;
	int sdf_function_index_ = 0;//for CubicLagrangeDiscreteGrid function query
	const Discregrid::CubicLagrangeDiscreteGrid* real_sdf_ = nullptr;
	// for quasi-static solver
	QuasiStaticFem solver_;
	int total_node_num_ = 0;
	int unknown_node_num_ = 0;// unknown node number = degree of freedom / 3 
	std::vector<std::pair<int, Eigen::Vector3d>> displacement_boundary_;// size = fixed node number
	std::vector<std::pair<int, Eigen::Vector3d>> force_boundary_;
	// notice elements inside have increasing order 
	std::vector<int> unknown_displace_node_id_;// usage e.g. _all_displacement[3 * unknown_displace_node_id_[iNode] + iCoordinate] = _unknown_displacement[3 * iNode + iCoordinate]

	// result
	Eigen::VectorXd all_displacement_;// size = total_node_num_ * 3
	Eigen::VectorXd nodal_force_;// (-1) * elastic force at all nodes; size = total_node_num_ * 3

	// for derivative info

	// size = 3*unknown_node number; 
	// elements in increasing order; 
	//_all_displacement[unknown_variables_global_id_[i]] = _unknown_displacement[i]
	std::vector<int> unknown_variables_global_id_;
	Eigen::SparseMatrix<double> stiffness_submatrix_unknown_nodes_;// K_Nf, this is also the hessian of elastic energy
	// size = total_node_num_ * 3
	// derivative_ordering_map_[i] gives a new index of the node variable i
	// it map fixed nodes to the top rows and free nodes to bottom rows in a vector
	// this ordering corresponds to the row ordering of Dx/Dx0 given by function: Derivative_Deformed_Original_Term1 and 2
	std::vector<int> derivative_ordering_map_;
	Eigen::SparseMatrix<double, Eigen::RowMajor> stiffness_matrix_row_major_;
};