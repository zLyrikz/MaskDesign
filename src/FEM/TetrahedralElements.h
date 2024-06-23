#pragma once
#include "../MeshProcessing/TetrahedralMesh.h"

#include <vector>
#include <Eigen/SparseCore>

//
// this is for compute and prodive information of each element of a mesh. 
// like:
// stifness matrix, volume force, surface force
// 
// we use linear polynomial element, uniform material
//
// strain matrix B= 6-by-12, L= 6-by-3, N= 3-by-12; 
// [B] = [L][N]; strain epsilon = B * displacement node delta
// here 6 is 6 strain functions, 3 is 3 displacement functions, 12 is degree of freedom of a linear tetrahedron element
//

//note!! check if volume sign consistantly positive
class TetrahedralElements
{
public:
	TetrahedralElements();
	~TetrahedralElements();

	// notice set all the parameters here in order before construct the linear system and solve
	void SetBaseMesh(const TetrahedralMesh* _mesh);
	void SetMaterial(double _E, double _nu);
	// If a node is not in the force and displacement boundary, it is considered as 0 force boundary condition
	void SetBoundaryCondition(const std::vector<std::pair<int, Eigen::Vector3d>>& _displacement_boundary,
		const std::vector<std::pair<int, Eigen::Vector3d>>& _force_boundary);
	// input force distribution, e.g. gravity would be density * g
	// default is 0 if function not called
	void SetUniformBodyForce(const Eigen::Vector3d& _body_force);

	void FindStiffnessMatrix();
	// construct equation include finding the stiffness matrix
	void ConstructEquation();// TODO make K and K_modify elements be 0 for those that should be 0
	void Solve();
	void ComputeBoundaryForceLoad();//find the node load on the displacement boundary
	void ComputeForceBoundaryForceLoad(std::vector<Eigen::Vector3d>& _load) const;//find the node load on the force boundary

	const Eigen::VectorXd* GetDisplacement() const;
	const std::vector<Eigen::Vector3d>* GetBoundaryForceLoad() const;// the force node load on the displacement boundary condition nodes
	const std::vector<std::pair<int, Eigen::Vector3d>>* GetDisplacementBoundary() const;
	const std::vector<std::pair<int, Eigen::Vector3d>>* GetForceBoundary() const;
	const Eigen::SparseMatrix<double>* GetStiffnessMatrix() const;

	// functions giving more information
	void computeVolumesGivenDisplacement(const Eigen::VectorXd& _displacement, std::vector<double>& _volumes) const;
	void computeElementEnergyGivenDisplacement(const Eigen::VectorXd& _displacement, std::vector<double>& _energy) const;// energy of each element
	bool IsVolumePositive() const;
	//void GetVolumes(std::vector<double>& _)

private:
	// assume each element has properly oriented vertices: i j k m
	// such that when compute volume with determinant, we get positive value
	void computeVolumes();// also find the "volume matrix"
	void AddUniformBodyForce(const Eigen::Vector3d& _body_force);//body force load need volume to compute
	void computeStrainMatrix();// this is represented as [B] = [L][N]; called after computeVolumes()
	void computeElementStiffnessMatrix();
	void assembleStiffnessMatrix();
	void computeForceLoad();// include the surface force boundary conditions 

	// setting stiffness matrix on the known boundary to 1;
	// move the known boundary with the stiffness to force load (right hand side)
	// then we can solve the linear systems KU=F
	void addDisplacementBoundaryCondition();

	// utility functions
	void ComputeBoundaryForceLoad(std::vector<Eigen::Vector3d>& _load, const std::vector<std::pair<int, Eigen::Vector3d>>& _boundary_id) const;

private:
	const TetrahedralMesh* tetrahedral_mesh_;

	// element info
	int totoal_dof_;// totoal degree of freedom = vertex number * 3
	int num_elements_;
	struct Element
	{
		//	  (1,x1,y1,z1
		//     1,x2,y2,z2
		//     1,x3,y3,z3
		//     1,x4,y4,z4)
		Eigen::Matrix4d volume_matrix_;
		double volume_ = 0.0;

		Eigen::MatrixXd B_;// strain matrix 6-by-12
		std::array<Eigen::Vector4d, 4> nodal_;// 4 coefficents of each linear nodal functions: (c1 c2 c3 c4) * (1 x y z)^T

		Eigen::MatrixXd Ke_;// stiffness matrix of this element 12-by-12
		Eigen::VectorXd body_force_ = Eigen::VectorXd::Zero(12);// equivalent nodal load of body force, size = element_dof_
	};
	std::vector<Element> elements_;
	Eigen::Vector3d body_force_density_;

	// linear system
	Eigen::SparseMatrix<double> K_;//stiffness matrix (without modified by adding boundary conditions; so it's singular)
	Eigen::SparseMatrix<double> K_modified_;// modified from K by setting 1 & 0 values according to displacement boundary conditions
	std::vector<std::pair<int, Eigen::Vector3d>> displacement_boundary_;// displacement boundary condition; node index and value pair
	std::vector<std::pair<int, Eigen::Vector3d>> force_boundary_;		// force boundary condition; the equivalent surface force load on the boundary node
	Eigen::VectorXd force_load_;// the F on the right hand side of the equilibrium
	
	// solutions
	Eigen::VectorXd displacement_;
	std::vector<Eigen::Vector3d> boundary_force_;// size = displacement boundary condition node number

	// material parameters
	double E_;// Young's modulus
	double nu_;// Poisson's ratio
	Eigen::MatrixXd D_;// constitutive matrix 6-by-6

	// constant values
	const int strain_dof_ = 6;		// strian degree of freedom = 6
	const int displacement_dof_ = 3;// displacement degree of freedom = 3
	const int element_dof_ = 12;	// degree of freedom of each element (3d * 4node) = 12
};

