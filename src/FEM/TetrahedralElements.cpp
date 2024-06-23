#include "TetrahedralElements.h"
#include <Eigen/LU>
#include <Eigen/QR>
#include<Eigen/IterativeLinearSolvers>
#include <iostream>
using std::cout;
using std::endl;
TetrahedralElements::TetrahedralElements()
{
	tetrahedral_mesh_ = nullptr;
	totoal_dof_ = 0;
	num_elements_ = 0;
	E_ = 0.0;
	nu_ = 0.0;
	body_force_density_ = Eigen::Vector3d::Zero();
}

TetrahedralElements::~TetrahedralElements()
{
}

void TetrahedralElements::SetBaseMesh(const TetrahedralMesh* _mesh)
{
	tetrahedral_mesh_ = _mesh;
	num_elements_ = _mesh->GetTopologyTetrahedra()->rows();
	totoal_dof_ = _mesh->GetVertex()->rows() * 3;
	elements_.resize(num_elements_);
}

void TetrahedralElements::SetMaterial(double _E, double _nu)
{
	if (E_ != _E || nu_ != _nu)
	{
		E_ = _E;
		nu_ = _nu;

		D_ = Eigen::MatrixXd::Zero(strain_dof_, strain_dof_);

		double common_coefficient = E_ / (1 + nu_) / (1 - 2 * nu_);
		D_(0, 0) = common_coefficient * (1 - nu_);
		D_(1, 1) = common_coefficient * (1 - nu_);
		D_(2, 2) = common_coefficient * (1 - nu_);
		D_(3, 3) = E_ / (1 + nu_) / 2.0;
		D_(4, 4) = E_ / (1 + nu_) / 2.0;
		D_(5, 5) = E_ / (1 + nu_) / 2.0;
		D_(0, 1) = common_coefficient * nu_;
		D_(1, 0) = common_coefficient * nu_;
		D_(0, 2) = common_coefficient * nu_;
		D_(2, 0) = common_coefficient * nu_;
		D_(1, 2) = common_coefficient * nu_;
		D_(2, 1) = common_coefficient * nu_;
	}
}

void TetrahedralElements::SetBoundaryCondition(const std::vector<std::pair<int, Eigen::Vector3d>>& _displacement_boundary, 
	const std::vector<std::pair<int, Eigen::Vector3d>>& _force_boundary)
{
	displacement_boundary_ = _displacement_boundary;
	force_boundary_ = _force_boundary;
}

void TetrahedralElements::AddUniformBodyForce(const Eigen::Vector3d& _body_force)
{
	for (auto& iElement : elements_)
	{
		for (int iDisplace = 0; iDisplace < displacement_dof_; ++iDisplace)
		{
			//note the integration of a linear node function on the element = volume_ / 4.0
			double force = _body_force(iDisplace) * iElement.volume_ / 4.0;
			for (int iNode = 0; iNode < 4; ++iNode)
			{		
				iElement.body_force_(3 * iNode + iDisplace) = force;
				//cout << "force " << iElement.body_force_(3 * iNode + iDisplace) << endl;
			}
		}
	}	
}

void TetrahedralElements::SetUniformBodyForce(const Eigen::Vector3d& _body_force)
{
	body_force_density_ = _body_force;
}

void TetrahedralElements::FindStiffnessMatrix()
{
	if (tetrahedral_mesh_ != nullptr)
	{
		computeVolumes();
		computeStrainMatrix();
		computeElementStiffnessMatrix();
		assembleStiffnessMatrix();
	}
	else
	{
		cout << "WARNING TetrahedralElements::FindStiffnessMatrix() MESH IS NULL" << endl;
	}
}

void TetrahedralElements::ConstructEquation()
{
	if (tetrahedral_mesh_ != nullptr)
	{
		FindStiffnessMatrix();

		AddUniformBodyForce(body_force_density_);
		computeForceLoad();
		addDisplacementBoundaryCondition();
	}
	else
	{
		cout << "WARNING TetrahedralElements::ConstructEquation() MESH IS NULL" << endl;
	}
}

void TetrahedralElements::Solve()
{
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> conjugate_gradient;
	conjugate_gradient.compute(K_modified_);

	if (conjugate_gradient.info() != Eigen::Success) {
		cout << "WARNING TetrahedralElements::Solve() decomposition failed" << endl;
		return;
	}

	assert(force_load_.rows() == totoal_dof_ && K_modified_.rows() == totoal_dof_ && K_modified_.cols() == totoal_dof_);

	// TODO By default the iterations start with x=0 as an initial guess of the solution. One can control the start using the solveWithGuess() method.
	displacement_ = conjugate_gradient.solve(force_load_);
	if (conjugate_gradient.info() != Eigen::Success) {
		cout << "WARNING TetrahedralElements::Solve() solving failed" << endl;
		return;
	}
}

void TetrahedralElements::ComputeBoundaryForceLoad()
{
	ComputeBoundaryForceLoad(boundary_force_, displacement_boundary_);
}

void TetrahedralElements::ComputeForceBoundaryForceLoad(std::vector<Eigen::Vector3d>& _load) const
{
	ComputeBoundaryForceLoad(_load, force_boundary_);
}

const Eigen::VectorXd* TetrahedralElements::GetDisplacement() const
{
	return &displacement_;
}

const std::vector<Eigen::Vector3d>* TetrahedralElements::GetBoundaryForceLoad() const
{
	return &boundary_force_;
}

const std::vector<std::pair<int, Eigen::Vector3d>>* TetrahedralElements::GetDisplacementBoundary() const
{
	return &displacement_boundary_;
}

const std::vector<std::pair<int, Eigen::Vector3d>>* TetrahedralElements::GetForceBoundary() const
{
	return &force_boundary_;
}

const Eigen::SparseMatrix<double>* TetrahedralElements::GetStiffnessMatrix() const
{
	return &K_;
}

void TetrahedralElements::computeVolumesGivenDisplacement(const Eigen::VectorXd& _displacement, std::vector<double>& _volumes) const
{
	if (_displacement.rows() == totoal_dof_)
	{
		TetrahedralMesh displaced_mesh(*tetrahedral_mesh_);
		displaced_mesh.TranslateIndividually(_displacement);
		_volumes.resize(num_elements_);
		for (int iTetrahedra = 0; iTetrahedra < num_elements_; ++iTetrahedra)
		{
			// get vertex
			std::array<Eigen::Vector3d, 4> vertex;
			displaced_mesh.GetTetrahedron(vertex, iTetrahedra);
			Eigen::Matrix4d volume_matrix;
			// get volume_matrix
			for (int iRow = 0; iRow < 4; ++iRow)
			{
				volume_matrix(iRow, 0) = 1.0;
			}
			for (int iRow = 0; iRow < 4; ++iRow)
			{
				for (int iCol = 1; iCol < 4; ++iCol)
				{
					volume_matrix(iRow, iCol) = vertex[iRow](iCol - 1);
				}
			}

			//get volume
			_volumes[iTetrahedra] = 1.0 / 6.0 * volume_matrix.determinant();
		}
	}

}

void TetrahedralElements::computeElementEnergyGivenDisplacement(const Eigen::VectorXd& _displacement, std::vector<double>& _energy) const
{
	_energy.resize(num_elements_);
	for (int iElement = 0; iElement < num_elements_; ++iElement)
	{
		// find the displace of current element
		Eigen::VectorXd displace_element(12);
		for (int iNode = 0; iNode < 4; ++iNode)
		{
			int global_node_id = (*tetrahedral_mesh_->GetTopologyTetrahedra())(iElement, iNode);

			for (int iCoordinate = 0; iCoordinate < 3; ++iCoordinate)
			{
				displace_element[3 * iNode + iCoordinate] = _displacement[3 * global_node_id + iCoordinate];
			}
		}
		_energy[iElement] = 0.5 * displace_element.transpose() * elements_[iElement].Ke_ * displace_element;
	}
}

bool TetrahedralElements::IsVolumePositive() const
{
	bool is_positive = true;
	for (auto& iElement : elements_)
	{
		if (iElement.volume_ < 0)
		{
			is_positive = false;
			break;
		}
	}

	return is_positive;
}

void TetrahedralElements::computeVolumes()
{
	// volume =
	// 1/6 * det(1,x1,y1,z1;
	//			 1,x2,y2,z2;
	//			 1,x3,y3,z3;
	//		     1,x4,y4,z4)
	// = 1/6 * det(x2-x1,y2-y1,z2-z1;
	//             x3-x1,y3-y1,z3-z1;
	//             x4-x1,y4-y1,z4-z1)

	assert(elements_.size() == num_elements_);
	//bool has_negative_volume = false;
	for (int iTetrahedra = 0; iTetrahedra < num_elements_; ++iTetrahedra)
	{
		// get vertex
		std::array<Eigen::Vector3d, 4> vertex;
		tetrahedral_mesh_->GetTetrahedron(vertex, iTetrahedra);

		// get volume_matrix
		for (int iRow = 0; iRow < 4; ++iRow)
		{
			elements_[iTetrahedra].volume_matrix_(iRow, 0) = 1.0;
		}
		for (int iRow = 0; iRow < 4; ++iRow)
		{
			for (int iCol = 1; iCol < 4; ++iCol)
			{
				elements_[iTetrahedra].volume_matrix_(iRow, iCol) = vertex[iRow](iCol - 1);
			}
		}

		//get volume
		elements_[iTetrahedra].volume_ = 1.0 / 6.0 * elements_[iTetrahedra].volume_matrix_.determinant();

		if (elements_[iTetrahedra].volume_ < 0)
		{
			//has_negative_volume = true;
			cout <<"iTetrahedra="<< iTetrahedra << ", volume negative=" << elements_[iTetrahedra].volume_ << endl;
			//cout << elements_[iTetrahedra].volume_matrix_ << endl;

		}

	}

}

void TetrahedralElements::computeStrainMatrix()
{
	for (auto& iElement : elements_)
	{
		// find nodal function coefficients
		Eigen::ColPivHouseholderQR<Eigen::Matrix4d> decomposed_A(iElement.volume_matrix_);
		for (int iNodal = 0; iNodal < 4; ++iNodal)
		{
			// solve linear system AX=node value
			// then X is the coefficient of iNodal node function
			// see Eigen linear system ref https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html

			Eigen::Vector4d node_value(0.0, 0.0, 0.0, 0.0);
			node_value(iNodal) = 1.0;

			// solve 
			iElement.nodal_[iNodal] = decomposed_A.solve(node_value);
		}
			
		// see B formula at FEM-Zhang Nian Mei
		iElement.B_ = Eigen::MatrixXd::Zero(strain_dof_, element_dof_);
		for (int iNode = 0; iNode < 4; ++iNode)
		{
			//construct B
			iElement.B_(0, 3 * iNode + 0) = iElement.nodal_[iNode](1);
			iElement.B_(1, 3 * iNode + 1) = iElement.nodal_[iNode](2);
			iElement.B_(2, 3 * iNode + 2) = iElement.nodal_[iNode](3);

			iElement.B_(3, 3 * iNode + 0) = iElement.nodal_[iNode](2);
			iElement.B_(3, 3 * iNode + 1) = iElement.nodal_[iNode](1);

			iElement.B_(4, 3 * iNode + 1) = iElement.nodal_[iNode](3);
			iElement.B_(4, 3 * iNode + 2) = iElement.nodal_[iNode](2);

			iElement.B_(5, 3 * iNode + 0) = iElement.nodal_[iNode](3);
			iElement.B_(5, 3 * iNode + 2) = iElement.nodal_[iNode](1);
		}
			
	}
}

void TetrahedralElements::computeElementStiffnessMatrix()
{
	for (auto& iElement : elements_)
	{
		iElement.Ke_ = iElement.volume_ * iElement.B_.transpose() * D_ * iElement.B_;
	}
}

void TetrahedralElements::assembleStiffnessMatrix()
{
	K_ = Eigen::SparseMatrix<double>();
	K_.resize(totoal_dof_, totoal_dof_);
	std::vector<Eigen::Triplet<double>> triplets;
	for (int iElement = 0; iElement < num_elements_; ++iElement)
	{
		// put each element stiffness matrix Ke to K
		// Ke is also consisted of 4 node blocks

		for (int iRow = 0; iRow < 4; ++iRow)
		{
			for (int iCol = 0; iCol < 4; ++iCol)
			{
				// now in the element, we have iRow node and iCol node
				// find the global index the two nodes
				int global_id_row = (*tetrahedral_mesh_->GetTopologyTetrahedra())(iElement, iRow);
				int global_id_col = (*tetrahedral_mesh_->GetTopologyTetrahedra())(iElement, iCol);

				// we have 3 displacement for each node
				// fill them one by one into K
				// this is the same as put the block in Ke to the block in K (i.e. view it as block operation)
				for (int iDisplacementRow = 0; iDisplacementRow < displacement_dof_; ++iDisplacementRow)
				{
					for (int iDisplacementCol = 0; iDisplacementCol < displacement_dof_; ++iDisplacementCol)
					{
						int row_K = displacement_dof_ * global_id_row + iDisplacementRow;
						int col_K = displacement_dof_ * global_id_col + iDisplacementCol;
						int row_Ke = displacement_dof_ * iRow + iDisplacementRow;
						int col_Ke = displacement_dof_ * iCol + iDisplacementCol;

						triplets.push_back(Eigen::Triplet<double>(row_K, col_K, elements_[iElement].Ke_(row_Ke, col_Ke)));
					}
				}
			}
		}
	}
	K_.setFromTriplets(triplets.begin(), triplets.end());
}

void TetrahedralElements::computeForceLoad()
{
	force_load_ = Eigen::VectorXd::Zero(totoal_dof_);

	// add body force
	for (int iElement = 0; iElement < num_elements_; ++iElement)
	{
		for (int iNode = 0; iNode < 4; ++iNode)
		{
			int node_id = (*tetrahedral_mesh_->GetTopologyTetrahedra())(iElement, iNode);
			// each node has 3 displacement degree of freedom
			for (int iDisplace = 0; iDisplace < displacement_dof_; ++iDisplace)
			{
				int row_id = displacement_dof_ * node_id + iDisplace;
				force_load_(row_id) += elements_[iElement].body_force_(displacement_dof_ * iNode + iDisplace);
			}
		}
	}

	//add surface force of the boundary force condition
	for (const std::pair<int, Eigen::Vector3d>& iBoundary : force_boundary_)
	{
		// node index is iBoundary.first
		for (int iDisplace = 0; iDisplace < displacement_dof_; ++iDisplace)
		{
			int row_id = displacement_dof_ * iBoundary.first + iDisplace;
			force_load_(row_id) += iBoundary.second(iDisplace);
		}
	}
}

void TetrahedralElements::addDisplacementBoundaryCondition()
{
	K_modified_ = K_;

	// change force load and construct K_modified_
	for (const std::pair<int, Eigen::Vector3d>& iBoundary : displacement_boundary_)
	{
		for (int iDisplace = 0; iDisplace < 3; ++iDisplace)
		{
			int column_id = displacement_dof_ * iBoundary.first + iDisplace;
			for (int iRow = 0; iRow < totoal_dof_; ++iRow)
			{
				double k = K_.coeff(iRow, column_id);
				if (abs(k) > 1e-10)
				{
					force_load_(iRow) -= k * iBoundary.second(iDisplace);

					// set all related columns and rows to 0
					K_modified_.coeffRef(iRow, column_id) = 0.0;
					K_modified_.coeffRef(column_id, iRow) = 0.0;
				}
			}
		}
	}
	// set modified K diagnal element to 1
	// set corresponding force load to the known displacement 
	for (const std::pair<int, Eigen::Vector3d>& iBoundary : displacement_boundary_)
	{
		for (int iDisplace = 0; iDisplace < 3; ++iDisplace)
		{
			int column_id = displacement_dof_ * iBoundary.first + iDisplace;

			force_load_(column_id) = iBoundary.second(iDisplace);
			K_modified_.coeffRef(column_id, column_id) = 1.0;
		}
	}
}

void TetrahedralElements::ComputeBoundaryForceLoad(std::vector<Eigen::Vector3d>& _load, const std::vector<std::pair<int, Eigen::Vector3d>>& _boundary_id) const
{
	_load.clear();
	_load.reserve(_boundary_id.size());
	// just use K*U to get the load
	for (const std::pair<int, Eigen::Vector3d>& iBoundary : _boundary_id)
	{
		Eigen::Vector3d node_load;
		for (int iDisplace = 0; iDisplace < displacement_dof_; ++iDisplace)
		{
			int index = displacement_dof_ * iBoundary.first + iDisplace;

			// load = K*U
			node_load(iDisplace) = K_.row(index) * displacement_;
		}
		_load.push_back(std::move(node_load));
	}
}
