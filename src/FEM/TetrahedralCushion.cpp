#include "TetrahedralCushion.h"
#include <iostream>
using std::cout;
using std::endl;
TetrahedralCushion::TetrahedralCushion()
{
}

TetrahedralCushion::~TetrahedralCushion()
{
}

void TetrahedralCushion::Tetrahedralize(const CushionSurface* _cushion, int _num_layer, double _height)
{
	if (TopologyTetrahedralize(_cushion, _num_layer))
	{
		height_ = _height;
		FindNodes();
		tet_cushion_.SetMesh_PossibleSameTopology(nodes_, topology_tetrahedra_, topology_boundary_repeated_);
		//cout << "number tetrahedra=" << topology_tetrahedra_.rows() << endl;

	}
}

void TetrahedralCushion::GetTetrahedralCushion(TetrahedralMesh* _tet_cushion) const
{
	*_tet_cushion = tet_cushion_;
}

const TetrahedralMesh* TetrahedralCushion::GetTetrahedralCushion() const
{
	return &tet_cushion_;
}

void TetrahedralCushion::GetFixConnectorBoundaryCondition(std::vector<std::pair<int, Eigen::Vector3d>>& _displacement_boundary)
{
	FindFixConnectorBoundaryCondition();
	_displacement_boundary = displacement_boundary_;
}

void TetrahedralCushion::GetCushionSurfaceNodes(std::vector<int>& _nodes)
{
	// update cushion surface nodes first
	int surface_node_num = num_v_ * num_u_;
	if (cushion_surface_nodes_.size() != surface_node_num)
	{
		cushion_surface_nodes_.resize(surface_node_num);
		int count_node = 0;

		for (int iU = 0; iU < num_u_; ++iU)
		{
			for (int iV = 0; iV < num_v_; ++iV)
			{
				cushion_surface_nodes_[count_node] = getGlobalIdFromNodeId(iU, iV, 0);
				++count_node;
			}
		}		
	}

	_nodes = cushion_surface_nodes_;
}

const std::vector<std::vector<int>>* TetrahedralCushion::GetFacesOfNodes() const
{
	return &faces_of_nodes_;
}

void TetrahedralCushion::GetNodesOfSurfaceFace(Eigen::Vector3d& _node1, Eigen::Vector3d& _node2, Eigen::Vector3d& _node3, int _face) const
{
	_node1 = nodes_.row(topology_boundary_(_face, 0));
	_node2 = nodes_.row(topology_boundary_(_face, 1));
	_node3 = nodes_.row(topology_boundary_(_face, 2));
}

void TetrahedralCushion::GetNodesIdOfSurfaceFace(int& _node1, int& _node2, int& _node3, int _face) const
{
	_node1 = topology_boundary_(_face, 0);
	_node2 = topology_boundary_(_face, 1);
	_node3 = topology_boundary_(_face, 2);
}

int TetrahedralCushion::GetUNumber() const
{
	return num_u_;
}

int TetrahedralCushion::GetVNumber() const
{
	return num_v_;
}

int TetrahedralCushion::GetTriangleNumber() const
{
	return topology_boundary_.rows();
}

int TetrahedralCushion::GetNodeNumber() const
{
	return nodes_.rows();
}

void TetrahedralCushion::GetNode(int _id, Eigen::Vector3d& _node) const
{
	if (_id < nodes_.rows())
	{
		_node = nodes_.row(_id);
	}
}

void TetrahedralCushion::GetDerivative_Nodes_CrossSectionKeyParameters(Eigen::MatrixXd& _Dx0_Dp)
{
	std::vector<std::vector<std::vector<std::array<Eigen::Matrix3Xd, 2>>>> half1_rest;
	std::vector<std::vector<std::vector<std::array<Eigen::Matrix3Xd, 2>>>> half2_rest;
	std::vector<std::vector<Eigen::Matrix3Xd>> half1_point1_length;
	std::vector<std::vector<Eigen::Matrix3Xd>> tangent_angle;
	std::vector<std::vector<Eigen::Matrix3Xd>> half2_point1_length;
	std::vector<std::vector<std::array<Eigen::Matrix3Xd, 2>>> point0;
	std::vector<std::vector<std::vector<std::array<Eigen::Matrix3Xd, 2>>>> normal_half1_rest;
	std::vector<std::vector<std::vector<std::array<Eigen::Matrix3Xd, 2>>>> normal_half2_rest;
	std::vector<std::vector<Eigen::Matrix3Xd>> normal_half1_point1_length;
	std::vector<std::vector<Eigen::Matrix3Xd>> normal_tangent_angle;
	std::vector<std::vector<Eigen::Matrix3Xd>> normal_half2_point1_length;
	std::vector<std::vector<std::array<Eigen::Matrix3Xd, 2>>> normal_point0;
	cushion_->getDerivative_DiscreteCushionAndCrossSectionNormal_KeyControlPoints(
		half1_rest, half2_rest, half1_point1_length, tangent_angle, half2_point1_length, point0,
		normal_half1_rest, normal_half2_rest, normal_half1_point1_length, normal_tangent_angle, normal_half2_point1_length, normal_point0);

	// (without point0) DoF of each key cross-section = 2 * (_half1_degree - 1) + 2 * (_half2_degree - 2) + 3
	// include point0, DoF += 2 -> DoF =  2 * (_half1_degree - 1) + 2 * (_half2_degree - 2) + 5
	int num_key_cross_section = cushion_->getKeyCrossSectionsNumber();
	int dim_each_cross_section = 2 * half1_rest.size() + 2 * half2_rest.size() + 5;
	int dim_parameter = dim_each_cross_section * num_key_cross_section;
	_Dx0_Dp.resize(3 * num_nodes_, dim_parameter);

	// here, we need to arrange parameters into one vector, assume it knows the order for now
	// order should be same with CrossSectionOptimizeInfo::GetXFromKeyCrossSections
	if (base_cushion_ == nullptr)// not vr version
	{
		int vertex_count = 0;
		for (int iLayer = 0; iLayer <= num_layer_; ++iLayer)
		{
			double depth = double(iLayer) / double(num_layer_) * height_;
			for (int iU = 0; iU < num_u_; ++iU)
			{
				for (int iV = 0; iV < num_v_; ++iV)
				{
					auto& x0_current_node = _Dx0_Dp.block(3 * vertex_count, 0, 3, dim_parameter);
					for (int iSection = 0; iSection < num_key_cross_section; ++iSection)
					{
						int base_index = iSection * dim_each_cross_section;

						x0_current_node.col(base_index) = half1_point1_length[iU][iV].col(iSection) + depth * normal_half1_point1_length[iU][iV].col(iSection);
						x0_current_node.col(base_index + 1) = tangent_angle[iU][iV].col(iSection) + depth * normal_tangent_angle[iU][iV].col(iSection);
						x0_current_node.col(base_index + 2) = half2_point1_length[iU][iV].col(iSection) + depth * normal_half2_point1_length[iU][iV].col(iSection);
						int start_index_half1 = base_index + 3;
						for (int iRest1 = 0; iRest1 < half1_rest.size(); ++iRest1)
						{
							x0_current_node.col(start_index_half1 + 2 * iRest1) = half1_rest[iRest1][iU][iV][0].col(iSection) + depth * normal_half1_rest[iRest1][iU][iV][0].col(iSection);
							x0_current_node.col(start_index_half1 + 2 * iRest1 + 1) = half1_rest[iRest1][iU][iV][1].col(iSection) + depth * normal_half1_rest[iRest1][iU][iV][1].col(iSection);
						}
						int start_index_half2 = start_index_half1 + 2 * half1_rest.size();
						for (int iRest2 = 0; iRest2 < half2_rest.size(); ++iRest2)
						{
							x0_current_node.col(start_index_half2 + 2 * iRest2) = half2_rest[iRest2][iU][iV][0].col(iSection) + depth * normal_half2_rest[iRest2][iU][iV][0].col(iSection);
							x0_current_node.col(start_index_half2 + 2 * iRest2 + 1) = half2_rest[iRest2][iU][iV][1].col(iSection) + depth * normal_half2_rest[iRest2][iU][iV][1].col(iSection);
						}

						// add point0 at the end
						int start_index_point0 = start_index_half2 + 2 * half2_rest.size();
						x0_current_node.col(start_index_point0) = point0[iU][iV][0].col(iSection) + depth * normal_point0[iU][iV][0].col(iSection);
						x0_current_node.col(start_index_point0 + 1) = point0[iU][iV][1].col(iSection) + depth * normal_point0[iU][iV][1].col(iSection);

					}

					++vertex_count;
				}
			}
		}
	}
	else
	{
		// vr version: x = (1-depth) * x_cushion_surface + depth * base_surface
		int vertex_count = 0;
		for (int iLayer = 0; iLayer <= num_layer_; ++iLayer)
		{
			double depth = double(iLayer) / double(num_layer_) * height_;
			for (int iU = 0; iU < num_u_; ++iU)
			{
				for (int iV = 0; iV < num_v_; ++iV)
				{
					auto& x0_current_node = _Dx0_Dp.block(3 * vertex_count, 0, 3, dim_parameter);
					for (int iSection = 0; iSection < num_key_cross_section; ++iSection)
					{
						int base_index = iSection * dim_each_cross_section;

						x0_current_node.col(base_index) = (1.0 - depth) * half1_point1_length[iU][iV].col(iSection);
						x0_current_node.col(base_index + 1) = (1.0 - depth) * tangent_angle[iU][iV].col(iSection);
						x0_current_node.col(base_index + 2) = (1.0 - depth) * half2_point1_length[iU][iV].col(iSection);
						int start_index_half1 = base_index + 3;
						for (int iRest1 = 0; iRest1 < half1_rest.size(); ++iRest1)
						{
							x0_current_node.col(start_index_half1 + 2 * iRest1) = (1.0 - depth) * half1_rest[iRest1][iU][iV][0].col(iSection);
							x0_current_node.col(start_index_half1 + 2 * iRest1 + 1) = (1.0 - depth) * half1_rest[iRest1][iU][iV][1].col(iSection);
						}
						int start_index_half2 = start_index_half1 + 2 * half1_rest.size();
						for (int iRest2 = 0; iRest2 < half2_rest.size(); ++iRest2)
						{
							x0_current_node.col(start_index_half2 + 2 * iRest2) = (1.0 - depth) * half2_rest[iRest2][iU][iV][0].col(iSection);
							x0_current_node.col(start_index_half2 + 2 * iRest2 + 1) = (1.0 - depth) * half2_rest[iRest2][iU][iV][1].col(iSection);
						}

						// add point0 at the end
						int start_index_point0 = start_index_half2 + 2 * half2_rest.size();
						x0_current_node.col(start_index_point0) = (1.0 - depth) * point0[iU][iV][0].col(iSection);
						x0_current_node.col(start_index_point0 + 1) = (1.0 - depth) * point0[iU][iV][1].col(iSection);

					}

					++vertex_count;
				}
			}
		}
	}
}


void TetrahedralCushion::ChangeNodes(Eigen::VectorXd& _nodes)
{
	if (_nodes.rows() == num_nodes_ * 3)
	{
		for (int iNode = 0; iNode < num_nodes_; ++iNode)
		{
			nodes_(iNode, 0) = _nodes(3 * iNode + 0);
			nodes_(iNode, 1) = _nodes(3 * iNode + 1);
			nodes_(iNode, 2) = _nodes(3 * iNode + 2);
		}
	}
}

void TetrahedralCushion::FindNodes()
{
	const std::vector<DiscreteCushionPiece>* cushion_pieces = cushion_->getDiscreteCushion();
	std::vector<std::vector<Eigen::Vector3d>> node_normals;
	cushion_->FindSurfaceNormal(node_normals);//TODO no need to re-compute cross-section derivative, just use the stored values

	nodes_.resize(num_nodes_, Eigen::NoChange);
	int vertex_count = 0;
	for (int iLayer = 0; iLayer <= num_layer_; ++iLayer)
	{
		double depth = double(iLayer) / double(num_layer_) * height_;
		for (int iU = 0; iU < num_u_; ++iU)
		{
			for (int iV = 0; iV < num_v_; ++iV)
			{
				nodes_.row(vertex_count) = cushion_pieces->at(iU).discrete_cross_section_[iV] + depth * node_normals[iU][iV];
				++vertex_count;
			}
		}
	}
}

void TetrahedralCushion::FindTetrahedaTopology()
{
	int num_hexahedra = num_layer_ * num_u_ * (num_v_ - 1);

	topology_tetrahedra_.resize(5 * num_hexahedra, Eigen::NoChange);
	int hexahedra_num_v = num_v_ - 1;

	int hexahedra_count = 0;
	for (int iLayer = 0; iLayer < num_layer_; ++iLayer)
	{
		for (int iU = 0; iU < num_u_; ++iU)
		{
			for (int iV = 0; iV < hexahedra_num_v; ++iV)
			{
				Eigen::MatrixX4i topology_tetrahedra;// 5 tetrahedra in one hexahedron
				FindTetrahedaTopology_InHexahedron(iU, iV, iLayer, topology_tetrahedra);
				topology_tetrahedra_.block(5 * hexahedra_count, 0, 5, 4) = topology_tetrahedra;
				hexahedra_count++;
			}
		}
	}
	
}

void TetrahedralCushion::FindSurfaceTopology()
{
	int hexahedra_num_v = num_v_ - 1;
	int num_surface_quads = 2 * (num_u_ * (num_v_ - 1)) + 2 * (num_u_ * num_layer_);
	topology_boundary_.resize(2 * num_surface_quads, Eigen::NoChange);

	// surface at first layer (cushion surface) layer = 0
	int quad_count = 0;
	for (int iU = 0; iU < num_u_; ++iU)
	{
		for (int iV = 0; iV < hexahedra_num_v; ++iV)
		{
			// quad 1234 in hexahedron
			int node1 = 0;
			int node2 = 0;
			int node3 = 0;
			int node4 = 0;
			HexahedronNode1234GlobalId(iU, iV, 0, node1, node2, node3, node4);
			if (DecidePattern(iU, iV, 0))
			{
				// pattern 1
				// triangle node 421  243
				topology_boundary_(2 * quad_count, 0) = node4;
				topology_boundary_(2 * quad_count, 1) = node2;
				topology_boundary_(2 * quad_count, 2) = node1;

				topology_boundary_(2 * quad_count + 1, 0) = node2;
				topology_boundary_(2 * quad_count + 1, 1) = node4;
				topology_boundary_(2 * quad_count + 1, 2) = node3;
			}
			else
			{
				// pattern 2
				// triangle node 213  314
				topology_boundary_(2 * quad_count, 0) = node2;
				topology_boundary_(2 * quad_count, 1) = node1;
				topology_boundary_(2 * quad_count, 2) = node3;

				topology_boundary_(2 * quad_count + 1, 0) = node3;
				topology_boundary_(2 * quad_count + 1, 1) = node1;
				topology_boundary_(2 * quad_count + 1, 2) = node4;
			}
			++quad_count;
		}
	}

	// surface at bottom layer
	for (int iU = 0; iU < num_u_; ++iU)
	{
		for (int iV = 0; iV < hexahedra_num_v; ++iV)
		{
			// quad 5678 in hexahedron
			int node5 = 0;
			int node6 = 0;
			int node7 = 0;
			int node8 = 0;
			HexahedronNode5678GlobalId(iU, iV, num_layer_ - 1, node5, node6, node7, node8);
			if (DecidePattern(iU, iV, num_layer_ - 1))
			{
				// pattern 1
				// triangle node 567  785
				topology_boundary_(2 * quad_count, 0) = node5;
				topology_boundary_(2 * quad_count, 1) = node6;
				topology_boundary_(2 * quad_count, 2) = node7;

				topology_boundary_(2 * quad_count + 1, 0) = node7;
				topology_boundary_(2 * quad_count + 1, 1) = node8;
				topology_boundary_(2 * quad_count + 1, 2) = node5;
			}
			else
			{
				// pattern 2
				// triangle node 568  678
				topology_boundary_(2 * quad_count, 0) = node5;
				topology_boundary_(2 * quad_count, 1) = node6;
				topology_boundary_(2 * quad_count, 2) = node8;

				topology_boundary_(2 * quad_count + 1, 0) = node6;
				topology_boundary_(2 * quad_count + 1, 1) = node7;
				topology_boundary_(2 * quad_count + 1, 2) = node8;
			}
			++quad_count;
		}
	}

	// surface at inner side (v=0)
	for (int iU = 0; iU < num_u_; ++iU)
	{
		for (int iLayer = 0; iLayer < num_layer_; ++iLayer)
		{
			// quad 1265 in hexahedron
			int node1 = 0;
			int node2 = 0;
			int node6 = 0;
			int node5 = 0;
			HexahedronNode1265GlobalId(iU, 0, iLayer, node1, node2, node6, node5);
			if (DecidePattern(iU, 0, iLayer))
			{
				// pattern 1
				// triangle node 125  652
				topology_boundary_(2 * quad_count, 0) = node1;
				topology_boundary_(2 * quad_count, 1) = node2;
				topology_boundary_(2 * quad_count, 2) = node5;

				topology_boundary_(2 * quad_count + 1, 0) = node6;
				topology_boundary_(2 * quad_count + 1, 1) = node5;
				topology_boundary_(2 * quad_count + 1, 2) = node2;
			}
			else
			{
				// pattern 2
				// triangle node 126  165
				topology_boundary_(2 * quad_count, 0) = node1;
				topology_boundary_(2 * quad_count, 1) = node2;
				topology_boundary_(2 * quad_count, 2) = node6;

				topology_boundary_(2 * quad_count + 1, 0) = node1;
				topology_boundary_(2 * quad_count + 1, 1) = node6;
				topology_boundary_(2 * quad_count + 1, 2) = node5;
			}
			++quad_count;
		}
	}

	// surface at outer(connector) side (v=num_v-2)
	for (int iU = 0; iU < num_u_; ++iU)
	{
		for (int iLayer = 0; iLayer < num_layer_; ++iLayer)
		{
			// quad 3487 in hexahedron
			int node3 = 0;
			int node4 = 0;
			int node8 = 0;
			int node7 = 0;
			HexahedronNode3487GlobalId(iU, num_v_ - 2, iLayer, node3, node4, node8, node7);
			if (DecidePattern(iU, num_v_ - 2, iLayer))
			{
				// pattern 1
				// triangle node 347  874
				topology_boundary_(2 * quad_count, 0) = node3;
				topology_boundary_(2 * quad_count, 1) = node4;
				topology_boundary_(2 * quad_count, 2) = node7;

				topology_boundary_(2 * quad_count + 1, 0) = node8;
				topology_boundary_(2 * quad_count + 1, 1) = node7;
				topology_boundary_(2 * quad_count + 1, 2) = node4;
			}
			else
			{
				// pattern 2
				// triangle node 483  738
				topology_boundary_(2 * quad_count, 0) = node4;
				topology_boundary_(2 * quad_count, 1) = node8;
				topology_boundary_(2 * quad_count, 2) = node3;

				topology_boundary_(2 * quad_count + 1, 0) = node7;
				topology_boundary_(2 * quad_count + 1, 1) = node3;
				topology_boundary_(2 * quad_count + 1, 2) = node8;
			}
			++quad_count;
		}
	}

	assert(quad_count == num_surface_quads);
}

void TetrahedralCushion::FindFixConnectorBoundaryCondition()
{
	int connector_node_num = (num_layer_ + 1) * num_u_;
	if (displacement_boundary_.size() != connector_node_num)
	{
		displacement_boundary_.resize(connector_node_num);
		int count_node = 0;
		for (int iLayer = 0; iLayer <= num_layer_; ++iLayer)
		{
			for (int iU = 0; iU < num_u_; ++iU)
			{
				int node_global_id = getGlobalIdFromNodeId(iU, num_v_ - 1, iLayer);
				displacement_boundary_[count_node] = std::make_pair(node_global_id , Eigen::Vector3d::Zero());
				++count_node;
			}
		}
	}
}

void TetrahedralCushion::FindFixBaseCushionBoundaryCondition()
{
	int fixed_node_num = num_v_ * num_u_;
	if (displacement_boundary_.size() != fixed_node_num)
	{
		displacement_boundary_.resize(fixed_node_num);
		int count_node = 0;
		for (int iU = 0; iU < num_u_; ++iU)
		{
			for (int iV = 0; iV < num_v_; ++iV)
			{
				int node_global_id = getGlobalIdFromNodeId(iU, iV, num_layer_);
				displacement_boundary_[count_node] = std::make_pair(node_global_id, Eigen::Vector3d::Zero());
				++count_node;
			}
		}	
	}
}

void TetrahedralCushion::FindFacesOfNodes()
{
	faces_of_nodes_.clear();
	faces_of_nodes_.resize(num_nodes_);
	for (int iFace = 0; iFace < topology_boundary_.rows(); ++iFace)
	{
		for (int iNode = 0; iNode < 3; ++iNode)
		{
			faces_of_nodes_[topology_boundary_(iFace, iNode)].push_back(iFace);
		}
	}
}

void TetrahedralCushion::FindSurfaceNodesOfNodes()
{
	surface_nodes_of_nodes_.resize(num_nodes_);
	for (int iNode = 0; iNode < faces_of_nodes_.size(); ++iNode)
	{
		std::vector<int>& adjacent_nodes = surface_nodes_of_nodes_[iNode];
		// loop over adjacent faces of this node
		for (int iFace = 0; iFace < faces_of_nodes_[iNode].size(); ++iFace)
		{
			// loop over nodes of this face
			for (int iFaceNode = 0; iFaceNode < 3; ++iFaceNode)
			{
				int current_adjacent_node = topology_boundary_(faces_of_nodes_[iNode][iFace], iFaceNode);

				// if node not stored in the list, store it
				if (std::find(adjacent_nodes.begin(), adjacent_nodes.end(), current_adjacent_node) != adjacent_nodes.end())
				{

					adjacent_nodes.push_back(current_adjacent_node);
				}

			}
		}
	}
}

void TetrahedralCushion::FindTetrahedaTopology_InHexahedron(int _u, int _v, int _layer, Eigen::MatrixX4i& _topology_tetrahedra) const
{
	int node1 = 0;
	int node2 = 0;
	int node3 = 0;
	int node4 = 0;
	int node5 = 0;
	int node6 = 0;
	int node7 = 0;
	int node8 = 0;
	HexahedronNodeGlobalId(_u, _v, _layer, node1, node2, node3, node4, node5, node6, node7, node8);
	if (DecidePattern(_u, _v, _layer))
	{
		FindTetrahedaTopology_InHexahedron_Pattern1(node1, node2, node3, node4, node5, node6, node7, node8, _topology_tetrahedra);
	}
	else
	{
		FindTetrahedaTopology_InHexahedron_Pattern2(node1, node2, node3, node4, node5, node6, node7, node8, _topology_tetrahedra);

	}
}

void TetrahedralCushion::HexahedronNodeGlobalId(int _u, int _v, int _layer,
	int& _node1, int& _node2, int& _node3, int& _node4, int& _node5, int& _node6, int& _node7, int& _node8) const
{
	// a hexahedron of eight local nodes
	// apply right hand side rule to the orientation of node 1234, the thumb should point to the layer of 5678
	// node 1          2		       3               4         
	//      (u,v,l)    (u+1,v,l)	  (u+1,v+1,l)     (u,v+1,l)
	//      5          6		       7               8        
	//      (u,v,l+1)  (u+1,v,l+1)	  (u+1,v+1,l+1)   (u,v+1,l+1)
	// 
	int next_u = (_u + 1) % num_u_;
	_node1 = getGlobalIdFromNodeId(_u, _v, _layer);
	_node2 = getGlobalIdFromNodeId(next_u, _v, _layer);
	_node3 = getGlobalIdFromNodeId(next_u, _v + 1, _layer);
	_node4 = getGlobalIdFromNodeId(_u, _v + 1, _layer);
	_node5 = getGlobalIdFromNodeId(_u, _v, _layer + 1);
	_node6 = getGlobalIdFromNodeId(next_u, _v, _layer + 1);
	_node7 = getGlobalIdFromNodeId(next_u, _v + 1, _layer + 1);
	_node8 = getGlobalIdFromNodeId(_u, _v + 1, _layer + 1);
}

void TetrahedralCushion::HexahedronNode1234GlobalId(int _u, int _v, int _layer, int& _node1, int& _node2, int& _node3, int& _node4) const
{
	int next_u = (_u + 1) % num_u_;
	_node1 = getGlobalIdFromNodeId(_u, _v, _layer);
	_node2 = getGlobalIdFromNodeId(next_u, _v, _layer);
	_node3 = getGlobalIdFromNodeId(next_u, _v + 1, _layer);
	_node4 = getGlobalIdFromNodeId(_u, _v + 1, _layer);
}

void TetrahedralCushion::HexahedronNode5678GlobalId(int _u, int _v, int _layer, int& _node5, int& _node6, int& _node7, int& _node8) const
{
	int next_u = (_u + 1) % num_u_;
	_node5 = getGlobalIdFromNodeId(_u, _v, _layer + 1);
	_node6 = getGlobalIdFromNodeId(next_u, _v, _layer + 1);
	_node7 = getGlobalIdFromNodeId(next_u, _v + 1, _layer + 1);
	_node8 = getGlobalIdFromNodeId(_u, _v + 1, _layer + 1);
}

void TetrahedralCushion::HexahedronNode1265GlobalId(int _u, int _v, int _layer, int& _node1, int& _node2, int& _node6, int& _node5) const
{
	int next_u = (_u + 1) % num_u_;
	_node1 = getGlobalIdFromNodeId(_u, _v, _layer);
	_node2 = getGlobalIdFromNodeId(next_u, _v, _layer);
	_node5 = getGlobalIdFromNodeId(_u, _v, _layer + 1);
	_node6 = getGlobalIdFromNodeId(next_u, _v, _layer + 1);
}

void TetrahedralCushion::HexahedronNode3487GlobalId(int _u, int _v, int _layer, int& _node3, int& _node4, int& _node8, int& _node7) const
{
	int next_u = (_u + 1) % num_u_;
	_node3 = getGlobalIdFromNodeId(next_u, _v + 1, _layer);
	_node4 = getGlobalIdFromNodeId(_u, _v + 1, _layer);
	_node7 = getGlobalIdFromNodeId(next_u, _v + 1, _layer + 1);
	_node8 = getGlobalIdFromNodeId(_u, _v + 1, _layer + 1);
}

void TetrahedralCushion::FindTetrahedaTopology_InHexahedron_Pattern1(
	int _node1, int _node2, int _node3, int _node4, int _node5, int _node6, int _node7, int _node8, 
	Eigen::MatrixX4i& _topology_tetrahedra) const
{
	// Pattern1:
	// tet 0       1      2       3      4
	//     1245   3427    6572    8754	 5724

	_topology_tetrahedra.resize(5, Eigen::NoChange);
	_topology_tetrahedra(0, 0) = _node1;
	_topology_tetrahedra(0, 1) = _node2;
	_topology_tetrahedra(0, 2) = _node4;
	_topology_tetrahedra(0, 3) = _node5;

	_topology_tetrahedra(1, 0) = _node3;
	_topology_tetrahedra(1, 1) = _node4;
	_topology_tetrahedra(1, 2) = _node2;
	_topology_tetrahedra(1, 3) = _node7;

	_topology_tetrahedra(2, 0) = _node6;
	_topology_tetrahedra(2, 1) = _node5;
	_topology_tetrahedra(2, 2) = _node7;
	_topology_tetrahedra(2, 3) = _node2;

	_topology_tetrahedra(3, 0) = _node8;
	_topology_tetrahedra(3, 1) = _node7;
	_topology_tetrahedra(3, 2) = _node5;
	_topology_tetrahedra(3, 3) = _node4;

	_topology_tetrahedra(4, 0) = _node5;
	_topology_tetrahedra(4, 1) = _node7;
	_topology_tetrahedra(4, 2) = _node2;
	_topology_tetrahedra(4, 3) = _node4;
}

void TetrahedralCushion::FindTetrahedaTopology_InHexahedron_Pattern2(
	int _node1, int _node2, int _node3, int _node4, int _node5, int _node6, int _node7, int _node8, 
	Eigen::MatrixX4i& _topology_tetrahedra) const
{

	// Pattern2:
	// tet 0       1      2       3      4
	//     2316   4138    5861    7683	 1638
	_topology_tetrahedra.resize(5, Eigen::NoChange);
	_topology_tetrahedra(0, 0) = _node2;
	_topology_tetrahedra(0, 1) = _node3;
	_topology_tetrahedra(0, 2) = _node1;
	_topology_tetrahedra(0, 3) = _node6;

	_topology_tetrahedra(1, 0) = _node4;
	_topology_tetrahedra(1, 1) = _node1;
	_topology_tetrahedra(1, 2) = _node3;
	_topology_tetrahedra(1, 3) = _node8;

	_topology_tetrahedra(2, 0) = _node5;
	_topology_tetrahedra(2, 1) = _node8;
	_topology_tetrahedra(2, 2) = _node6;
	_topology_tetrahedra(2, 3) = _node1;

	_topology_tetrahedra(3, 0) = _node7;
	_topology_tetrahedra(3, 1) = _node6;
	_topology_tetrahedra(3, 2) = _node8;
	_topology_tetrahedra(3, 3) = _node3;

	_topology_tetrahedra(4, 0) = _node1;
	_topology_tetrahedra(4, 1) = _node6;
	_topology_tetrahedra(4, 2) = _node3;
	_topology_tetrahedra(4, 3) = _node8;
}

int TetrahedralCushion::getGlobalIdFromNodeId(int _u, int _v, int _layer) const
{
	return _layer * (num_u_ * num_v_) + _u * num_v_ + _v;
}

bool TetrahedralCushion::DecidePattern(int _u, int _v, int _layer) const
{
	return ((_u + _v + _layer) % 2) == 0;
}

bool TetrahedralCushion::TopologyTetrahedralize(const CushionSurface* _cushion, int _num_layer)
{
	cushion_ = _cushion;
	const std::vector<DiscreteCushionPiece>* cushion_pieces = cushion_->getDiscreteCushion();
	int num_u = cushion_pieces->size();
	if (num_u % 2 != 0)
	{
		cout << "!!!ERROR!!! FROM TetrahedralCushion input u number not even" << endl;
		return false;
	}
	else
	{
		int num_v = cushion_pieces->at(0).discrete_cross_section_.size();
		if (num_u != num_u_ || num_v != num_v_ || num_layer_ != _num_layer)
		{
			// changing the topology
			num_u_ = num_u;
			num_v_ = num_v;
			num_layer_ = _num_layer;
			num_nodes_ = (num_layer_ + 1) * num_u_ * num_v_;

			FindTetrahedaTopology();
			FindSurfaceTopology();

			// get SurfaceTopology doubled and construct the tetrahedral mesh
			topology_boundary_repeated_.resize(2 * topology_boundary_.rows(), Eigen::NoChange);
			for (int iTrangle = 0; iTrangle < topology_boundary_.rows(); ++iTrangle)
			{
				// first triangle
				topology_boundary_repeated_.row(2 * iTrangle) = topology_boundary_.row(iTrangle);
				// reversed oriented triangle
				topology_boundary_repeated_(2 * iTrangle + 1, 0) = topology_boundary_(iTrangle, 2);
				topology_boundary_repeated_(2 * iTrangle + 1, 1) = topology_boundary_(iTrangle, 1);
				topology_boundary_repeated_(2 * iTrangle + 1, 2) = topology_boundary_(iTrangle, 0);
			}

			// more topology information

			FindFacesOfNodes();


			//FindSurfaceNodesOfNodes();
		}
		//height_ = _height;
		//FindNodes();
		//tet_cushion_.SetMesh_PossibleSameTopology(nodes_, topology_tetrahedra_, topology_boundary_repeated_);

		//cout << "---------------number of tetrahedra = " << topology_tetrahedra_.rows() <<"-----------------------------" << endl;;
	}
	return true;
}
