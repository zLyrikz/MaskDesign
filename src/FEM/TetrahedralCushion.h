#pragma once
#include "../MeshProcessing/TetrahedralMesh.h"
#include "../MaskInterface/CushionSurface.h"

//
// 
//   8 ______7
//   /      /|
// 5/_____6/ |
//  |      | |3
//  |______|/
//  1      2
// 
// a hexahedron of eight local nodes
// apply right hand side rule to the orientation of node 1234, the thumb should point to the layer of 5678
// hexahedron (u,v,l)
// node 1          2		       3               4         
//      (u,v,l)    (u+1,v,l)	  (u+1,v+1,l)     (u,v+1,l)
//      5          6		       7               8        
//      (u,v,l+1)  (u+1,v,l+1)	  (u+1,v+1,l+1)   (u,v+1,l+1)
// 
// 
// they give 5 tetrahedra:
// Pattern1:
// tet  0       1      2       3      4
// node 1245   3427    6572    8754	 5724
//
// Pattern2:
// tet  0       1      2       3      4
// node 2316   4138    5861    7683	 1638
//
// Same patterns CANNOT be adjacent to each other. Only different patterns can.
// so (-1)^(u+v+l) can be an indicator of which pattern to use (or (u+v+l)%2 )
// NOTICE that num_u have to be even! since cushion is periodic on u

class TetrahedralCushion
{
public:
	TetrahedralCushion();
	~TetrahedralCushion();

	// num_layer is how many layers of hexahedra along the height
	void Tetrahedralize(const CushionSurface* _cushion, int _num_layer, double _height = 3.0);
	void GetTetrahedralCushion(TetrahedralMesh* _tet_cushion) const;// copy TetrahderalMesh
	const TetrahedralMesh* GetTetrahedralCushion() const;// only pass the pointer
	void GetFixConnectorBoundaryCondition(std::vector<std::pair<int, Eigen::Vector3d>>& _displacement_boundary);
	void GetCushionSurfaceNodes(std::vector<int>& _nodes);
	const std::vector<std::vector<int>>* GetFacesOfNodes() const;
	void GetNodesOfSurfaceFace(Eigen::Vector3d& _node1, Eigen::Vector3d& _node2, Eigen::Vector3d& _node3, int _face) const;
	void GetNodesIdOfSurfaceFace(int& _node1, int& _node2, int& _node3, int _face) const;
	int GetUNumber() const;
	int GetVNumber() const;
	int GetTriangleNumber() const;
	int GetNodeNumber() const;
	void GetNode(int _id, Eigen::Vector3d& _node) const;
	//int GetTetrahedraNumber() const;
	// derivative

	// jacobian matrix Dx0/Dp
	// row=3*num nodes; col=dim of p
	void GetDerivative_Nodes_CrossSectionKeyParameters(Eigen::MatrixXd& _Dx0_Dp);

public:
	// additional function
	// only change nodes to input ones 
	void ChangeNodes(Eigen::VectorXd& _nodes);

private:
	void FindNodes();
	void FindTetrahedaTopology();
	void FindSurfaceTopology();// surface mesh, boundary of the tetrahedral mesh

	void FindFixConnectorBoundaryCondition();
	void FindFixBaseCushionBoundaryCondition();
	void FindFacesOfNodes();
	void FindSurfaceNodesOfNodes();

	void FindTetrahedaTopology_InHexahedron(int _u, int _v, int _layer, Eigen::MatrixX4i& _topology_tetrahedra) const;
	void HexahedronNodeGlobalId(int _u, int _v, int _layer,
		int& _node1, int& _node2, int& _node3, int& _node4, int& _node5, int& _node6, int& _node7, int& _node8) const;
	void HexahedronNode1234GlobalId(int _u, int _v, int _layer,
		int& _node1, int& _node2, int& _node3, int& _node4) const;
	void HexahedronNode5678GlobalId(int _u, int _v, int _layer,
		int& _node5, int& _node6, int& _node7, int& _node8) const;
	void HexahedronNode1265GlobalId(int _u, int _v, int _layer,
		int& _node1, int& _node2, int& _node6, int& _node5) const;
	void HexahedronNode3487GlobalId(int _u, int _v, int _layer,
		int& _node3, int& _node4, int& _node8, int& _node7) const;
	void FindTetrahedaTopology_InHexahedron_Pattern1(
		int _node1, int _node2, int _node3, int _node4, int _node5, int _node6, int _node7, int _node8, 
		Eigen::MatrixX4i& _topology_tetrahedra) const;
	void FindTetrahedaTopology_InHexahedron_Pattern2(
		int _node1, int _node2, int _node3, int _node4, int _node5, int _node6, int _node7, int _node8, 
		Eigen::MatrixX4i& _topology_tetrahedra) const;
	int getGlobalIdFromNodeId(int _u, int _v, int _layer) const;
	bool DecidePattern(int _u, int _v, int _layer) const;// splitting pattern of the hexahedron(u,v,l); true = pattern 1; false = pattern 2

	// ultility function
	bool TopologyTetrahedralize(const CushionSurface* _cushion, int _num_layer);


private:
	const CushionSurface* cushion_ = nullptr;
	const CushionSurface* base_cushion_ = nullptr;
	TetrahedralMesh tet_cushion_;

	int num_u_ = 0;
	int num_v_ = 0;
	int num_layer_ = 0;
	double height_ = 0;

	// node id = layer * (num_u * num_v) + u * num_v + v
	int num_nodes_ = 0.0;
	Eigen::MatrixX3d nodes_;
	Eigen::MatrixX4i topology_tetrahedra_;
	Eigen::MatrixX3i topology_boundary_;// not repeat a triangle twice
	Eigen::MatrixX3i topology_boundary_repeated_;// suit the .mesh file convention

	// info for FEM evaluation
	std::vector<std::pair<int, Eigen::Vector3d>> displacement_boundary_;// connector nodes id and zero displacement
	std::vector<int> cushion_surface_nodes_;
	// TODO noly need topology around surface_nodes rather than all nodes
	std::vector<std::vector<int>> faces_of_nodes_;// faces_of_nodes_[i] is a list of triangle faces on the surface adjacent to the node i
	std::vector<std::vector<int>> surface_nodes_of_nodes_;// surface_nodes_of_nodes_[i] is a list of nodes on the surface adjacent to the node i
};

