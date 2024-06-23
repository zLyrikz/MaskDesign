#include "TetrahedralMesh.h"
#include <Eigen/Geometry>
#include <igl/readMESH.h>
#include <igl/writeMESH.h>
#include <iostream>
#include <fstream>
#include <numeric>
#include <unordered_set>
//#include <unordered_map>
#include <map>
using std::cout;
using std::endl;
TetrahedralMesh::TetrahedralMesh()
{
}

TetrahedralMesh::~TetrahedralMesh()
{
}

void TetrahedralMesh::SetVertex(const Eigen::MatrixX3d& _vertex)
{
	if (vertex_.rows() != _vertex.rows())
	{
		cout << "WARNING!!!!!!!!!!!!!!!!!!!!!!!TetrahedralMesh::SetVertex vertex number changed when reset vertex" << endl;
	}
	vertex_ = _vertex;
}

void TetrahedralMesh::SetMesh(const std::string& _mesh_file)
{
	igl::readMESH(_mesh_file, vertex_, topology_tetrahedra_, topology_boundary_);
}

void TetrahedralMesh::SetMesh(const Eigen::MatrixX3d& _vertex, const Eigen::MatrixX4i& _topology_tetrahedra, const Eigen::MatrixX3i& _topology_boundary)
{
	vertex_ = _vertex;
	topology_tetrahedra_ = _topology_tetrahedra;
	topology_boundary_ = _topology_boundary;
}

void TetrahedralMesh::SetMesh_PossibleSameTopology(const Eigen::MatrixX3d& _vertex, const Eigen::MatrixX4i& _topology_tetrahedra, const Eigen::MatrixX3i& _topology_boundary)
{
	vertex_ = _vertex;
	if (topology_tetrahedra_.rows() != _topology_tetrahedra.rows())
	{
		topology_tetrahedra_ = _topology_tetrahedra;
	}
	if (topology_boundary_.rows() != _topology_boundary.rows())
	{
		topology_boundary_ = _topology_boundary;
	}
}

void TetrahedralMesh::SetMeshFromInpFile(const std::string& _inp_file)
{
	std::ifstream file(_inp_file);
	std::string line;

	bool node_section = false;
	bool element_section = false;
	std::vector<Eigen::Vector3d> vertex;
	std::vector<int> vertex_node_id;
	std::vector<Eigen::Vector4i> tetrahedra;

	if (file.is_open())
	{
		while (std::getline(file, line)) {
			if (line.find("*NODE") != std::string::npos) {
				node_section = true;
				element_section = false;
				continue;
			}
			if (line.find("*ELEMENT") != std::string::npos) {
				node_section = false;
				element_section = true;
				continue;
			}

			if (node_section) {
				if (line.empty()) {
					continue;
				}
				if (line[0] == '*')
				{
					node_section = false;
					continue;
				}
				std::istringstream iss(line);
				int node_id;
				double x, y, z;
				char comma;
				iss >> node_id >> comma >> x >> comma >> y >> comma >> z;
				vertex_node_id.push_back(node_id);
				vertex.push_back(Eigen::Vector3d(x, y, z));
			}

			if (element_section) {
				if (line.empty()) {
					continue;
				}
				if (line[0] == '*')
				{
					element_section = false;
					continue;
				}
				std::istringstream iss(line);
				int element_id, node1, node2, node3, node4;
				char comma;
				iss >> element_id >> comma >> node1 >> comma >> node2 >> comma >> node3 >> comma >> node4;

				auto iter_node1 = std::find(vertex_node_id.begin(), vertex_node_id.end(), node1);
				int index_node1 = std::distance(vertex_node_id.begin(), iter_node1);
				auto iter_node2 = std::find(vertex_node_id.begin(), vertex_node_id.end(), node2);
				int index_node2 = std::distance(vertex_node_id.begin(), iter_node2);
				auto iter_node3 = std::find(vertex_node_id.begin(), vertex_node_id.end(), node3);
				int index_node3 = std::distance(vertex_node_id.begin(), iter_node3);
				auto iter_node4 = std::find(vertex_node_id.begin(), vertex_node_id.end(), node4);
				int index_node4 = std::distance(vertex_node_id.begin(), iter_node4);
				tetrahedra.push_back(Eigen::Vector4i(index_node1, index_node2, index_node3, index_node4));
			}
		}
		file.close();
	}
	else
	{
		cout << "inp file not open" << endl;
	}
	vertex_.resize(vertex.size(), 3);
	for (int i = 0; i < vertex.size(); ++i) {
		vertex_.row(i) = vertex[i];
		
	}

	topology_tetrahedra_.resize(tetrahedra.size(), 4);
	for (int i = 0; i < tetrahedra.size(); ++i) {
		topology_tetrahedra_.row(i) = tetrahedra[i];
	}

	FindBoundaryTopology(topology_tetrahedra_, topology_boundary_, vertex_);

}

void TetrahedralMesh::WriteMesh(const std::string& _mesh_file) const
{
	igl::writeMESH(_mesh_file, vertex_, topology_tetrahedra_, topology_boundary_);
	// fix a bug in igl
	// write and "End" to the mesh file
	std::ofstream file(_mesh_file, std::ios::app);
	if (!file.is_open()) {
		std::cerr << "TetrahedralMesh::WriteMesh cannot open file" << std::endl;
		return;
	}
	file << "End" << std::endl;
	file.close();
}

void TetrahedralMesh::ComputeNormals()
{
	const int number_tetrahedron = topology_tetrahedra_.rows();
	tetrahedra_normal_.resize(number_tetrahedron);
	for (int iTetrahedron = 0; iTetrahedron < number_tetrahedron; ++iTetrahedron)
	{
		const Eigen::Vector3d& vertex0 = vertex_.row(topology_tetrahedra_(iTetrahedron, 0)).transpose();
		const Eigen::Vector3d& vertex1 = vertex_.row(topology_tetrahedra_(iTetrahedron, 1)).transpose();
		const Eigen::Vector3d& vertex2 = vertex_.row(topology_tetrahedra_(iTetrahedron, 2)).transpose();
		const Eigen::Vector3d& vertex3 = vertex_.row(topology_tetrahedra_(iTetrahedron, 3)).transpose();
		// opposite face of vertex 0
		// 123
		tetrahedra_normal_[iTetrahedron][0] = (vertex2 - vertex1).cross(vertex3 - vertex1).normalized();
		// opposite face of vertex 1
		// 203
		tetrahedra_normal_[iTetrahedron][1] = (vertex0 - vertex2).cross(vertex3 - vertex2).normalized();
		// opposite face of vertex 2
		// 130
		tetrahedra_normal_[iTetrahedron][2] = (vertex3 - vertex1).cross(vertex0 - vertex1).normalized();
		// opposite face of vertex 3
		// 102
		tetrahedra_normal_[iTetrahedron][3] = (vertex0 - vertex1).cross(vertex2 - vertex1).normalized();
	}

	const int number_boundary_face = topology_boundary_.rows();
	boundary_normal_.resize(number_boundary_face);
	for (int iBoundary = 0; iBoundary < number_boundary_face; ++iBoundary)
	{
		const Eigen::Vector3d& vertex0 = vertex_.row(topology_boundary_(iBoundary, 0)).transpose();
		const Eigen::Vector3d& vertex1 = vertex_.row(topology_boundary_(iBoundary, 1)).transpose();
		const Eigen::Vector3d& vertex2 = vertex_.row(topology_boundary_(iBoundary, 2)).transpose();
		boundary_normal_[iBoundary] = (vertex1 - vertex0).cross(vertex2 - vertex0).normalized();
	}
}

void TetrahedralMesh::ComputeSurfaceVertexNormal()
{
	boundary_vertex_normal_.clear();

	// a simple method, using OpenMesh
	ConstructSurfaceMesh();
	boundary_mesh_.request_face_normals();
	boundary_mesh_.request_vertex_normals();
	boundary_mesh_.update_face_normals();
	boundary_mesh_.update_vertex_normals();
	
	boundary_vertex_normal_.reserve(boundary_vertex_.size());
	for (int iVertex = 0; iVertex < boundary_vertex_.size(); ++iVertex)
	{
		const Mesh::Point& normal = boundary_mesh_.normal(boundary_mesh_.vertex_handle(iVertex));
		boundary_vertex_normal_.push_back(normal);
	}
}

void TetrahedralMesh::FindBoundaryVertex()
{
	// indicate wether the vertex is at boundary
	std::vector<bool> is_boundary_vertex(vertex_.rows(), false);
	for (int iFace = 0; iFace < topology_boundary_.rows(); ++iFace)
	{
		for (int iVertex = 0; iVertex < 3; ++iVertex)
		{
			is_boundary_vertex[topology_boundary_(iFace, iVertex)] = true;
		}
	}

	// collect the boundary vertex
	boundary_vertex_.clear();
	for (int iVertex = 0; iVertex < vertex_.rows(); ++iVertex)
	{
		if (is_boundary_vertex[iVertex])
		{
			boundary_vertex_.push_back(iVertex);
		}
	}
}

void TetrahedralMesh::FindBoundaryTopology(const Eigen::MatrixX4i& _topology_tetrahedra, Eigen::MatrixX3i& _topology_boundary, const Eigen::MatrixX3d& _vertex) const
{
	// idea from https://stackoverflow.com/questions/66607716/how-to-extract-surface-triangles-from-a-tetrahedral-mesh

	// std::unordered_set<int> is a face in the tetrahedron, 
	// bool is whether this face shows only once in all tetrahedra
	// std::vector<int> is the actual vertex order of the face in the tetrahedron
	// the following vectors have same size
	std::vector<std::unordered_set<int>> unordered_face_envelope;
	std::vector<std::vector<bool>> face_envelope;// each element is a vector of all nodes, indicating whether it has the nodes
	//std::vector<std::vector<int>> face_envelope_on_vertex;// element i-th is the vector of face index for vertex i
	std::vector<bool> face_is_boundary;
	std::vector<std::array<int, 3>> ordered_face_envelope;
	int count_boundary_face = 0;
	cout << "START find tetrahedral boundary...";
	for (int iTet = 0; iTet < _topology_tetrahedra.rows(); ++iTet) 
	{
		//cout << "iTet=" << iTet << endl;
		std::vector<std::vector<int>> faces{
			std::vector<int>{_topology_tetrahedra(iTet, 0), _topology_tetrahedra(iTet, 1), _topology_tetrahedra(iTet, 2)},
				std::vector<int>{_topology_tetrahedra(iTet, 1), _topology_tetrahedra(iTet, 0), _topology_tetrahedra(iTet, 3)},
				std::vector<int>{_topology_tetrahedra(iTet, 0), _topology_tetrahedra(iTet, 2), _topology_tetrahedra(iTet, 3)},
				std::vector<int>{_topology_tetrahedra(iTet, 2), _topology_tetrahedra(iTet, 1), _topology_tetrahedra(iTet, 3)}
		};

		int count_face = 0;
		for (const auto& iFace : faces)
		{
			//cout << "count_face=" << count_face << endl;
			std::unordered_set<int> unordered_face = {iFace[0], iFace[1], iFace[2]};
			//auto iterator = std::find(unordered_face_envelope.begin(), unordered_face_envelope.end(), unordered_face);

			// find the face
			int face_id = -1;
			for (int iAll = 0; iAll < face_envelope.size(); ++iAll)
			{
				if (face_envelope[iAll][iFace[0]] == true &&
					face_envelope[iAll][iFace[1]] == true &&
					face_envelope[iAll][iFace[2]] == true)
				{
					face_id = iAll;
					break;
				}
			}

			if (/*iterator != unordered_face_envelope.end()*/ face_id != -1)
			{
				// current face exists
				//int index = std::distance(unordered_face_envelope.begin(), iterator);
				int index = face_id;
				face_is_boundary[index] = false;
				count_boundary_face -= 1;// notice a face can only be repeatedly found at most twice, for it is shared by at most two tetrahedra
			}
			else
			{
				// current face encounter for the first time
				unordered_face_envelope.emplace_back(unordered_face);
				face_envelope.emplace_back(std::vector<bool>(_vertex.rows(), false));
				face_envelope.back()[iFace[0]] = true;
				face_envelope.back()[iFace[1]] = true;
				face_envelope.back()[iFace[2]] = true;

				face_is_boundary.emplace_back(true);
				ordered_face_envelope.emplace_back(std::array<int, 3>{iFace[0], iFace[1], iFace[2]});
				count_boundary_face += 1;
			}
			count_face++;
		}
	}
	cout << "done" << endl;

	// store the faces labeled true (twice as what seems in .mesh file)
	_topology_boundary.resize(2 * count_boundary_face, 3);
	int index = 0;
	for (int iFace = 0; iFace < face_is_boundary.size(); ++iFace)
	{
		if (face_is_boundary[iFace])
		{
			_topology_boundary(2 * index, 0) = ordered_face_envelope[iFace][0];
			_topology_boundary(2 * index, 1) = ordered_face_envelope[iFace][1];
			_topology_boundary(2 * index, 2) = ordered_face_envelope[iFace][2];

			_topology_boundary(2 * index + 1, 0) = ordered_face_envelope[iFace][2];
			_topology_boundary(2 * index + 1, 1) = ordered_face_envelope[iFace][1];
			_topology_boundary(2 * index + 1, 2) = ordered_face_envelope[iFace][0];

			++index;
		}
	}

}

void TetrahedralMesh::ConstructSurfaceMesh()
{
	boundary_mesh_.clear();
	// add vertex
	std::vector<Mesh::VertexHandle> vertex_handles;
	for (int iVertex = 0; iVertex < boundary_vertex_.size(); ++iVertex)
	{
		const Mesh::VertexHandle& vertex_handle = 
			boundary_mesh_.add_vertex(Mesh::Point(vertex_(boundary_vertex_[iVertex], 0), vertex_(boundary_vertex_[iVertex], 1), vertex_(boundary_vertex_[iVertex], 2)));
		vertex_handles.emplace_back(vertex_handle);
	}

	// inverse map, from vertex global id to boudary id
	std::map<int, int> boundary_vertex_global2local;
	for (int iVertex = 0; iVertex < boundary_vertex_.size(); ++iVertex)
	{
		boundary_vertex_global2local[boundary_vertex_[iVertex]] = iVertex;
	}

	// add face
	for (int iFace = 0; iFace < topology_boundary_.rows() / 2; ++iFace)
	{
		Mesh::VertexHandle vertex0 = vertex_handles[boundary_vertex_global2local[topology_boundary_(2 * iFace, 0)]];
		Mesh::VertexHandle vertex1 = vertex_handles[boundary_vertex_global2local[topology_boundary_(2 * iFace, 1)]];
		Mesh::VertexHandle vertex2 = vertex_handles[boundary_vertex_global2local[topology_boundary_(2 * iFace, 2)]];
		boundary_mesh_.add_face(vertex0, vertex1, vertex2);
	}
}

void TetrahedralMesh::Translate(const Eigen::Vector3d& _translate)
{
	for (int iVertex = 0; iVertex < vertex_.rows(); ++iVertex)
	{
		vertex_.row(iVertex) += _translate;
	}
}

void TetrahedralMesh::TranslateIndividually(const Eigen::VectorXd& _translate)
{
	int num_vertex = vertex_.rows();
	assert(num_vertex == _translate.rows() / 3);

	for (int iVertex = 0; iVertex < num_vertex; ++iVertex)
	{
		for (int iDimension = 0; iDimension < 3; ++iDimension)
		{
			vertex_(iVertex, iDimension) += _translate(3 * iVertex + iDimension);
		}
	}
}

void TetrahedralMesh::InitializeRenderInformation()
{
	is_vertex_show_.resize(vertex_.rows(), true);
	boundary_show_index_.resize(topology_boundary_.rows());
	std::iota(boundary_show_index_.begin(), boundary_show_index_.end(), 0);
	tetrahedra_show_index_.clear();
	ComputeNormals();
}

void TetrahedralMesh::SliceVertex(const Eigen::Vector3d& _plane_normal, const Eigen::Vector3d& _point_on_plane)
{
	assert(is_vertex_show_.size() == vertex_.rows());

	// check the vertex visible status
	for (int iVertex = 0; iVertex < vertex_.rows(); ++iVertex)
	{
		Eigen::Vector3d direction = vertex_.row(iVertex).transpose() - _point_on_plane;
		if (direction.dot(_plane_normal) < 0)
		{
			is_vertex_show_[iVertex] = false;
		}
		else
		{
			is_vertex_show_[iVertex] = true;
		}
	}

	// find the visible tetrahedra (instead of all vertices visible, we choose those on the plane)
	tetrahedra_show_index_.resize(topology_tetrahedra_.rows());
	int num_visible = 0;// record the number of visible tetrahedra 
	for (int iTetrahedra = 0; iTetrahedra < topology_tetrahedra_.rows(); ++iTetrahedra)
	{
		// check visibility (it has both visible and invisible vertex)
		bool is_visible = false;
		bool exist_visible_vertex = false;
		bool exist_invisible_vertex = false;
		for (int iVertex = 0; iVertex < 4; ++iVertex)
		{
			if (is_vertex_show_[topology_tetrahedra_(iTetrahedra, iVertex)])
			{
				exist_visible_vertex = true;
			}
			else
			{
				exist_invisible_vertex = true;
			}

			if (exist_visible_vertex && exist_invisible_vertex)
			{
				is_visible = true;
				break;
			}
		}

		if (is_visible)
		{
			assert(num_visible < topology_tetrahedra_.rows());

			tetrahedra_show_index_[num_visible] = iTetrahedra;
			++num_visible;
		}
	}
	tetrahedra_show_index_.resize(num_visible);

	// find the visible boundry (all vertices visible)
	boundary_show_index_.resize(topology_boundary_.rows());
	int num_visible_face = 0;// record the number of visible boundary
	for (int iBoundary = 0; iBoundary < topology_boundary_.rows(); ++iBoundary)
	{
		bool is_visible = true;
		// check visibility
		for (int iVertex = 0; iVertex < 3; ++iVertex)
		{
			if (is_vertex_show_[topology_boundary_(iBoundary, iVertex)] == false)
			{
				is_visible = false;
				break;
			}
		}

		if (is_visible)
		{
			assert(num_visible_face < topology_boundary_.rows());

			boundary_show_index_[num_visible_face] = iBoundary;
			++num_visible_face;
		}
	}
	boundary_show_index_.resize(num_visible_face);
}

int TetrahedralMesh::GetVertexNumber() const
{
	return vertex_.rows();
}

void TetrahedralMesh::GetVertex(int _index, Eigen::Vector3d& _vertex) const
{
	assert(_index < GetVertexNumber());
	_vertex = vertex_.row(_index);
}

const Eigen::MatrixX3d* TetrahedralMesh::GetVertex() const
{
	return &vertex_;
}

const Eigen::MatrixX4i* TetrahedralMesh::GetTopologyTetrahedra() const
{
	return &topology_tetrahedra_;
}

const Eigen::MatrixX3i* TetrahedralMesh::GetTopologyBoundary() const
{
	return &topology_boundary_;
}

void TetrahedralMesh::GetTetrahedron(std::array<Eigen::Vector3d, 4>& _vertex, int _element_id) const
{
	for (int iVertex = 0; iVertex < 4; ++iVertex)
	{
		_vertex[iVertex] = vertex_.row(topology_tetrahedra_(_element_id, iVertex))/*.transpose()*/;
	}
}

const std::vector<std::array<Eigen::Vector3d, 4>>* TetrahedralMesh::GetTetrahedraNormal() const
{
	return &tetrahedra_normal_;
}

const std::vector<Eigen::Vector3d>* TetrahedralMesh::GetBoundaryNormal() const
{
	return &boundary_normal_;
}

const std::vector<Mesh::Point>* TetrahedralMesh::GetBoundaryVertexNormal() const
{
	return &boundary_vertex_normal_;
}

const std::vector<int>* TetrahedralMesh::GetBoundaryVertex() const
{
	return &boundary_vertex_;
}

const std::vector<int>* TetrahedralMesh::GetVisibleTetrahedraId() const
{
	return &tetrahedra_show_index_;
}

const std::vector<int>* TetrahedralMesh::GetVisibleBoundaryFaceId() const
{
	return &boundary_show_index_;
}

const std::vector<bool>* TetrahedralMesh::GetVisibleVertexId() const
{
	return &is_vertex_show_;
}

void TetrahedralMesh::GetBoundaryMesh(Mesh& _surface_mesh) const
{
	//_surface_mesh.clear();
	//std::vector<Mesh::VertexHandle> vertices(boundary_vertex_.size());
	//std::map<int, int> boundary_vertex_global2local_id;// reversed map of the boundary_vertex_
	//for (int iVertex = 0; iVertex < boundary_vertex_.size(); ++iVertex)
	//{
	//	vertices[iVertex] = _surface_mesh.add_vertex
	//	(Mesh::Point(vertex_(boundary_vertex_[iVertex], 0), vertex_(boundary_vertex_[iVertex], 1), vertex_(boundary_vertex_[iVertex], 2)));

	//	boundary_vertex_global2local_id.insert(std::pair{boundary_vertex_[iVertex], iVertex});
	//}
	//
	//for (int iFace = 0; iFace < topology_boundary_.rows(); ++iFace)
	//{

	//	//_surface_mesh.add_face
	//	//	(vertices[boundary_vertex_global2local_id.at(topology_boundary_(iFace, 0))], 
	//	//	vertices[boundary_vertex_global2local_id.at(topology_boundary_(iFace, 1))], 
	//	//	vertices[boundary_vertex_global2local_id.at(topology_boundary_(iFace, 2))]);

	//	auto iter_vertex0 = std::find(boundary_vertex_.begin(), boundary_vertex_.end(), topology_boundary_(iFace, 0));
	//	auto iter_vertex1 = std::find(boundary_vertex_.begin(), boundary_vertex_.end(), topology_boundary_(iFace, 1));
	//	auto iter_vertex2 = std::find(boundary_vertex_.begin(), boundary_vertex_.end(), topology_boundary_(iFace, 2));
	//	int id0 = std::distance(boundary_vertex_.begin(), iter_vertex0);
	//	int id1 = std::distance(boundary_vertex_.begin(), iter_vertex1);
	//	int id2 = std::distance(boundary_vertex_.begin(), iter_vertex2);
	//	_surface_mesh.add_face
	//		(vertices[id0], 
	//		vertices[id1], 
	//		vertices[id2]);
	//}

	_surface_mesh.clear();
	std::vector<Mesh::VertexHandle> vertices(vertex_.rows());
	for (int iVertex = 0; iVertex < boundary_vertex_.size(); ++iVertex)
	{
		vertices[boundary_vertex_[iVertex]] = _surface_mesh.add_vertex
			(Mesh::Point(vertex_(boundary_vertex_[iVertex], 0), vertex_(boundary_vertex_[iVertex], 1), vertex_(boundary_vertex_[iVertex], 2)));
	}

	// in topology_boundary_, each triangle is repeated twice in different direction
	// e.g.
	// 5 3 2
	// 2 3 5
	// 3 1 2
	// 2 1 3
	for (int iFace = 0; iFace < topology_boundary_.rows() / 2; ++iFace)
	{

		_surface_mesh.add_face
			(vertices[topology_boundary_(2 * iFace + 1, 0)], 
			vertices[topology_boundary_(2 * iFace + 1, 1)],
			vertices[topology_boundary_(2 * iFace + 1, 2)]);

	}

	//cout << "topology_boundary_.rows()" << topology_boundary_.rows() << endl;
}

void TetrahedralMesh::GetVertexBehindPlane(const Eigen::Vector3d& _plane_normal, const Eigen::Vector3d& _point_on_plane,
	const std::vector<int>& _vertex_group, std::vector<int>& _selected) const
{
	for (const int& iVertex : _vertex_group)
	{
		Eigen::Vector3d direction = vertex_.row(iVertex).transpose() - _point_on_plane;
		if (direction.dot(_plane_normal) < 0)
		{
			// behind the plane
			_selected.push_back(iVertex);
		}
	}
}

void TetrahedralMesh::GetVertexNearMesh(const tmd::TriangleMeshDistance& _mesh_sdf,
	const std::vector<int>& _vertex_group, std::vector<int>& _selected, double _epsilon) const
{
	for (const int& iVertex : _vertex_group)
	{
		tmd::Result result = _mesh_sdf.signed_distance({ vertex_.row(iVertex)[0], vertex_.row(iVertex)[1], vertex_.row(iVertex)[2] });
		if (abs(result.distance) <= _epsilon)
		{
			// near the mesh with the given distance (epsilon)
			_selected.push_back(iVertex);
		}
	}
}

void TetrahedralMesh::MoveVertex(int _index, const Eigen::Vector3d& _position)
{
	assert(_index < GetVertexNumber());
	vertex_.row(_index) = _position;
}

void TetrahedralMesh::ChangeVertex(const Eigen::VectorXd& _vertex)
{
	if (_vertex.rows() == vertex_.rows() * 3)
	{
		for (int iNode = 0; iNode < vertex_.rows(); ++iNode)
		{
			vertex_(iNode, 0) = _vertex(3 * iNode + 0);
			vertex_(iNode, 1) = _vertex(3 * iNode + 1);
			vertex_(iNode, 2) = _vertex(3 * iNode + 2);
		}
	}
}
