#include "SubMesh.h"
#include <list>
using std::vector;
using std::list;
using std::cout;
using std::endl;

SubMesh::SubMesh()
{
	whole_mesh_ = nullptr;
	sub_mesh_ = nullptr;
}

SubMesh::~SubMesh()
{
	if (sub_mesh_ != nullptr)
	{
		delete sub_mesh_;
		sub_mesh_ = nullptr;
	}
}

void SubMesh::setData(const Mesh* _whole_mesh, std::vector<int>& _face_area)
{
	whole_mesh_ = _whole_mesh;
	face_area_ = _face_area;
}

void SubMesh::ConstructSubMesh()
{
	if (sub_mesh_ != nullptr)
	{
		delete sub_mesh_;
		sub_mesh_ = nullptr;
	}
	sub_mesh_ = new Mesh();

	// find all vertices in the choosen area
	vector<Mesh::VertexHandle> old_vertex_in_area;
	vector<Mesh::VertexHandle> new_vertex_in_area;
	for (int iFace = 0; iFace < face_area_.size(); ++iFace)
	{
		Mesh::FaceHandle face_handle = whole_mesh_->face_handle(face_area_[iFace]);
		vector<Mesh::VertexHandle> new_face;
		new_face.reserve(3);

		for (Mesh::ConstFaceVertexCCWIter iterVertex = whole_mesh_->cfv_ccwiter(face_handle); iterVertex.is_valid(); ++iterVertex)
		{
			// if this vertex is not stored in [vertex_in_area], store it.
			Mesh::VertexHandle current_vertex = *iterVertex;
			vector<Mesh::VertexHandle>::iterator iter_find = std::find(old_vertex_in_area.begin(), old_vertex_in_area.end(), current_vertex);
			if (iter_find == old_vertex_in_area.end())
			{		
				old_vertex_in_area.push_back(current_vertex);

				Mesh::VertexHandle new_vertex = sub_mesh_->add_vertex(whole_mesh_->point(current_vertex));
				new_vertex_in_area.push_back(new_vertex);

				new_face.push_back(new_vertex);
			}
			else // already added this vertex
			{
				int vertex_id = std::distance(old_vertex_in_area.begin(), iter_find);
				Mesh::VertexHandle new_vertex = new_vertex_in_area[vertex_id];
				new_face.push_back(new_vertex);
			}
		}

		sub_mesh_->add_face(new_face);
	}

	sub_mesh_->request_vertex_status();
	sub_mesh_->request_edge_status();
	sub_mesh_->request_face_status();

	sub_mesh_->request_face_normals();
	sub_mesh_->request_vertex_normals();

	sub_mesh_->update_face_normals();
	sub_mesh_->update_vertex_normals();
}

void SubMesh::Mesh2EigenMatrix(const Mesh& _mesh, Eigen::MatrixXd& _vertices, Eigen::MatrixXi& _faces) const
{
	int num_vertex = _mesh.n_vertices();
	int num_face = _mesh.n_faces();
	_vertices.resize(num_vertex, 3);
	_faces.resize(num_face, 3);
	for (int iVertex = 0; iVertex < num_vertex; ++iVertex)
	{
		const Mesh::Point& vertex = _mesh.point(_mesh.vertex_handle(iVertex));
		_vertices(iVertex, 0) = vertex[0];
		_vertices(iVertex, 1) = vertex[1];
		_vertices(iVertex, 2) = vertex[2];
	}
	for (int iFace = 0; iFace < num_face; ++iFace)
	{
		const Mesh::FaceHandle& face = _mesh.face_handle(iFace);

		int iVertex = 0;
		for (Mesh::VertexHandle& iVertexHandle :_mesh.fv_ccw_range(face))
		{
			if (iVertex == 3)
			{
				cout << "[WARNING] SubMesh::Mesh2EigenMatrix fail, face not trianlge!" << endl;
			}
			_faces(iFace, iVertex) = iVertexHandle.idx();
			++iVertex;
		}
		assert(iVertex == 3);// triangle face has 3 vertices
	}
}

void SubMesh::MeshFromEigenMatrix(Mesh& _mesh, const Eigen::MatrixXd& _vertices, const Eigen::MatrixXi& _faces) const
{
	_mesh.clear();
	vector<Mesh::VertexHandle> vertices(_vertices.rows());
	for (int iVertex = 0; iVertex < _vertices.rows(); ++iVertex)
	{
		vertices[iVertex] = _mesh.add_vertex
			(Mesh::Point(_vertices(iVertex, 0), _vertices(iVertex, 1), _vertices(iVertex, 2)));
	}
	for (int iFace = 0; iFace < _faces.rows(); ++iFace)
	{
		_mesh.add_face
			(vertices[_faces(iFace, 0)], vertices[_faces(iFace, 1)], vertices[_faces(iFace, 2)]);
	}
}

void SubMesh::MeshFromVerticesAndTriangles(Mesh& _result_mesh,
	const std::vector<std::array<double, 3>>& _vertices, const std::vector<std::array<int, 3>>& _triangles) const
{
	_result_mesh = Mesh();
	vector<Mesh::VertexHandle> vertices(_vertices.size());
	for (int iVertex = 0; iVertex < _vertices.size(); ++iVertex)
	{
		vertices[iVertex] = _result_mesh.add_vertex
			(Mesh::Point(_vertices[iVertex][0], _vertices[iVertex][1], _vertices[iVertex][2]));
	}
	for (int iFace = 0; iFace < _triangles.size(); ++iFace)
	{
		_result_mesh.add_face
			(vertices[_triangles[iFace][0]], vertices[_triangles[iFace][1]], vertices[_triangles[iFace][2]]);
	}
}

void SubMesh::MeshFromUVSurface(const std::vector<std::vector<Mesh::Point>>& _uv_points, Mesh& _surface_mesh)
{
	int num_u = _uv_points.size();
	if (num_u > 0)
	{
		_surface_mesh.clear();
		int num_v = _uv_points[0].size();

		std::vector<std::vector<Mesh::VertexHandle>> uv_points;
		//uv_points.reserve(num_u);
		//std::vector<Mesh::VertexHandle> v_points(num_v);

		//for (int iU = 0; iU < num_u; ++iU)
		//{
		//	for (int iV = 0; iV < num_v; ++iV)
		//	{
		//		v_points[iV] = _surface_mesh.add_vertex(_uv_points[iU][iV]);
		//	}
		//	uv_points.push_back(v_points);
		//}

		MeshAddVertexFromUVPoints(_uv_points, _surface_mesh, uv_points);

		MeshAddFaceFromUVPoints(_surface_mesh, uv_points);
		//for (int iU = 0; iU < num_u - 1; ++iU)
		//{
		//	for (int iV = 0; iV < num_v - 1; ++iV)
		//	{
		//		_surface_mesh.add_face(uv_points[iU][iV], uv_points[iU + 1][iV + 1], uv_points[iU][iV + 1]);
		//		_surface_mesh.add_face(uv_points[iU][iV], uv_points[iU + 1][iV], uv_points[iU + 1][iV + 1]);
		//	}
		//}
	}
}

void SubMesh::MeshFromUVSurface_PeriodU(const std::vector<std::vector<Mesh::Point>>& _uv_points, Mesh& _surface_mesh)
{
	int num_u = _uv_points.size();
	if (num_u > 0)
	{
		_surface_mesh.clear();
		int num_v = _uv_points[0].size();

		std::vector<std::vector<Mesh::VertexHandle>> uv_points;
		//uv_points.reserve(num_u);
		//std::vector<Mesh::VertexHandle> v_points(num_v);

		//for (int iU = 0; iU < num_u; ++iU)
		//{
		//	for (int iV = 0; iV < num_v; ++iV)
		//	{
		//		v_points[iV] = _surface_mesh.add_vertex(_uv_points[iU][iV]);
		//	}
		//	uv_points.push_back(v_points);
		//}

		//for (int iU = 0; iU < num_u - 1; ++iU)
		//{
		//	for (int iV = 0; iV < num_v - 1; ++iV)
		//	{
		//		_surface_mesh.add_face(uv_points[iU][iV], uv_points[iU + 1][iV + 1], uv_points[iU][iV + 1]);
		//		_surface_mesh.add_face(uv_points[iU][iV], uv_points[iU + 1][iV], uv_points[iU + 1][iV + 1]);
		//	}
		//}
		MeshAddVertexFromUVPoints(_uv_points, _surface_mesh, uv_points);

		MeshAddFaceFromUVPoints_PeriodicU(_surface_mesh, uv_points);
		//// period case, connect last with first
		//for (int iV = 0; iV < num_v - 1; ++iV)
		//{
		//	_surface_mesh.add_face(uv_points[num_u - 1][iV], uv_points[0][iV + 1], uv_points[num_u - 1][iV + 1]);
		//	_surface_mesh.add_face(uv_points[num_u - 1][iV], uv_points[0][iV], uv_points[0][iV + 1]);
		//}
	}
}

void SubMesh::MeshFromUVSurface_PeriodUV(const std::vector<std::vector<Mesh::Point>>& _uv_points, Mesh& _surface_mesh, bool _flip_normal)
{
	int num_u = _uv_points.size();
	if (num_u > 0)
	{
		_surface_mesh.clear();
		int num_v = _uv_points[0].size();

		std::vector<std::vector<Mesh::VertexHandle>> uv_points;
		MeshAddVertexFromUVPoints(_uv_points, _surface_mesh, uv_points);

		if (_flip_normal)
		{
			MeshAddFaceFromUVPoints_PeriodicUV_FlipNormal(_surface_mesh, uv_points);
		}
		else
		{
			MeshAddFaceFromUVPoints_PeriodicUV(_surface_mesh, uv_points);
		}
	}
}

void SubMesh::CuboidTubeMeshFrom4Polylines_Period(const std::vector<Mesh::Point>& _polyline1,
	const std::vector<Mesh::Point>& _polyline2, const std::vector<Mesh::Point>& _polyline3, const std::vector<Mesh::Point>& _polyline4, 
	Mesh& _result)
{
	assert(_polyline1.size() == _polyline2.size());
	assert(_polyline3.size() == _polyline2.size());
	assert(_polyline4.size() == _polyline3.size());

	int num_point = _polyline1.size();
	std::vector<Mesh::VertexHandle> polyline1;//vertex handle for the result
	std::vector<Mesh::VertexHandle> polyline2;
	std::vector<Mesh::VertexHandle> polyline3;
	std::vector<Mesh::VertexHandle> polyline4;
	polyline1.reserve(num_point);
	polyline2.reserve(num_point);
	polyline3.reserve(num_point);
	polyline4.reserve(num_point);
	for (int iPoint = 0; iPoint < num_point; ++iPoint)
	{
		Mesh::VertexHandle vertex1 = _result.add_vertex(_polyline1[iPoint]);
		Mesh::VertexHandle vertex2 = _result.add_vertex(_polyline2[iPoint]);
		Mesh::VertexHandle vertex3 = _result.add_vertex(_polyline3[iPoint]);
		Mesh::VertexHandle vertex4 = _result.add_vertex(_polyline4[iPoint]);
		polyline1.push_back(vertex1);
		polyline2.push_back(vertex2);
		polyline3.push_back(vertex3);
		polyline4.push_back(vertex4);
	}
	for (int iPoint = 0; iPoint < num_point; ++iPoint)
	{
		int next_point_id = (iPoint + 1) % num_point;

		_result.add_face(polyline1[iPoint], polyline1[next_point_id], polyline2[iPoint]);
		_result.add_face(polyline2[iPoint], polyline1[next_point_id], polyline2[next_point_id]);

		_result.add_face(polyline1[next_point_id], polyline1[iPoint], polyline4[next_point_id]);
		_result.add_face(polyline4[next_point_id], polyline1[iPoint], polyline4[iPoint]);

		_result.add_face(polyline3[next_point_id], polyline2[iPoint], polyline2[next_point_id]);
		_result.add_face(polyline2[iPoint], polyline3[next_point_id], polyline3[iPoint]);

		_result.add_face(polyline4[next_point_id], polyline3[iPoint], polyline3[next_point_id]);
		_result.add_face(polyline3[iPoint], polyline4[next_point_id], polyline4[iPoint]);
	}
}

void SubMesh::CylinricalMeshFromPolylines_Period(const std::vector<const std::vector<Mesh::Point>*>& _polylines, Mesh& _result)
{
	if (_polylines.size() > 0)
	{
		int num_points = _polylines[0]->size();

		// sanity check
		for (int iPolyline = 1; iPolyline < _polylines.size(); ++iPolyline)
		{
			if (_polylines[iPolyline]->size() != num_points)
			{
				cout << "ERROR FROM SubMesh::CylinricalMeshFromPolylines_Period polyline points not same" << endl;
				return;
			}
		}

		std::vector<std::vector<Mesh::VertexHandle>> polylines(_polylines.size(), std::vector<Mesh::VertexHandle>(num_points));//vertex handle for the result
		for (int iPoint = 0; iPoint < num_points; ++iPoint)
		{
			for (int iPolyline = 0; iPolyline < _polylines.size(); ++iPolyline)
			{
				polylines[iPolyline][iPoint] = _result.add_vertex(_polylines[iPolyline]->at(iPoint));
				//cout << _polylines[iPolyline]->at(iPoint) << endl;
			}
		}

		for (int iPoint = 0; iPoint < num_points; ++iPoint)
		{
			int next_point_id = (iPoint + 1) % num_points;

			for (int iPolyline = 0; iPolyline < _polylines.size(); ++iPolyline)
			{
				int next_polyline_id = (iPolyline + 1) % _polylines.size();

				_result.add_face(polylines[iPolyline][iPoint], polylines[next_polyline_id][iPoint], polylines[iPolyline][next_point_id]);
				_result.add_face(polylines[next_polyline_id][iPoint], polylines[next_polyline_id][next_point_id], polylines[iPolyline][next_point_id]);

			}
		}
		
	}
}

void SubMesh::PointCloudMeshFromPoints(const std::vector<Mesh::Point>& _vertices, Mesh& _point_cloud) const
{
	//_point_cloud.clear();
	for (int iPoint = 0; iPoint < _vertices.size(); ++iPoint)
	{
		_point_cloud.add_vertex(_vertices[iPoint]);
	}
}

void SubMesh::PointCloudMeshFromPoints(const std::vector<Eigen::Vector3d>& _vertices, Mesh& _point_cloud) const
{
	for (int iPoint = 0; iPoint < _vertices.size(); ++iPoint)
	{
		Mesh::Point point(_vertices[iPoint][0], _vertices[iPoint][1], _vertices[iPoint][2]);
		_point_cloud.add_vertex(point);
	}
}
void SubMesh::PointCloudMeshToEigenVectorX(const Mesh& _point_cloud, std::vector<Eigen::VectorXd>& _points) const
{
	_points.clear();
	int num_point = _point_cloud.n_vertices();
	_points.reserve(num_point);
	for (Mesh::VertexHandle iVertex : _point_cloud.vertices())
	{
		Mesh::Point point = _point_cloud.point(iVertex);
		_points.push_back(Eigen::Map<Eigen::VectorXd>(point.data(), 3));
	}
}

void SubMesh::PointCloudMeshToPoints(const Mesh& _point_cloud, std::vector<Mesh::Point>& _points) const
{
	_points.clear();
	int num_point = _point_cloud.n_vertices();
	_points.reserve(num_point);
	for (Mesh::VertexHandle iVertex : _point_cloud.vertices())
	{
		_points.push_back(_point_cloud.point(iVertex));
	}
}

void SubMesh::MergeMesh(Mesh& _main, const Mesh& _merge, std::vector<Mesh::VertexHandle>& _added_vertex) const
{
	// make the added vertex index in this vector same as the vertex id of merge mesh (for convenience of adding faces)
	_added_vertex.resize(_merge.n_vertices());

	// add vertices
	for (const Mesh::VertexHandle& iVertex : _merge.vertices())
	{
		Mesh::VertexHandle new_vertex = _main.add_vertex(_merge.point(iVertex));
		_added_vertex[iVertex.idx()] = new_vertex;
	}

	// add faces
	for (const Mesh::FaceHandle& iFace : _merge.faces())
	{
		std::vector<Mesh::VertexHandle> face;
		face.reserve(3);
		for (const Mesh::VertexHandle& iVertex : _merge.fv_ccw_range(iFace))
		{
			// thanks to the same index we kept
			face.push_back(_added_vertex[iVertex.idx()]);
		}
		_main.add_face(face);
	}
}

void SubMesh::MeshAddVertexFromUVPoints(const std::vector<std::vector<Mesh::Point>>& _uv_points, Mesh& _surface_mesh, std::vector<std::vector<Mesh::VertexHandle>>& _uv_handle) const
{
	assert(!_uv_points.empty());
	assert(_uv_handle.empty());

	int num_u = _uv_points.size();
	int num_v = _uv_points[0].size();

	_uv_handle.reserve(num_u);
	std::vector<Mesh::VertexHandle> v_points(num_v);

	for (int iU = 0; iU < num_u; ++iU)
	{
		for (int iV = 0; iV < num_v; ++iV)
		{
			v_points[iV] = _surface_mesh.add_vertex(_uv_points[iU][iV]);
		}
		_uv_handle.push_back(v_points);
	}
}
void SubMesh::MeshAddFaceFromUVPoints(Mesh& _surface_mesh, const std::vector<std::vector<Mesh::VertexHandle>>& _uv_handle) const
{
	for (int iU = 0; iU < _uv_handle.size() - 1; ++iU)
	{
		for (int iV = 0; iV < _uv_handle[0].size() - 1; ++iV)
		{
			_surface_mesh.add_face(_uv_handle[iU][iV], _uv_handle[iU + 1][iV + 1], _uv_handle[iU][iV + 1]);
			_surface_mesh.add_face(_uv_handle[iU][iV], _uv_handle[iU + 1][iV], _uv_handle[iU + 1][iV + 1]);
		}
	}
}
void SubMesh::MeshAddFaceFromUVPoints_PeriodicU(Mesh& _surface_mesh, const std::vector<std::vector<Mesh::VertexHandle>>& _uv_handle) const
{
	MeshAddFaceFromUVPoints(_surface_mesh, _uv_handle);

	assert(!_uv_handle.empty());
	int num_u = _uv_handle.size();
	int num_v = _uv_handle[0].size();
	// period case, connect last with first
	if (num_u >= 1)
	{
		for (int iV = 0; iV < num_v - 1; ++iV)
		{
			_surface_mesh.add_face(_uv_handle[num_u - 1][iV], _uv_handle[0][iV + 1], _uv_handle[num_u - 1][iV + 1]);
			_surface_mesh.add_face(_uv_handle[num_u - 1][iV], _uv_handle[0][iV], _uv_handle[0][iV + 1]);
		}
	}

}
void SubMesh::MeshAddFaceFromUVPoints_PeriodicUV(Mesh& _surface_mesh, const std::vector<std::vector<Mesh::VertexHandle>>& _uv_handle) const
{
	MeshAddFaceFromUVPoints_PeriodicU(_surface_mesh, _uv_handle);

	assert(!_uv_handle.empty());
	int num_u = _uv_handle.size();
	int num_v = _uv_handle[0].size();

	if (num_u >= 1 && num_v >= 1)
	{
		for (int iU = 0; iU < num_u - 1; ++iU)
		{
			_surface_mesh.add_face(_uv_handle[iU][num_v - 1], _uv_handle[iU + 1][0], _uv_handle[iU][0]);
			_surface_mesh.add_face(_uv_handle[iU][num_v - 1], _uv_handle[iU + 1][num_v - 1], _uv_handle[iU + 1][0]);
		}

		_surface_mesh.add_face(_uv_handle[num_u - 1][num_v - 1], _uv_handle[0][0], _uv_handle[num_u - 1][0]);
		_surface_mesh.add_face(_uv_handle[num_u - 1][num_v - 1], _uv_handle[0][num_v - 1], _uv_handle[0][0]);
	}

}
void SubMesh::MeshAddFaceFromUVPoints_FlipNormal(Mesh& _surface_mesh, const std::vector<std::vector<Mesh::VertexHandle>>& _uv_handle) const
{
	for (int iU = 0; iU < _uv_handle.size() - 1; ++iU)
	{
		for (int iV = 0; iV < _uv_handle[0].size() - 1; ++iV)
		{
			_surface_mesh.add_face(_uv_handle[iU][iV], _uv_handle[iU][iV + 1], _uv_handle[iU + 1][iV + 1]);
			_surface_mesh.add_face(_uv_handle[iU][iV], _uv_handle[iU + 1][iV + 1], _uv_handle[iU + 1][iV]);
		}
	}
}
void SubMesh::MeshAddFaceFromUVPoints_PeriodicUV_FlipNormal(Mesh& _surface_mesh, const std::vector<std::vector<Mesh::VertexHandle>>& _uv_handle) const
{
	MeshAddFaceFromUVPoints_FlipNormal(_surface_mesh, _uv_handle);

	assert(!_uv_handle.empty());
	int num_u = _uv_handle.size();
	int num_v = _uv_handle[0].size();

	// period case, connect last with first
	if (num_u >= 1 && num_v >= 1)
	{
		for (int iV = 0; iV < num_v - 1; ++iV)
		{
			_surface_mesh.add_face(_uv_handle[num_u - 1][iV], _uv_handle[num_u - 1][iV + 1], _uv_handle[0][iV + 1]);
			_surface_mesh.add_face(_uv_handle[num_u - 1][iV], _uv_handle[0][iV + 1], _uv_handle[0][iV]);
		}
		for (int iU = 0; iU < num_u - 1; ++iU)
		{
			_surface_mesh.add_face(_uv_handle[iU][num_v - 1], _uv_handle[iU][0], _uv_handle[iU + 1][0]);
			_surface_mesh.add_face(_uv_handle[iU][num_v - 1], _uv_handle[iU + 1][0], _uv_handle[iU + 1][num_v - 1]);
		}
		_surface_mesh.add_face(_uv_handle[num_u - 1][num_v - 1], _uv_handle[num_u - 1][0], _uv_handle[0][0]);
		_surface_mesh.add_face(_uv_handle[num_u - 1][num_v - 1], _uv_handle[0][0], _uv_handle[0][num_v - 1]);
	}
}
//*/

void SubMesh::getSubMesh(Mesh& _sub_mesh) const
{
	if (sub_mesh_ != nullptr)
	{
		_sub_mesh = *sub_mesh_;
	}
}
