#include "MeshOffset.h"
#include "MeshInfo.h"
#include "SubMesh.h"
#include "../Utility/VectorOperations.h"

#include <igl/signed_distance.h>
#include <igl/offset_surface.h>

#include <iostream>
using std::cout;
using std::endl;
MeshOffset::MeshOffset()
{
}

MeshOffset::~MeshOffset()
{
}

void MeshOffset::OffsetSurface(const Mesh& _original, Mesh& _offset, double _distance) const
{
	_offset = _original;
	ComputeMeshNormal(_offset);

	for (OpenMesh::SmartVertexHandle iVertex : _offset.vertices())
	{
		Mesh::Point normal = _offset.normal(iVertex);
		Mesh::Point point = _offset.point(iVertex);
		_offset.set_point(iVertex, point + _distance * normal);

		//cout << "normal\n" << normal << endl;
		//cout << "normal.norm()" << normal.norm() << endl;
		assert(abs(normal.norm()-1) < 1e-10);
	}

	//FlipNormal(_offset);
}

void MeshOffset::OffsetSurface(const std::vector<std::vector<Mesh::Point>>& _vu_points, std::vector<std::vector<Mesh::Point>>& _offset_vu_points, double _distance) const
{
	int num_v = _vu_points.size();
	if (num_v > 0)
	{
		// we will get the offset direction from the mesh normal
		// so we construct a mesh from the original uv points
		// we need the natural vertex index correspoindence from the construction
		Mesh original;
		SubMesh construct_mesh;
		construct_mesh.MeshFromUVSurface(_vu_points, original);
		ComputeMeshNormal(original);

		int num_u = _vu_points[0].size();
		_offset_vu_points.resize(num_v, std::vector<Mesh::Point>(num_u));
		assert(num_u * num_v == original.n_vertices());
		for (int iVertex = 0; iVertex < original.n_vertices(); ++iVertex)
		{
			//iVertex = v_idx*num_u + u_idx; (v_idx, u_idx)
			int u_idx = iVertex % num_u;
			int v_idx = iVertex / num_u;

			Mesh::Point point = original.point(original.vertex_handle(iVertex));
			assert((_vu_points[v_idx][u_idx] - point).sqrnorm() < 1e-10);// make sure the correspondence is correct

			Mesh::Point normal = original.normal(original.vertex_handle(iVertex));

			_offset_vu_points[v_idx][u_idx] = point + _distance * normal;
		}
	}
}

void MeshOffset::ThickenSurface(const Mesh& _original, Mesh& _watertight, double _distance) const
{
	_watertight = _original;
	ComputeMeshNormal(_watertight);

	//before offset, get initial boundary halfedge for traverse on boundary later
	std::vector<Mesh::HalfedgeHandle> initial_boundaries;
	MeshInfo find_boundary(&_watertight);
	find_boundary.getInitialBoundaryHalfedge(initial_boundaries);

	// add offset surface
	// traverse faces of the original surface, add offset vertex and record the face vertex in clockwise order
	// after the construction, offset_vertex[i] is the offset of vertex i
	std::vector<Mesh::VertexHandle> offset_vertex(_watertight.n_vertices());
	std::vector<bool> is_vertex_added(_watertight.n_vertices(), false);
	std::vector<std::vector<Mesh::VertexHandle>> reoriented_faces;
	reoriented_faces.reserve(_watertight.n_faces());
	for (OpenMesh::FaceHandle iFace : _watertight.faces())
	{
		std::vector<Mesh::VertexHandle> reoriented_face;// clockwise
		reoriented_face.reserve(3);
		for (Mesh::VertexHandle iVertex : _watertight.fv_cw_range(iFace))
		{
			int vertex_id = iVertex.idx();// OpenMesh must not change previous vertex idx if we add new vertex to the mesh
			assert(vertex_id < _original.n_vertices());
			if (vertex_id >= _original.n_vertices())
			{
				cout << "[WARNING] from MeshOffset::ThickenSurface unexpected behavior" << endl;
			}
			// add new offset vertex if not added
			if (is_vertex_added[vertex_id] == false)
			{
				Mesh::Point normal = _watertight.normal(iVertex);
				Mesh::Point point = _watertight.point(iVertex);
				Mesh::VertexHandle new_vertex = _watertight.add_vertex(point + _distance * normal);
				offset_vertex[vertex_id] = new_vertex;

				is_vertex_added[vertex_id] = true;
			}

			reoriented_face.push_back(offset_vertex[vertex_id]);
		}
		reoriented_faces.push_back(reoriented_face);
	}
	// add new faces
	for (auto& iNewFace : reoriented_faces)
	{
		_watertight.add_face(iNewFace);
	}

	// connect the boundary of offset and the original
	std::vector<std::vector<Mesh::VertexHandle>> connector_faces;
	// get the connecting faces
	for (Mesh::HalfedgeHandle& iInitialBoundary : initial_boundaries)
	{
		{
			Mesh::HalfedgeHandle current_boundary = iInitialBoundary;
			do
			{
				Mesh::VertexHandle vertex_from = _watertight.from_vertex_handle(current_boundary);
				Mesh::VertexHandle vertex_to = _watertight.to_vertex_handle(current_boundary);
				int from_id = vertex_from.idx();
				int to_id = vertex_to.idx();
				Mesh::VertexHandle offset_from = offset_vertex[from_id];
				Mesh::VertexHandle offset_to = offset_vertex[to_id];

				std::vector<Mesh::VertexHandle> triangle1{ vertex_from , vertex_to, offset_to };
				std::vector<Mesh::VertexHandle> triangle2{ offset_to , offset_from, vertex_from };
				connector_faces.push_back(triangle1);
				connector_faces.push_back(triangle2);

				// update for next loop
				current_boundary = _watertight.next_halfedge_handle(current_boundary);
			} while (current_boundary != iInitialBoundary);

		}
	}
	// add the connecting faces
	for (std::vector<Mesh::VertexHandle>& iFace : connector_faces)
	{
		_watertight.add_face(iFace);
	}
}

void MeshOffset::ThickenSurface_Libigl(const Mesh& _original, Mesh& _watertight, double _distance) const
{
	Eigen::MatrixXd V, SV;
	Eigen::MatrixXi F, SF;
	Eigen::MatrixXd GV;
	Eigen::VectorXd S;
	Eigen::RowVector3i side;

	double isolevel = _distance;
	int s = 200;
	igl::SignedDistanceType signed_distance_type = igl::SignedDistanceType::SIGNED_DISTANCE_TYPE_UNSIGNED;

	SubMesh convert_mesh_type;
	convert_mesh_type.Mesh2EigenMatrix(_original, V, F);

	igl::offset_surface(V, F, isolevel, s, signed_distance_type, SV, SF, GV, side, S);
	convert_mesh_type.MeshFromEigenMatrix(_watertight, SV, SF);
}

void MeshOffset::ThickenSurface(const std::vector<std::vector<Mesh::Point>>& _vu_points, std::vector<std::vector<Mesh::Point>>& _merged_vu_points, double _distance) const
{
	std::vector<std::vector<Mesh::Point>> offset_vu_points;
	OffsetSurface(_vu_points, offset_vu_points, _distance);
	_merged_vu_points = _vu_points;
	for (int iV = 0; iV < offset_vu_points.size(); ++iV)
	{
		VectorOperations append_points;
		append_points.AppendBack(_merged_vu_points[iV], offset_vu_points[iV]);
	}
}

void MeshOffset::WatertightFromOffset(const Mesh& _original, const Mesh& _offset, Mesh& _watertight) const
{
	_watertight = _original;

	std::vector<Mesh::HalfedgeHandle> initial_boundaries;
	MeshInfo find_boundary(&_original);
	find_boundary.getInitialBoundaryHalfedge(initial_boundaries);


	std::vector<Mesh::VertexHandle> added_vertex;
	SubMesh merge_it;
	merge_it.MergeMesh(_watertight, _offset, added_vertex);

	std::vector<std::vector<Mesh::VertexHandle>> connector_faces;
	// get the connecting faces
	for (Mesh::HalfedgeHandle& iInitialBoundary : initial_boundaries)
	{
		{
			Mesh::HalfedgeHandle current_boundary = iInitialBoundary;
			do
			{
				Mesh::VertexHandle vertex_from = _watertight.from_vertex_handle(_watertight.halfedge_handle(current_boundary.idx()));
				Mesh::VertexHandle vertex_to = _watertight.to_vertex_handle(_watertight.halfedge_handle(current_boundary.idx()));
				int from_id = vertex_from.idx();
				int to_id = vertex_to.idx();
				Mesh::VertexHandle offset_from = added_vertex[from_id];
				Mesh::VertexHandle offset_to = added_vertex[to_id];

				std::vector<Mesh::VertexHandle> triangle1{ vertex_from , vertex_to, offset_to };
				std::vector<Mesh::VertexHandle> triangle2{ offset_to , offset_from, vertex_from };
				connector_faces.push_back(triangle1);
				connector_faces.push_back(triangle2);

				// update for next loop
				current_boundary = _watertight.next_halfedge_handle(current_boundary);
			} while (current_boundary != iInitialBoundary);

		}
	}
	// add the connecting faces
	for (std::vector<Mesh::VertexHandle>& iFace : connector_faces)
	{
		_watertight.add_face(iFace);
	}

}

void MeshOffset::FlipNormal(Mesh& _mesh) const
{
	Mesh fliped;
	// keep vertex index the same
	for (int iVertex = 0; iVertex < _mesh.n_vertices(); ++iVertex)
	{
		fliped.add_vertex(_mesh.point(_mesh.vertex_handle(iVertex)));
	}
	std::vector<std::vector<Mesh::VertexHandle>> reoriented_faces;
	reoriented_faces.reserve(_mesh.n_faces());
	for (OpenMesh::SmartFaceHandle iFace : _mesh.faces())
	{
		std::vector<Mesh::VertexHandle> reoriented_face;
		reoriented_face.reserve(3);
		for (Mesh::VertexHandle iVertex : iFace.vertices_cw())
		{
			reoriented_face.push_back(fliped.vertex_handle(iVertex.idx()));// this works because vertex index is the same
		}
		reoriented_faces.push_back(reoriented_face);
	}
	for (auto& iFace : reoriented_faces)
	{
		fliped.add_face(iFace);
	}
	_mesh = fliped;
}

void MeshOffset::ComputeMeshNormal(Mesh& _mesh) const
{
	if (!_mesh.has_vertex_normals())
	{
		//_mesh.request_vertex_status();

		_mesh.request_face_normals();
		_mesh.request_vertex_normals();
		_mesh.update_face_normals();
		_mesh.update_vertex_normals();

		//_mesh.update_vertex_normals();
	}
}
