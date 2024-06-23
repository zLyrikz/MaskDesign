#include "MeshInfo.h"
#include <vector>
#include <iostream>
using std::cout;
using std::endl;

MeshInfo::MeshInfo(const Mesh* _mesh)
{
	mesh_ = _mesh;
}

void MeshInfo::getBoundingBox(Mesh::Point& _bbMin, Mesh::Point& _bbMax, bool _add_epsilon) const
{
	Mesh::VertexIter vIt = mesh_->vertices_begin();
	Mesh::VertexIter vEnd = mesh_->vertices_end();
	_bbMin = _bbMax = OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_->point(vIt));

	for (; vIt != vEnd; ++vIt)
	{
		_bbMin.minimize(OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_->point(vIt)));
		_bbMax.maximize(OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh_->point(vIt)));
	}

	if (_add_epsilon)
	{
		_bbMin -= Mesh::Point(1e-5);
		_bbMax += Mesh::Point(1e-5);
	}
}

double MeshInfo::getFaceMaxCoordinate(int _face_id, int _coordinate) const
{
	double max = -DBL_MAX;
	Mesh::Point point;

	for (auto iVertex : mesh_->fv_range(mesh_->face_handle(_face_id)))
	{
		point = mesh_->point(iVertex);
		if (max < point[_coordinate])
		{
			max = point[_coordinate];
		}
	}
	return max;
}

double MeshInfo::getFaceMinCoordinate(int _face_id, int _coordinate) const
{
	double min = DBL_MAX;
	Mesh::Point point;

	for (auto iVertex : mesh_->fv_range(mesh_->face_handle(_face_id)))
	{
		point = mesh_->point(iVertex);
		if (min > point[_coordinate])
		{
			min = point[_coordinate];
		}
	}
	return min;
}

double MeshInfo::getAverageEdgeLength() const
{
	double average_length = 0.0;
	double num_edge = 0.0;// to be a divisor
	for (auto iEdge : mesh_->edges())
	{
		average_length += (mesh_->point(iEdge.v0()) - mesh_->point(iEdge.v1())).norm();
		num_edge += 1.0;
	}
	average_length /= num_edge;

	return average_length;
}

void MeshInfo::getVerticesAndTriangles(std::vector<std::array<double, 3>>& _vertices, std::vector<std::array<int, 3>>& _triangles) const
{
	_vertices.resize(mesh_->n_vertices());
	_triangles.resize(mesh_->n_faces());
	std::vector<bool> is_vertex_added(mesh_->n_vertices(), false);
	{
		int count_face = 0;
		for (auto iFace : mesh_->faces())
		{
			int count = 0;
			for (auto iVertex : iFace.vertices())
			{
				_triangles[count_face][count] = iVertex.idx();
				
				if (is_vertex_added[iVertex.idx()] == false)
				{
					Mesh::Point point = mesh_->point(iVertex);
					_vertices[iVertex.idx()][0] = point[0];
					_vertices[iVertex.idx()][1] = point[1];
					_vertices[iVertex.idx()][2] = point[2];
					is_vertex_added[iVertex.idx()] = true;
				}
				++count;
			}
			++count_face;
		}
	}

	//int count_t = 0;
	//for (auto iT : _triangles)
	//{
	//	cout << "trinagle " << count_t << endl;
	//	for (auto iTV : iT)
	//	{
	//		cout << iTV << endl;
	//	}
	//	++count_t;
	//}
}

void MeshInfo::getInitialBoundaryHalfedge(std::vector<Mesh::HalfedgeHandle>& _initial_halfedges) const
{
	size_t maxItr(mesh_->n_halfedges());

	std::vector<bool> is_checked(mesh_->n_halfedges(), false);

	// check all halfedges to find boundary
	for (auto he_itr = mesh_->halfedges_begin(); he_itr != mesh_->halfedges_end(); ++he_itr)
	{
		int current_halfedge_idx = he_itr.handle().idx();
		if (is_checked[current_halfedge_idx] == true)
			continue;

		is_checked[current_halfedge_idx] = true;

		if (false == mesh_->is_boundary(he_itr))
			continue;

		// boundry found                
		_initial_halfedges.push_back(*he_itr);
		size_t counter = 1;
		// traverse this boundary loop, tick them as checked
		auto next_he = mesh_->next_halfedge_handle(he_itr);
		while ((next_he != he_itr) && counter < maxItr)
		{
			assert(mesh_->is_boundary(next_he));
			assert(false == is_checked[next_he.idx()]);

			is_checked[next_he.idx()] = true;
			next_he = mesh_->next_halfedge_handle(next_he);
			counter++;
		}

		//std::cout << "MeshInfo::getInitialBoundaryHalfedge: Found cycle of length " << counter << std::endl;
		if (counter >= maxItr)
		{
			std::cout << "[WARNING] MeshInfo::getInitialBoundaryHalfedge:  Failed to close boundry loop." << std::endl;
			assert(false);
		}
	}
}
