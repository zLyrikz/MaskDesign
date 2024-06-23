#include "MeshPlaneIntersection.h"
#include "../Utility/VectorOperations.h"

#include<iostream>
using std::cout;
using std::endl;
using std::vector;
MeshPlaneIntersection::MeshPlaneIntersection()
{
	Init();
}

void MeshPlaneIntersection::setData(const Mesh* _mesh, const Mesh::Point& _normal, double _D)
{
	mesh_ = _mesh;
	normal_ = _normal;
	D_ = _D;
}

void MeshPlaneIntersection::setData(const Mesh* _mesh, const Mesh::Point& _normal, const Mesh::Point& _point_on_plane)
{
	double D = _point_on_plane.dot(_normal);
	setData(_mesh, _normal, D);
}

void MeshPlaneIntersection::FindIntersectFaces()
{
	intersect_elements_.clear();
	FindIntersectBoundaryFaces();

	FindIntersectFaces(intersect_boundary_face_, intersect_elements_);
}

void MeshPlaneIntersection::FindIntersectFaces_AllParts()
{
	intersect_elements_all_parts_.clear();
	FindIntersectBoundaryFaces();


	if (intersect_boundary_face_.empty())// a watertight mesh with no boundary, but for surface genus >=1, there still can be multiple parts
	{
		std::vector<IntersectFace> intersect_elements_all_parts;
		//FindIntersectFaces(intersect_boundary_face_, intersect_elements_all_parts);

		Mesh::FaceHandle intersect_face;
		while (FindOneIntersectFaceExceptSome(intersect_elements_all_parts, intersect_face))
		{
			// search an intersect part start from this intersect face
			std::vector<IntersectFace> intersect_part;
			FindIntersectFaces(std::vector<Mesh::FaceHandle>{intersect_face}, intersect_part);

			intersect_elements_all_parts_.push_back(intersect_part);

			VectorOperations append;
			append.AppendBack(intersect_elements_all_parts, intersect_part);
		}

	}
	else
	{
		std::vector<Mesh::FaceHandle> intersect_boundary_copy = intersect_boundary_face_;

		while (!intersect_boundary_copy.empty())
		{
			std::vector<IntersectFace> intersect_elements;
			FindIntersectFaces(intersect_boundary_copy, intersect_elements);
			intersect_elements_all_parts_.push_back(intersect_elements);

			// std::remove reference http://c.biancheng.net/view/6846.html
			// remove first and last faces in the intersect_elements, because they must be boundary and already been found
			auto end_iter = std::remove(intersect_boundary_copy.begin(), intersect_boundary_copy.end(), intersect_elements[0].face_);
			end_iter = std::remove(intersect_boundary_copy.begin(), end_iter, intersect_elements.back().face_);
			intersect_boundary_copy.erase(end_iter, intersect_boundary_copy.end());
		}
	}
}

bool MeshPlaneIntersection::FindLocalIntersect(const Mesh::FaceHandle& _start_face, const Mesh::FaceHandle& _end_face, 
	const Mesh::Point& _direction, std::vector<Mesh::Point>& _intersect_polyline) const
{
	//// get start face intersect information
	//std::vector<IntersectFace> element(1);
	//element[0].face_ = _start_face;
	//FindIntersectPolyline(element);
	//assert(element[0].halfedges_points_.size() == 2);

	//// find which side of the halfedge to go
	//Mesh::Point direction = element[0].halfedges_points_[1].second - element[0].halfedges_points_[0].second;
	//if (direction.dot(_direction) < 0)
	//{
	//	// guess a direction to search
	//	// change the halfedge stored order to search order
	//	std::swap(element[0].halfedges_points_[0], element[0].halfedges_points_[1]);
	//}
	
	// get start face intersect information
	std::vector<IntersectFace> element;
	PrepareStartFace(_start_face, _direction, element);

	// given a start face and a adjacent halfedge, search on until meet the end face
	std::vector<IntersectFace> all_elements = element;
	bool hit_end = FindLocalIntersect_WithStartInfo(_start_face, _end_face, all_elements);
	if (!hit_end)
	{
		// try another direction
		std::swap(element[0].halfedges_points_[0], element[0].halfedges_points_[1]);

		all_elements = element;// refresh the intersect elements
		hit_end = FindLocalIntersect_WithStartInfo(_start_face, _end_face, all_elements);
		if (!hit_end)
		{
			cout << "[WARINING MeshPlaneIntersection::FindLocalIntersect] find intersect fail! cannot hit the end face!" << endl;
			return false;
		}
	}

	getIntersectPolyline(all_elements, _intersect_polyline);

	return true;
}

bool MeshPlaneIntersection::FindLocalIntersect(const Mesh::FaceHandle& _start_face, const Mesh::Point& _start_point,  double _end_length, const Mesh::Point& _direction,
	std::vector<Mesh::Point>& _intersect_polyline, double& _actual_length) const
{
	// get start face intersect information
	std::vector<IntersectFace> element;
	PrepareStartFace(_start_face, _direction, element);

	//assert(_intersect_element.size() == 1);//have the start face info

	bool reach_length = true;// will be updated to false if cannot hit the given end face
	std::vector<IntersectFace> final_all_element;
	{
		final_all_element = element;
		// initialize the intersect information on the start face
		reach_length = FindLocalIntersect_WithStartInfo(_start_face, _start_point, _end_length, final_all_element, _actual_length);

		// check if the final polyline also fits the search direction
		Mesh::Point direction = final_all_element.back().halfedges_points_[1].second - final_all_element[0].halfedges_points_[0].second;
		bool correct_direction = direction.dot(_direction) > -1e-5;

		if (!reach_length || !correct_direction)
		{
			// try another direction
			std::swap(element[0].halfedges_points_[0], element[0].halfedges_points_[1]);

			// re-initialize these stuff
			std::vector<IntersectFace> all_element = element;
			double actual_length = 0.0;
			reach_length = FindLocalIntersect_WithStartInfo(_start_face, _start_point, _end_length, all_element, actual_length);

			// check if the final polyline also fits the search direction
			direction = all_element.back().halfedges_points_[1].second - all_element[0].halfedges_points_[0].second;
			correct_direction = direction.dot(_direction) > -1e-5;

			if (!reach_length || !correct_direction)
			{
				//direction still not coordinate 
				//cout << "[WARNING FROM MeshPlaneIntersection::FindLocalIntersect] search direction is hard to tell!" << endl;

			}
			else
			{
				// update the final state to this direction case
				final_all_element = std::move(all_element);
				_actual_length = actual_length;
			}
		}

	}

	getIntersectPolyline(final_all_element, _intersect_polyline);


	return reach_length;

}

void MeshPlaneIntersection::FindIntersectPolyline()
{
	FindIntersectFaces();

	//FindIntersectPolyline(intersect_elements_);
}

void MeshPlaneIntersection::FindIntersectPolyline_AllParts()
{
	FindIntersectFaces_AllParts();

	//for (int iPart = 0; iPart < intersect_elements_all_parts_.size(); ++iPart)
	//{
	//	FindIntersectPolyline(intersect_elements_all_parts_[iPart]);
	//}
}

void MeshPlaneIntersection::Init()
{
	D_ = 0.0;
}

void MeshPlaneIntersection::FindIntersectBoundaryFaces()
{
	intersect_boundary_face_.clear();

	//
	// find the intersection faces that's also on the mesh boundary
	//
	for (Mesh::FaceIter iterFace = mesh_->faces_begin(); iterFace != mesh_->faces_end(); ++iterFace)
	{
		if (mesh_->is_boundary(iterFace))
		{
			//// previous code, didn't check if it's the boundary edge on this face that has intersection with the plane
			//bool is_intersect = OneFaceIntersection(iterFace);
			//if (is_intersect == true)
			//{
			//	intersect_boundary_face_.push_back(iterFace);
			//}



			// check if the plane intersect with the boundary edge on this face
			std::vector<IntersectFace> this_face(1);
			this_face[0].face_ = *iterFace;
			FindIntersectPolyline(this_face);
			if (!this_face[0].halfedges_points_.empty())//do find intersection with some edges
			{
				assert(this_face[0].halfedges_points_.size() == 2);

				// there exists a intersect edge is a boundary; note remember we need to check the oppesite halfedge for boundary
				if (mesh_->is_boundary(this_face[0].halfedges_points_[0].first.opp()) ||
					mesh_->is_boundary(this_face[0].halfedges_points_[1].first.opp()))
				{
					intersect_boundary_face_.push_back(iterFace);
				}
			}
			



		}
	}
}

bool MeshPlaneIntersection::OneFaceIntersection(Mesh::FaceHandle _face) const
{
	// check if this face has intersection with the mesh
	// It is not intersected equivalent to sign(face_vertex.dot(plane_normal) - D) are the same for all vertices of this face
	
	bool is_intersect = false;

	Mesh::FaceVertexIter vertex = mesh_->cfv_iter(_face);
	Mesh::Point point = mesh_->point(vertex);
	bool sign0 = ((point.dot(normal_) - D_) < 0);// sign of the first vertex

	for (; vertex.is_valid(); ++vertex) // remember: We'd better write code that's easy to understand(straightforward).
	{
		point = mesh_->point(vertex);
		bool sign = ((point.dot(normal_) - D_) < 0);// sign of the current vertex
		if (sign != sign0)
		{
			is_intersect = true;
			break;
		}
	}

	return is_intersect;
}

void MeshPlaneIntersection::FindIntersectFaces(const std::vector<Mesh::FaceHandle>& _intersect_boundary_face,
	std::vector<IntersectFace>& _intersect_elements) const
{
	vector<bool> is_stored(mesh_->n_faces(), false);// indicate whether a face of the mesh is stored in intersect_face_ or not
	//cout << "mesh_->n_faces() = " << mesh_->n_faces() << endl;
	{
		// If 0 < number of intersected faces on the boundary , then the search starts from one of it
		// else the search starts from one of the intersected face
		Mesh::FaceHandle face;
		if (_intersect_boundary_face.size() != 0)
		{
			face = _intersect_boundary_face[0];
		}
		else
		{
			// find one intersected face
			for (auto iterFace = mesh_->faces_begin(); iterFace != mesh_->faces_end(); ++iterFace)
			{
				Mesh::FaceHandle current_face = *iterFace;
				bool is_intersect = OneFaceIntersection(current_face);
				if (is_intersect)
				{
					face = current_face;
					break;
				}
			}
		}

		if (face.is_valid())
		{
			bool is_end = false; // indicate whether the search should end
			do
			{
				is_end = true; // Later, change is_end to be false, if we find a new intersected face. 
				// so the while loop will stop if we can't find a new intersected face

				//cout << "face.idx() = " << face.idx() << endl;

				// store current face if it's not been stored
				if (!is_stored[face.idx()])
				{
					std::vector<IntersectFace> element(1);
					element[0].face_ = face;
					FindIntersectPolyline(element);

					_intersect_elements.push_back(element[0]);


					is_stored[face.idx()] = true;
				}

				// traverse the edge ajacent faces of the current face
				// break the traverse if we find an ajacent face that has intersection with the plane
				// and then go to the next while loop
				for (int iEdge = 0; iEdge < 2; ++iEdge)
				{
					Mesh::FaceHandle ajacent_face = _intersect_elements.back().halfedges_points_[iEdge].first.opp().face();
					if (ajacent_face.is_valid())
					{
						bool is_intersect = OneFaceIntersection(ajacent_face);
						//cout << "ajacent_face.idx() = " << ajacent_face.idx() << endl;

						if (is_intersect && !is_stored[ajacent_face.idx()])
						{
							face = ajacent_face; // update a new face
							is_end = false;
							break;
						}
					}
				}
				//for (Mesh::FaceHalfedgeIter iterHalfedge = mesh_->cfh_iter(face); iterHalfedge.is_valid(); ++iterHalfedge)
				//{
				//	Mesh::FaceHandle ajacent_face = (*iterHalfedge).opp().face();
				//	if (ajacent_face.is_valid())
				//	{
				//		bool is_intersect = OneFaceIntersection(ajacent_face);
				//		//cout << "ajacent_face.idx() = " << ajacent_face.idx() << endl;

				//		if (is_intersect && !is_stored[ajacent_face.idx()])
				//		{
				//			face = ajacent_face; // update a new face
				//			is_end = false;
				//			break;
				//		}
				//	}
				//}
			} while (!is_end);
		}
	}
}

void MeshPlaneIntersection::FindIntersectPolyline(std::vector<IntersectFace>& _intersect_elements) const
{
	// this method is used by FindIntersectFace

	for (int iFace = 0; iFace < _intersect_elements.size(); ++iFace)
	{
		// find intersect halfedges
		_intersect_elements[iFace].halfedges_points_.reserve(2);

		for (Mesh::FaceHalfedgeIter iterHalfedge = mesh_->cfh_iter(_intersect_elements[iFace].face_);
			iterHalfedge.is_valid(); ++iterHalfedge)
		{
			OpenMesh::SmartHalfedgeHandle halfedge = *iterHalfedge;
			//Mesh::Point from_vertex = mesh_->point(halfedge.from());
			//Mesh::Point to_vertex = mesh_->point(halfedge.to());
			//bool sign_from = ((from_vertex.dot(normal_) - D_) < 0);
			//bool sign_to = ((to_vertex.dot(normal_) - D_) < 0);

			//if (sign_from != sign_to) // this halfedge intersects with the plane
			//{
			//	// find the intersect point
			//	Mesh::Point line = to_vertex - from_vertex;
			//	double translation = (D_ - from_vertex.dot(normal_)) / line.dot(normal_);
			//	Mesh::Point point_on_plane = from_vertex + translation * line;

			//	std::pair<OpenMesh::SmartHalfedgeHandle, Mesh::Point> intersects_pair(halfedge, point_on_plane);
			//	_intersect_elements[iFace].halfedges_points_.push_back(intersects_pair);
			//}

			if (HalfedgeIntersectPlane(halfedge, normal_, D_))
			{
				Mesh::Point from_vertex = mesh_->point(halfedge.from());
				Mesh::Point to_vertex = mesh_->point(halfedge.to());
				Mesh::Point point_on_plane;
				LineIntersectPlane(from_vertex, to_vertex, normal_, D_, point_on_plane);

				std::pair<OpenMesh::SmartHalfedgeHandle, Mesh::Point> intersects_pair(halfedge, point_on_plane);
				_intersect_elements[iFace].halfedges_points_.push_back(intersects_pair);
			}
		}
	}
}

void MeshPlaneIntersection::getIntersectPoints(const std::vector<IntersectFace>& _intersect_elements, std::vector<Mesh::Point>& _points) const
{
	int face_num = _intersect_elements.size();
	_points = std::vector<Mesh::Point>(face_num);
	for (int i = 0; i < face_num; ++i)
	{
		if (_intersect_elements[i].halfedges_points_.size() != 2)
		{
			cout << "WARNING unexpected! intersect edge number on a face != 2" << endl;
		}

		_points[i] = 0.5 * (_intersect_elements[i].halfedges_points_[0].second +
			_intersect_elements[i].halfedges_points_[1].second);
	}
}

bool MeshPlaneIntersection::FindOneIntersectFaceExceptSome(const std::vector<IntersectFace>& _intersect_elements, Mesh::FaceHandle& _face) const
{
	// get the face handles from the elements
	std::vector<Mesh::FaceHandle> intersect_faces;
	intersect_faces.reserve(_intersect_elements.size());
	for (auto& iFace : _intersect_elements)
	{
		intersect_faces.push_back(iFace.face_);
	}

	// find one intersected face except the given intersect_elements
	for (auto iterFace = mesh_->faces_begin(); iterFace != mesh_->faces_end(); ++iterFace)
	{
		Mesh::FaceHandle current_face = *iterFace;
		bool is_intersect = OneFaceIntersection(current_face);
		if (is_intersect)
		{
			auto iter_find = std::find(intersect_faces.begin(), intersect_faces.end(), current_face);
			if (iter_find == intersect_faces.end())// find a new one
			{
				_face = current_face;
				return true;
			}
		}
	}

	return false;
}

bool MeshPlaneIntersection::HalfedgeIntersectPlane(const OpenMesh::SmartHalfedgeHandle& _halfedge, const Mesh::Point& _normal, double _D) const
{
	Mesh::Point from_vertex = mesh_->point(_halfedge.from());
	Mesh::Point to_vertex = mesh_->point(_halfedge.to());
	bool sign_from = ((from_vertex.dot(_normal) - _D) < 0);
	bool sign_to = ((to_vertex.dot(_normal) - _D) < 0);


	return (sign_to != sign_from);
}

void MeshPlaneIntersection::LineIntersectPlane(const Mesh::Point& _start, const Mesh::Point& _end, const Mesh::Point& _normal, double _D, Mesh::Point& _intersect) const
{
	// find the intersect point
	Mesh::Point line = _end - _start;
	double translation = (_D - _start.dot(_normal)) / line.dot(_normal);
	_intersect = _start + translation * line;
}

bool MeshPlaneIntersection::FindLocalIntersect_WithStartInfo(const Mesh::FaceHandle& _start_face, const Mesh::FaceHandle& _end_face,
	std::vector<IntersectFace>& _intersect_element) const
{
	assert(_intersect_element.size() == 1);//have the start face info

	bool hit_end = true;// will be updated to false if cannot hit the given end face
	{
		Mesh::FaceHandle current_face = _start_face;
		OpenMesh::SmartHalfedgeHandle current_halfedge = _intersect_element.back().halfedges_points_[1].first;

		while (current_face != _end_face && hit_end)
		{

			bool next_face_is_valid = SearchNextFace(current_face, current_halfedge, _intersect_element);
			if (!next_face_is_valid || current_face == _start_face)
			{
				hit_end = false;// still not hit end face even when no more valid adjacent face exists
				// this cause while loop ends
			}


			//Mesh::Point adjacent_point = _intersect_element.back().halfedges_points_[1].second;
			//OpenMesh::SmartHalfedgeHandle adjacent_halfedge = current_halfedge.opp();
			//
			//// update the next intersect face info
			//current_face = adjacent_halfedge.face();
			//if (current_face.is_valid())
			//{
			//	_intersect_element.push_back(IntersectFace());
			//	_intersect_element.back().face_ = current_face;
			//	_intersect_element.back().halfedges_points_.push_back(std::make_pair(adjacent_halfedge, adjacent_point));

			//	// find the other intersect half edge, update current halfedge to this one
			//	for (Mesh::FaceHalfedgeIter iterHalfedge = mesh_->cfh_iter(current_face); iterHalfedge.is_valid(); ++iterHalfedge)
			//	{
			//		OpenMesh::SmartHalfedgeHandle halfedge = *iterHalfedge;

			//		if (halfedge != adjacent_halfedge && HalfedgeIntersectPlane(halfedge, normal_, D_))
			//		{
			//			// find the intersect point on the halfedge
			//			Mesh::Point from_vertex = mesh_->point(halfedge.from());
			//			Mesh::Point to_vertex = mesh_->point(halfedge.to());
			//			Mesh::Point point_on_plane;
			//			LineIntersectPlane(from_vertex, to_vertex, normal_, D_, point_on_plane);

			//			// store the intersect information
			//			_intersect_element.back().halfedges_points_.push_back(std::make_pair(halfedge, point_on_plane));

			//			// update halfedge for next search
			//			current_halfedge = std::move(halfedge);
			//		}
			//	}
			//	assert(current_halfedge != adjacent_halfedge.opp());// a new halfedge is sure to be updated
			//}
			//else
			//{
			//	assert(mesh_->is_boundary(adjacent_halfedge));
			//	hit_end = false;// still not hit end face even when no more valid adjacent face exists
			//	// this cause while loop ends
			//}
		}
	}

	return hit_end;
}

bool MeshPlaneIntersection::FindLocalIntersect_WithStartInfo(const Mesh::FaceHandle& _start_face, const Mesh::Point& _start_point, double _end_length,
	std::vector<IntersectFace>& _intersect_element, double& _actual_length) const
{
	assert(_intersect_element.size() == 1);//have the start face info

	Mesh::FaceHandle current_face = _start_face;
	OpenMesh::SmartHalfedgeHandle current_halfedge = _intersect_element.back().halfedges_points_[1].first;
	_actual_length = (_start_point - _intersect_element.back().halfedges_points_[1].second).norm();
	bool reach_length = true;
	while (_actual_length < _end_length && reach_length)
	{
		bool next_face_is_valid = SearchNextFace(current_face, current_halfedge, _intersect_element);
		if (!next_face_is_valid || current_face == _start_face)
		{
			reach_length = false;// still not reach certain length even when no more valid adjacent face exists
			// this cause while loop ends
		}
		else
		{
			_actual_length += ComputeIntersectLineSegmentLength(_intersect_element.back());// add the new intersect line segment to the polyline length
		}
	}

	return reach_length;
}

void MeshPlaneIntersection::PrepareStartFace(const Mesh::FaceHandle& _start_face, const Mesh::Point& _direction, std::vector<IntersectFace>& _intersect_element) const
{
	// get start face intersect information
	_intersect_element.resize(1);
	_intersect_element[0].face_ = _start_face;
	FindIntersectPolyline(_intersect_element);
	assert(_intersect_element[0].halfedges_points_.size() == 2);

	// find which side of the halfedge to go
	Mesh::Point direction = _intersect_element[0].halfedges_points_[1].second - _intersect_element[0].halfedges_points_[0].second;
	if (direction.dot(_direction) < 0)
	{
		// guess a direction to search
		// change the halfedge stored order to search order
		std::swap(_intersect_element[0].halfedges_points_[0], _intersect_element[0].halfedges_points_[1]);
	}
}

bool MeshPlaneIntersection::SearchNextFace(Mesh::FaceHandle& _current_face, OpenMesh::SmartHalfedgeHandle& _current_halfedge, std::vector<IntersectFace>& _intersect_element) const
{
	const Mesh::Point& adjacent_point = _intersect_element.back().halfedges_points_[1].second;
	const OpenMesh::SmartHalfedgeHandle& adjacent_halfedge = _current_halfedge.opp();

	// update the next intersect face info
	if (adjacent_halfedge.face().is_valid())
	{
		_current_face = adjacent_halfedge.face();

		_intersect_element.push_back(IntersectFace());
		_intersect_element.back().face_ = _current_face;
		_intersect_element.back().halfedges_points_.push_back(std::make_pair(adjacent_halfedge, adjacent_point));

		// find the other intersect half edge, update current halfedge to this one
		for (Mesh::FaceHalfedgeIter iterHalfedge = mesh_->cfh_iter(_current_face); iterHalfedge.is_valid(); ++iterHalfedge)
		{
			OpenMesh::SmartHalfedgeHandle halfedge = *iterHalfedge;

			if (halfedge != adjacent_halfedge && HalfedgeIntersectPlane(halfedge, normal_, D_))
			{
				// find the intersect point on the halfedge
				Mesh::Point from_vertex = mesh_->point(halfedge.from());
				Mesh::Point to_vertex = mesh_->point(halfedge.to());
				Mesh::Point point_on_plane;
				LineIntersectPlane(from_vertex, to_vertex, normal_, D_, point_on_plane);

				// store the intersect information
				_intersect_element.back().halfedges_points_.push_back(std::make_pair(halfedge, point_on_plane));

				// update halfedge for next search
				_current_halfedge = std::move(halfedge);
				break;// should be no more other halfedge that has intersection!// try add this to an assertion			
			}
		}
		assert(_current_halfedge != adjacent_halfedge.opp());// a new halfedge is sure to be updated
	}
	else
	{
		assert(mesh_->is_boundary(adjacent_halfedge));
		return false;
	}


	return true;
}

double MeshPlaneIntersection::ComputeIntersectLineSegmentLength(const IntersectFace& _face) const
{
	assert(_face.halfedges_points_.size() == 2);
	return (_face.halfedges_points_[0].second - _face.halfedges_points_[1].second).norm();
}

void MeshPlaneIntersection::getIntersectPolyline(const std::vector<IntersectFace>& _intersect_elements, std::vector<Mesh::Point>& _polyline) const
{
	// get the intersect polyline from the ordered intersect elements
	_polyline.clear();
	_polyline.reserve(_intersect_elements.size() + 1);// number of intersect edges = number of intersec faces +1
	_polyline.push_back(_intersect_elements[0].halfedges_points_[0].second);
	for (auto& iIntersect : _intersect_elements)
	{
		_polyline.push_back(iIntersect.halfedges_points_[1].second);
	}
}



void MeshPlaneIntersection::getIntersectBoundaryFaces(std::vector<Mesh::FaceHandle>& _intersect_boundary_face) const
{
	_intersect_boundary_face = intersect_boundary_face_;
}

void MeshPlaneIntersection::getIntersectFaces(std::vector<Mesh::FaceHandle>& _intersect_face) const
{
	int face_num = intersect_elements_.size();
	_intersect_face = std::vector<Mesh::FaceHandle>(face_num);
	for (int i = 0; i < face_num; ++i)
	{
		_intersect_face[i] = intersect_elements_[i].face_;
	}
}

void MeshPlaneIntersection::getIntersectPoints(std::vector<Mesh::Point>& _points) const
{
	getIntersectPoints(intersect_elements_, _points);
}

void MeshPlaneIntersection::getIntersectPoints_ClosetToAPoint(const Mesh::Point& _point, std::vector<Mesh::Point>& _points) const
{
	// get all parts of the polyline
	int num_parts = intersect_elements_all_parts_.size();
	std::vector<std::vector<Mesh::Point>> points_all_parts(num_parts);
	//for (int iPart = 0; iPart < num_parts; ++iPart)
	//{
	//	getIntersectPoints(intersect_elements_all_parts_[iPart], points_all_parts[iPart]);
	//}
	getIntersectPoints_AllParts(points_all_parts);

	// find a point in all polylines closest to the given point
	int closest_part_idx = 0;// to be found

	if (num_parts > 1)
	{
		double min_distance = FLT_MAX;
		for (int iPart = 0; iPart < num_parts; ++iPart)
		{
			for (int iPoint = 0; iPoint < points_all_parts[iPart].size(); ++iPoint)
			{
				double distance = (points_all_parts[iPart][iPoint] - _point).norm();
				if (distance < min_distance)
				{
					min_distance = distance;
					closest_part_idx = iPart;// update best part
				}
			}
		}
	}
	if (points_all_parts.size() > 0)
	{
		_points = points_all_parts[closest_part_idx];
	}
}

void MeshPlaneIntersection::getIntersectPoints_AllParts(std::vector<std::vector<Mesh::Point>>& _points) const
{
	// get all parts of the polyline
	int num_parts = intersect_elements_all_parts_.size();
	_points.resize(num_parts);
	for (int iPart = 0; iPart < num_parts; ++iPart)
	{
		getIntersectPoints(intersect_elements_all_parts_[iPart], _points[iPart]);
	}
}
