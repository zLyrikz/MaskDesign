#include "ConnectFrontNBack.h"
#include "../MeshProcessing/MeshInfo.h"
#include "../MeshProcessing/MeshOffset.h"
#include "../MeshProcessing/AabbTree.h"
#include "../MeshProcessing/SubMesh.h"
#include "../MeshProcessing/FindNearestPoint.h"
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <iostream>

using std::cout;
using std::endl;
ConnectFrontNBack::ConnectFrontNBack()
{
}

ConnectFrontNBack::~ConnectFrontNBack()
{
}

void ConnectFrontNBack::CreateIntermediatePart(const std::vector<Mesh::Point>& _back_inside_boundary,
	const std::vector<Mesh::Point>& _back_outside_boundary, 
	const Mesh& _cushion_surface, double _distance, 
	Mesh& _intermediate, double _csg_tolerance) const
{
	//find the connecting boundaries
	std::vector<Mesh::Point> cushion_connecting_boundary;
	std::vector<Mesh::Point> offset_connecting_boundary;
	GetCushionConnectingBoundary(_cushion_surface, cushion_connecting_boundary, offset_connecting_boundary, _distance);

	// create intermediate part from the boundary
	CreateIntermediatePart(_back_inside_boundary, _back_outside_boundary, offset_connecting_boundary, cushion_connecting_boundary, _intermediate, _csg_tolerance);
}

void ConnectFrontNBack::CreateIntermediatePart(
	const std::vector<Mesh::Point>& _back_inside_boundary, const std::vector<Mesh::Point>& _back_outside_boundary, 
	const std::vector<Mesh::Point>& _cushion_inside_boundary, const std::vector<Mesh::Point>& _cushion_outside_boundary, 
	Mesh& _intermediate, double _csg_tolerance) const
{

	std::vector<Mesh::Point> cushion_outside_boundary = _cushion_outside_boundary;
	std::vector<Mesh::Point> cushion_inside_boundary = _cushion_inside_boundary;
	std::vector<Mesh::Point> correspond_back_outside;//resorting of back boundary points so that correspond to the cushion boundary
	std::vector<Mesh::Point> correspond_back_inside;


	GetBoundaryCorrespondenceAndAddTolerance(
		_back_inside_boundary, _back_outside_boundary,
		cushion_inside_boundary, cushion_outside_boundary,
		correspond_back_outside, correspond_back_inside, _csg_tolerance);

	SubMesh construct_intermediate;
	construct_intermediate.CuboidTubeMeshFrom4Polylines_Period(
		cushion_outside_boundary, cushion_inside_boundary, correspond_back_inside, correspond_back_outside, _intermediate);

}

void ConnectFrontNBack::ComputeIntermediatePartPolylineLayers(
	const std::vector<Mesh::Point>& _back_inside_boundary, const std::vector<Mesh::Point>& _back_outside_boundary, 
	const Mesh& _cushion_surface, double _distance, 
	std::vector<std::vector<Mesh::Point>>& _polyline_layers, int _num_layer, bool _include_boundary) const
{
	//find the connecting boundaries
	std::vector<Mesh::Point> cushion_connecting_boundary;
	std::vector<Mesh::Point> offset_connecting_boundary;
	GetCushionConnectingBoundary(_cushion_surface, cushion_connecting_boundary, offset_connecting_boundary, _distance);

	ComputeIntermediatePartPolylineLayers(_back_inside_boundary, _back_outside_boundary, 
		offset_connecting_boundary, cushion_connecting_boundary, _polyline_layers, _num_layer, _include_boundary);
}

void ConnectFrontNBack::ComputeIntermediatePartPolylineLayers(
	const std::vector<Mesh::Point>& _back_inside_boundary, const std::vector<Mesh::Point>& _back_outside_boundary, 
	const std::vector<Mesh::Point>& _cushion_inside_boundary, const std::vector<Mesh::Point>& _cushion_outside_boundary, 
	std::vector<std::vector<Mesh::Point>>& _polyline_layers, int _num_layer, bool _include_boundary) const
{
	std::vector<Mesh::Point> cushion_outside_boundary = _cushion_outside_boundary;
	std::vector<Mesh::Point> cushion_inside_boundary = _cushion_inside_boundary;
	std::vector<Mesh::Point> correspond_back_outside;//resorting of back boundary points so that correspond to the cushion boundary
	std::vector<Mesh::Point> correspond_back_inside;

	// don't add tolerance here, so no change on the cushion boundary
	GetBoundaryCorrespondenceAndAddTolerance(
		_back_inside_boundary, _back_outside_boundary,
		cushion_inside_boundary, cushion_outside_boundary,
		correspond_back_outside, correspond_back_inside, 0.0);

	assert(cushion_inside_boundary.size() == cushion_outside_boundary.size());
	assert(correspond_back_outside.size() == cushion_outside_boundary.size());
	assert(correspond_back_outside.size() == correspond_back_inside.size());
	int num_points = cushion_inside_boundary.size();
	//_polyline_layers.resize(_num_layer, std::vector<Mesh::Point>(2 * num_points));
	_polyline_layers.resize(_num_layer, std::vector<Mesh::Point>(num_points));
	for (int iLayer = 0; iLayer < _num_layer; ++iLayer)
	{
		// we have _num_layer+1 intervals
		double interval_length = 1.0 / (_num_layer + 1);
		double t = (iLayer + 1) * interval_length;
		if (_include_boundary && _num_layer >= 2)
		{
			interval_length = 1.0 / (_num_layer - 1);
			t = iLayer * interval_length;
		}

		// outside layer
		for (int iPoint = 0; iPoint < num_points; ++iPoint)
		{
			_polyline_layers[iLayer][iPoint] = cushion_outside_boundary[iPoint] + t * (correspond_back_outside[iPoint] - cushion_outside_boundary[iPoint]);
		}

		//// inside layer
		//for (int iPoint = 0; iPoint < num_points; ++iPoint)
		//{
		//	_polyline_layers[iLayer][iPoint + num_points] = cushion_inside_boundary[iPoint] + t * (correspond_back_inside[iPoint] - cushion_inside_boundary[iPoint]);
		//}
	}

}

void ConnectFrontNBack::CreateOverlap(const Mesh& _back, Mesh& _cushion_surface) const
{
	Mesh::HalfedgeHandle connecting_boundary;
	FindConnectingBoundary(_cushion_surface, connecting_boundary);

	AabbTree back_tree;
	back_tree.constructTree(&_back);

	{
		Mesh::HalfedgeHandle current_boundary = connecting_boundary;
		do
		{
			Mesh::VertexHandle current_vertexhandle = _cushion_surface.from_vertex_handle(current_boundary);
			Mesh::Point current_point = _cushion_surface.point(current_vertexhandle);
			Mesh::Point negative_z = Mesh::Point(0.0, 0.0, -1.0);
			Mesh::Point not_used_point;
			bool is_intersect = back_tree.RayMeshIntersection(current_point, negative_z, not_used_point);

			if (is_intersect == false)
			{
				// translate along positive z
				Mesh::Point intersect_point;
				bool is_intersect2 = back_tree.RayMeshIntersection(current_point, -negative_z, intersect_point);
				if (is_intersect2 == true)
				{
					Mesh::Point tolerance = - 0.1 * negative_z;
					_cushion_surface.set_point(current_vertexhandle, intersect_point + tolerance);
				}
				else
				{
					cout << "[WARNING] from ConnectFrontNBack::CreateOverlap, cannot translate to create overlap" << endl;
				}
			}
			// update for next loop
			current_boundary = _cushion_surface.next_halfedge_handle(current_boundary);
		} while (current_boundary != connecting_boundary);
	}

}

void ConnectFrontNBack::CreateWatertightOverlap(const Mesh& _back, const Mesh& _cushion_surface, double _distance, Mesh& _watertight) const
{
	Mesh original = _cushion_surface;
	CreateOverlap(_back, original);

	Mesh offset;
	MeshOffset offset_it;
	offset_it.OffsetSurface(original, offset, _distance);
	offset_it.FlipNormal(offset);
	CreateOverlap(_back, offset);


	MeshOffset get_watertight;
	get_watertight.WatertightFromOffset(original, offset, _watertight);
}

void ConnectFrontNBack::CsgConnect(const Mesh& _back, const Mesh& _cushion_watertight, Mesh& _csg_connected) const
{
	// convert mesh type to libigl

	Eigen::MatrixXd back_vertex;
	Eigen::MatrixXi back_face;
	Eigen::MatrixXd cushion_watertight_vertex;
	Eigen::MatrixXi cushion_watertight_face;
	SubMesh convert_mesh_type;
	convert_mesh_type.Mesh2EigenMatrix(_back, back_vertex, back_face);
	convert_mesh_type.Mesh2EigenMatrix(_cushion_watertight, cushion_watertight_vertex, cushion_watertight_face);

	// union
	Eigen::MatrixXd csg_connected_vertex;
	Eigen::MatrixXi csg_connected_face;
	igl::MeshBooleanType boolean_type(igl::MESH_BOOLEAN_TYPE_UNION);

	igl::copyleft::cgal::mesh_boolean(back_vertex, back_face, cushion_watertight_vertex, cushion_watertight_face, 
		boolean_type, csg_connected_vertex, csg_connected_face);

	// convert mesh type to openmesh
	convert_mesh_type.MeshFromEigenMatrix(_csg_connected, csg_connected_vertex, csg_connected_face);

}

void ConnectFrontNBack::CsgIntersect(const Mesh& _mesh1, const Mesh& _mesh2, Mesh& _csg) const
{
	// convert mesh type to libigl

	Eigen::MatrixXd mesh1_vertex;
	Eigen::MatrixXi mesh1_face;
	Eigen::MatrixXd mesh2_vertex;
	Eigen::MatrixXi mesh2_face;
	SubMesh convert_mesh_type;
	convert_mesh_type.Mesh2EigenMatrix(_mesh1, mesh1_vertex, mesh1_face);
	convert_mesh_type.Mesh2EigenMatrix(_mesh2, mesh2_vertex, mesh2_face);

	// union
	Eigen::MatrixXd csg_vertex;
	Eigen::MatrixXi csg_face;
	igl::MeshBooleanType boolean_type(igl::MESH_BOOLEAN_TYPE_INTERSECT);

	igl::copyleft::cgal::mesh_boolean(mesh1_vertex, mesh1_face, mesh2_vertex, mesh2_face,
		boolean_type, csg_vertex, csg_face);

	// convert mesh type to openmesh
	convert_mesh_type.MeshFromEigenMatrix(_csg, csg_vertex, csg_face);
}

void ConnectFrontNBack::GetBoundaryCorrespondenceInXOY(const std::vector<Mesh::Point>& _source_curve,
	const std::vector<Mesh::Point>& _target_curve, std::vector<Mesh::Point>& _correspond_source) const
{
	// find center point of target curve
	OpenMesh::Vec2d target_center(0, 0);
	for (const auto& iPoint : _target_curve)
	{
		target_center[0] += iPoint[0];
		target_center[1] += iPoint[1];
	}
	target_center /= double(_target_curve.size());
	//cout << "target_center=" << target_center << endl;

	// find the polar angle of the two curves, with origin at target_center
	std::vector<double> source_angle;
	GetPointsPolarAngleInXOY(source_angle, _source_curve, target_center);
	std::vector<double> target_angle;
	GetPointsPolarAngleInXOY(target_angle, _target_curve, target_center);

	// sort the source angles from small to big, instead of store the ranked value, we store the index
	// i.e. index_vector[i] gives the index of i-th biggest value in source_angles
	// check https://www.cnblogs.com/peimingzhang/p/13261843.html
	int num_source = source_angle.size();
	std::vector<int> source_rank_index(num_source);
	std::iota(source_rank_index.begin(), source_rank_index.end(), 0);
	std::sort(source_rank_index.begin(), source_rank_index.end(),
		[&source_angle](int a, int b) { return source_angle[a] < source_angle[b]; });

	_correspond_source.reserve(_target_curve.size());
	for (const auto& iTarget : target_angle)
	{
		// find the closest two angles in source from this angle in target


		int rank = 0;// how many angles in source is smaller than this target angle
		for (int iSource = 0; iSource < num_source; ++iSource)
		{
			//cout << "source_angle[source_rank_index[iSource]]=" << source_angle[source_rank_index[iSource]] << " iTarget=" << iTarget << endl;

			if (source_angle[source_rank_index[iSource]] < iTarget)
			{
				++rank;
			}
			else
			{
				// rest source angles surely also larget than this target angle
				break;
			}
		}
		int smaller_id = 0;
		int larger_id = 0;
		// note we need to handle period boundary case
		if (rank == 0 || rank == num_source)
		{
			smaller_id = source_rank_index[num_source - 1];
			larger_id = source_rank_index[0];
		}
		else
		{
			smaller_id = source_rank_index[rank - 1];
			larger_id = source_rank_index[rank];
		}
		//cout << "rank=" << rank << " smaller_id=" << smaller_id << " larger_id=" << larger_id << endl;


		// we take the corresponding point as
		// convex combination of the two source point
		// the conbination weight is approximated by the angle ratio
		// so nex we find the angle ratio
		double source1 = source_angle[smaller_id];
		double source2 = source_angle[larger_id];
		if (rank == 0) // handle period case
		{
			source1 -= 2 * M_PI;
		}
		else if (rank == num_source)
		{
			source2 += 2 * M_PI;
		}
		double angle_source1_to_2 = source2 - source1;
		double angle_source1_to_target = iTarget - source1;
		double angle_ratio = angle_source1_to_target / angle_source1_to_2;

		//Mesh::Point correspond_point;
		//OpenMesh::Vec2d source_vector1(_source_curve[smaller_id][0] - target_center[0], _source_curve[smaller_id][1] - target_center[1]);
		//OpenMesh::Vec2d source_vector2(_source_curve[larger_id][0] - target_center[0], _source_curve[larger_id][1] - target_center[1]);
		//double length = ((source_vector1 + source_vector2) / 2.0).norm();
		//correspond_point[0] = length * cos(iTarget) + target_center[0];
		//correspond_point[1] = length * sin(iTarget) + target_center[1];
		//correspond_point[2] = (_source_curve[smaller_id][2] + _source_curve[larger_id][2]) / 2.0;
		////cout << correspond_point << endl;
		//cout << "source_vector1=" << source_vector1 << endl;
		//cout << "source_vector2=" << source_vector2 << endl;
		//cout <<"length=" << length << " iTarget=" << iTarget << " cos=" << cos(iTarget) << " sin=" << sin(iTarget) << endl;

		Mesh::Point correspond_point = _source_curve[smaller_id] * (1.0 - angle_ratio) + angle_ratio * _source_curve[larger_id];
		_correspond_source.push_back(std::move(correspond_point));
	}

}

void ConnectFrontNBack::GetBoundaryCorrespondenceAndAddTolerance(
	const std::vector<Mesh::Point>& _back_inside_boundary, const std::vector<Mesh::Point>& _back_outside_boundary, 
	std::vector<Mesh::Point>& _cushion_inside_boundary, std::vector<Mesh::Point>& _cushion_outside_boundary,
	std::vector<Mesh::Point>& _correspond_back_outside, std::vector<Mesh::Point>& _correspond_back_inside,
	double _csg_tolerance) const
{
	Mesh::Point tolerance(0.0, 0.0, _csg_tolerance);
	// add a little tolerance for CSG
	assert(_cushion_outside_boundary.size() == _cushion_inside_boundary.size());
	int num_point = _cushion_outside_boundary.size();
	for (int iPoint = 0; iPoint < num_point; ++iPoint)
	{
		_cushion_outside_boundary[iPoint]  -= /*0.1 * */tolerance;
		_cushion_inside_boundary[iPoint] -= /*0.1 * */tolerance;
	}

	// find reasonable boundary correspondence (also add tolerance here for csg)
	// here we use closest point (ANN)
	_correspond_back_outside.clear();//resorting of mesh1 boundary points so that correspond to the cushion boundary
	_correspond_back_outside.reserve(_cushion_outside_boundary.size());
	FindNearestPoint find_from_outside;
	find_from_outside.setPoints(_back_outside_boundary);

	for (const Mesh::Point& iVertex : _cushion_outside_boundary)
	{
		int ann_point_id = find_from_outside.AnnInPoints(iVertex);
		_correspond_back_outside.push_back(_back_outside_boundary[ann_point_id] + tolerance);
	}

	// same process for the inside correspondence
	_correspond_back_inside.clear();//resorting of back boundary points so that correspond to the cushion boundary
	_correspond_back_inside.reserve(_cushion_inside_boundary.size());
	FindNearestPoint find_from_inside;
	find_from_inside.setPoints(_back_inside_boundary);
	for (Mesh::Point& iVertex : _cushion_inside_boundary)
	{
		int ann_point_id = find_from_inside.AnnInPoints(iVertex);
		_correspond_back_inside.push_back(_back_inside_boundary[ann_point_id] + tolerance);
	}
}

void ConnectFrontNBack::GetCushionConnectingBoundary(const Mesh& _cushion_surface, 
	std::vector<Mesh::Point>& _cushion_boundary, std::vector<Mesh::Point>& _offset_boundary, double _distance) const
{
	//find the connecting boundaries
	// 
	// offset cushion surface
	Mesh offset;
	MeshOffset offset_it;
	offset_it.OffsetSurface(_cushion_surface, offset, _distance);
	//offset_it.FlipNormal(offset);

	// get connecting boundary
	std::vector<Mesh::VertexHandle> cushion_connecting_boundary_handle;
	std::vector<Mesh::VertexHandle> offset_connecting_boundary_handle;
	FindConnectingBoundary(_cushion_surface, cushion_connecting_boundary_handle);
	FindConnectingBoundary(offset, offset_connecting_boundary_handle);
	assert(cushion_connecting_boundary_handle.size() == offset_connecting_boundary_handle.size());

	int num_point = cushion_connecting_boundary_handle.size();
	_cushion_boundary.resize(num_point);
	_offset_boundary.resize(num_point);
	for (int iPoint = 0; iPoint < num_point; ++iPoint)
	{
		_cushion_boundary[iPoint] = _cushion_surface.point(cushion_connecting_boundary_handle[iPoint]);
		_offset_boundary[iPoint] = offset.point(offset_connecting_boundary_handle[iPoint]);
	}
	// check if boundary orientation is same
	if (_cushion_boundary.size() > 1)
	{
		Mesh::Point direction1 = _cushion_boundary[1] - _cushion_boundary[0];
		Mesh::Point direction2 = _offset_boundary[1] - _offset_boundary[0];
		if (direction1.dot(direction2) < 0.0)// orientation
		{
			cout << "[WARNING FROM ConnectFrontNBack::CreateIntermediatePart] boundary incongruent orientation, will reorient them" << endl;
			// reverse offset boundary
			reverse(_offset_boundary.begin(), _offset_boundary.end());
		}
	}
}

void ConnectFrontNBack::GetPointsPolarAngleInXOY(std::vector<double>& _angles, const std::vector<Mesh::Point>& _points, const OpenMesh::Vec2d& _origin) const
{
	_angles.reserve(_points.size());
	for (const auto& iPoint : _points)
	{
		double vector_x = iPoint[0] - _origin[0];
		double vector_y = iPoint[1] - _origin[1];
		double norm = vector_x * vector_x + vector_y * vector_y;
		norm = sqrt(norm);
		double x = vector_x / norm;
		double angle = acos(x);
		if (vector_y < 0)
		{
			angle = 2 * M_PI - angle;
		}
		_angles.push_back(angle);
	}
}

void ConnectFrontNBack::FindConnectingBoundary(const Mesh& _cushion_surface, Mesh::HalfedgeHandle& _boundary) const
{
	std::vector<Mesh::HalfedgeHandle> initial_boundaries;
	MeshInfo find_boundary(&_cushion_surface);
	find_boundary.getInitialBoundaryHalfedge(initial_boundaries);
	// find the boundary with a larger z value(so this boundary should be connecting with the back part)
	{
		float vertex_z = -FLT_MAX;
		for (Mesh::HalfedgeHandle& iInitial : initial_boundaries)
		{
			// a point's z value of the boundary
			float z = _cushion_surface.point(_cushion_surface.from_vertex_handle(iInitial))[2];
			if (z > vertex_z)
			{
				vertex_z = z;
				_boundary = iInitial;
			}
		}
	}
}

void ConnectFrontNBack::FindConnectingBoundary(const Mesh& _cushion_surface, std::vector<Mesh::VertexHandle>& _boundary) const
{
	Mesh::HalfedgeHandle connecting_boundary;
	FindConnectingBoundary(_cushion_surface, connecting_boundary);
	
	{
		Mesh::HalfedgeHandle current_hafledge = connecting_boundary;
		do
		{
			Mesh::VertexHandle vertex = _cushion_surface.from_vertex_handle(current_hafledge);
			_boundary.push_back(vertex);

			// update for next loop
			current_hafledge = _cushion_surface.next_halfedge_handle(current_hafledge);
		} while (current_hafledge != connecting_boundary);
	}
}

void ConnectFrontNBack::FlattenPolylineLayers(const std::vector<std::vector<Mesh::Point>>& _polyline_layers, Mesh& _flattened_polyline) const
{
	for (auto& iLayer : _polyline_layers)
	{
		for (auto& iPolyline : iLayer)
		{
			_flattened_polyline.add_vertex(iPolyline);
		}
	}
}
