#include "CushionSurface.h"
#include "SelfIntersectionTest.h"
#include "FrameRelax.h"
#include "PullBackConnector.h"
#include "../MeshProcessing/MeshPlaneIntersection.h"
#include "../MeshProcessing/PolylinePlaneIntersection.h"
#include "../MeshProcessing/TransformMesh.h"
#include "../MeshProcessing/SubMesh.h"
#include "../CurveNSurface/CurveInfo.h"
#include "../CurveNSurface/ChangeFrame.h"
#include "../CurveNSurface/SamplePolylineFromCurve.h"
#include "../Utility/DivisionWithReminder.h"
#include "../Utility/ReadRhino.h"
#include <fstream>
#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>

using std::cout;
using std::endl;

CushionSurface::CushionSurface()
{
	head_tree_ = nullptr;
	cushion_ = nullptr;
	connector_pull_back_length_ = 0.0;
	cross_section_control_curve_degree_ = 3;
}

CushionSurface::~CushionSurface()
{
}

void CushionSurface::setHeadTree(const AabbTree* _head_tree)
{
	head_tree_ = _head_tree;
}

void CushionSurface::setCushionMesh(const Mesh* _cushion)
{
	cushion_ = _cushion;
}

void CushionSurface::freeCushionMesh()
{
	cushion_ = nullptr;
}

void CushionSurface::setCushionPolyline(const std::vector<Eigen::VectorXd>& _cushion_polyline)
{
	cushion_polyline_ = _cushion_polyline;
}

void CushionSurface::setRailByFittingCushionPolyline(int _segments, int _degree, int _smoothness)
{
	rail_.Fitting(cushion_polyline_, _segments, _degree, _smoothness);

	original_rail_ = rail_;
}

void CushionSurface::setRailByFittingCushionPolyline(const std::vector<int>& _parts_idx, int _degree, int _smoothness)
{
	rail_.Fitting(cushion_polyline_, _parts_idx, _degree, _smoothness);

	original_rail_ = rail_;
}

void CushionSurface::changeRailByFittingProjectedCushionPolylineOnFace(int _segments, int _degree, int _smoothness)
{
	std::vector<Eigen::VectorXd> projected_cushion_polyline;
	ProjectPolylineToFace(cushion_polyline_, projected_cushion_polyline);
	setRailByFitting(projected_cushion_polyline, _segments, _degree, _smoothness);
}

void CushionSurface::changeRailByFittingProjectedCushionPolylineOnFace(const std::vector<int>& _parts_idx, int _degree, int _smoothness)
{
	std::vector<Eigen::VectorXd> projected_cushion_polyline;
	ProjectPolylineToFace(cushion_polyline_, projected_cushion_polyline);
	setRailByFitting(projected_cushion_polyline, _parts_idx, _degree, _smoothness);
}

void CushionSurface::setRailByFitting(const std::vector<Eigen::VectorXd>& _polyline, int _segments, int _degree, int _smoothness)
{
	rail_.Fitting(_polyline, _segments, _degree, _smoothness);
}

void CushionSurface::setRailByFitting(const std::vector<Eigen::VectorXd>& _polyline, const std::vector<int>& _parts_idx, int _degree, int _smoothness)
{
	rail_.Fitting(_polyline, _parts_idx, _degree, _smoothness);

}

void CushionSurface::setRail(const BSpline& _rail)
{
	rail_.setRail(_rail);
}

void CushionSurface::setCrossSectionSurface(const BSpline& _half1_rest, const BSpline& _half1_point1_length, const BSpline& _tangent, const BSpline& _half2_point1_length, const BSpline& _half2_rest)
{
	cross_section_surface_.setControlCurve_Degree2(_half1_rest, _half1_point1_length, _tangent, _half2_point1_length, _half2_rest);
}

void CushionSurface::setLineSegmentAsCrossSectionSurface(double _line_radius)
{
	std::array<double, 2> interval{ 0, 1 };

	BSpline tangent_curve;
	tangent_curve.defineSpace(0, 0, 0, interval);
	Eigen::VectorXd tangent(2);
	tangent(0) = 1;
	tangent(1) = 0;
	tangent_curve.setControlPoints(std::vector<Eigen::VectorXd>{tangent});

	BSpline length_curve;
	length_curve.defineSpace(0, 0, 0, interval);
	Eigen::VectorXd length(1);
	length(0) = _line_radius;
	length_curve.setControlPoints(std::vector<Eigen::VectorXd>{length});

	BSpline half12_curve;
	half12_curve.defineSpace(0, 0, 0, interval);
	Eigen::VectorXd half12(2);
	half12(0) = -_line_radius;
	half12(1) = /*3*/0;
	half12_curve.setControlPoints(std::vector<Eigen::VectorXd>{half12});

	// this is the inside half 
	BSpline half22_curve;
	half22_curve.defineSpace(0, 0, 0, interval);
	Eigen::VectorXd half22(2);
	half22(0) = _line_radius;
	half22(1) = 0;
	half22_curve.setControlPoints(std::vector<Eigen::VectorXd>{half22});

	setCrossSectionSurface(half12_curve, length_curve, tangent_curve, length_curve, half22_curve);
}

void CushionSurface::FitCushion(const Mesh* _cushion_surface, 
	const std::vector<Eigen::VectorXd>& _rail_polygon, const std::vector<int>& _polygon_segment_parts_idx, const std::vector<double>& _key_cross_section_sample,
	int _sample_frames, int _sample_v, const Mesh& _connector_boundary,
	int _half1_degree, int _half2_degree)
{
	cout << "key cross section number = " << _key_cross_section_sample.size() << endl;

	// 1. rail		
	setCushionPolyline(_rail_polygon);
	int rail_segment = 20;
	int rail_degree = 3;
	int rail_smoothness = 2;
	//setRailByFittingCushionPolyline(_polygon_segment_parts_idx, rail_degree, rail_smoothness);// fit by parts so that control point is evenly spread
	setRailByFittingCushionPolyline(rail_segment, rail_degree, rail_smoothness);

	// 2. frames
	computeSampleFrames(20.0, _sample_frames);

	//3. cross section
	setCushionMesh(_cushion_surface);
	setKeyCrossSections_NewFrame_ChangableDegree(_key_cross_section_sample, _half1_degree, _half2_degree);

	int control_curve_degree = 3;
	ConstructCrossSectionSurface(control_curve_degree);
	freeCushionMesh();// original scanned cushion no longer needed

	// fix the boundary to connector
	setFixedBoundary(_connector_boundary, 100);

	//4. get cushion mesh
	//// directly construct a surface mesh without resample or do a resample if needed
	computeDiscreteCrossSections_NewFrames(_sample_v);
	constructSurfaceMesh();

}

void CushionSurface::computeSampleFrames(float _radius, int _num_sample_u)
{
	computeStandardFrames(_num_sample_u);
	FrameRelax relax_it;
	relax_it.GlobalCircleRelax(_radius, discrete_cushion_);

	// construct again with new frames
	computeNewFrames_AxisZ();
}

void CushionSurface::computeStandardSampleFrames(int _num_sample_u)
{
	computeStandardFrames(_num_sample_u);

	computeNewFrames_AxisZ();
}

void CushionSurface::computeInterpolatedFrameAtU(double _u, int _frame0_idx, int _frame1_idx, DiscreteCushionPiece& _cushion_piece) const
{
	_cushion_piece.u_ = _u;
	// interpolate the rotation
	double u0 = discrete_cushion_[_frame0_idx].u_;
	double u1 = discrete_cushion_[_frame1_idx].u_;
	if (_frame1_idx < _frame0_idx)// periodic case
	{
		assert(u0 > u1);
		if (_u >= u0 - 1e-10)
		{
			u1 += (rail_.getDomain()[1] - rail_.getDomain()[0]);// add domain length to u1, so that u in [u0,u1] 
		}
		else if (_u <= u1 + 1e-10)
		{
			u0 -= (rail_.getDomain()[1] - rail_.getDomain()[0]);// u0 subtract domain length, so that u in [u0,u1] 
		}
	}
	double t = getLocalLinearParameter(_u, u0, u1);
	assert(t > -1e-10 && t < 1.0 + 1e-10);
	//if (t < 0.0 || t > 1.0 || _frame0_idx == 999)
	//{
	//	cout << "t=" << t << endl;
	//	cout << "_u=" << _u << ", frame" << _frame0_idx << "_u=" << discrete_cushion_[_frame0_idx].u_ << ", frame" << _frame1_idx << "_u=" << discrete_cushion_[_frame1_idx].u_ << endl;
	//}
	_cushion_piece.rotation_ = discrete_cushion_[_frame0_idx].rotation_.slerp(t, discrete_cushion_[_frame1_idx].rotation_);

	Eigen::VectorXd rail_value;
	Eigen::VectorXd rail_tangent;
	rail_.getNormalPlane(_u, rail_tangent, rail_value);
	// standard frame
	_cushion_piece.frame_value_ = rail_value;
	_cushion_piece.frame_tangent_ = rail_tangent;
	_cushion_piece.frame_tangent_.z() = 0.0;
	_cushion_piece.frame_tangent_.normalize();
	getRailNormal_AxisZ(_cushion_piece.frame_tangent_, _cushion_piece.frame_normal_);
	_cushion_piece.frame_binormal_ = _cushion_piece.frame_tangent_.cross(_cushion_piece.frame_normal_);
	// relaxed frame
	_cushion_piece.rotated_tangent_ = _cushion_piece.rotation_ * _cushion_piece.frame_tangent_;
	getRailNormal_AxisZ(_cushion_piece.rotated_tangent_, _cushion_piece.new_frame_normal_);
	_cushion_piece.new_frame_binormal_ = _cushion_piece.rotated_tangent_.cross(_cushion_piece.new_frame_normal_);
}

void CushionSurface::ResampleFrames(int _num_sample_u)
{
	const std::array<double, 2>& domain_u = rail_.getDomain();
	double domain_length = domain_u[1] - domain_u[0];

	std::vector<DiscreteCushionPiece> discrete_cushion;// newly sampled cushion
	discrete_cushion.reserve(_num_sample_u);

	int frame_idx = 0;
	for (int iU = 0; iU < _num_sample_u; ++iU)
	{
		double u = (double)iU / (double)_num_sample_u * domain_length + domain_u[0];

		// after the while loop, we get frame_idx s.t. u between (frame_idx-1, frame_idx]
		while (frame_idx < discrete_cushion_.size() && discrete_cushion_[frame_idx].u_ < u - 1e-10)
		{
			++frame_idx;
		}

		// u between frame0, frame1
		int frame0 = frame_idx - 1;
		int frame1 = frame_idx;
		// mind the boundary case
		if (frame_idx == 0)
		{
			frame0 = discrete_cushion_.size() - 1;
		}
		else if (frame_idx == discrete_cushion_.size())
		{
			frame1 = 0;
		}

		//cout << "frame0, frame1 " << frame0 << "," << frame1 << endl;

		discrete_cushion.push_back(std::move(DiscreteCushionPiece()));
		computeInterpolatedFrameAtU(u, frame0, frame1, discrete_cushion[iU]);
	}

	discrete_cushion_ = std::move(discrete_cushion);
}

void CushionSurface::setFixedBoundary(const std::string& _original_boundary, int _num_interpolate_point, double _pull_back_connector)
{
	// read original boundary from a Rhino modelled file
	std::vector<Mesh::Point> outside_boundary;
	ReadRhino read_boundary;
	read_boundary.Read3DPoints_Brackets(_original_boundary, outside_boundary, 1);
	Mesh original_boundary;
	SubMesh convert_typer;
	convert_typer.PointCloudMeshFromPoints(outside_boundary, original_boundary);

	setFixedBoundary(original_boundary, _num_interpolate_point, _pull_back_connector);
}

void CushionSurface::setFixedHalf1Boundary(const std::string& _original_boundary, int _num_interpolate_point)
{
	// read original boundary from a Rhino modelled file
	std::vector<Mesh::Point> boundary;
	ReadRhino read_boundary;
	read_boundary.Read3DPoints_Brackets(_original_boundary, boundary, 1);

	// -------------------------------------------------- following code copied from setFixedBoundary -------------------------------------------
	
	// NOTE this sampling not creating exactly _num_interpolate_point of points......
	std::vector<int> sample_frames;
	std::vector<Mesh::Point> intersect_points;
	FindConnectorAndSampleFrameIntersectPoints(boundary, _num_interpolate_point, sample_frames, intersect_points);


	// set boundary

	// get the interpolation points and corresponding u values

	std::vector<Eigen::VectorXd> interpolate_points;
	std::vector<double> u;
	interpolate_points.reserve(sample_frames.size());
	u.reserve(sample_frames.size());
	for (int iSample = 0; iSample < sample_frames.size(); ++iSample)
	{
		u.push_back(discrete_cushion_[sample_frames[iSample]].u_);

		// 
		// 2. get the point in 2d local frame : similar to (CrossSection::compute2DLocalCoordinateOf3DPolyline)
		//
		// find transformation matrix
		Eigen::Matrix4d local_frame;
		getNewFrames(sample_frames[iSample], local_frame);
		Eigen::Matrix4d world_frame(Eigen::Matrix4d::Identity());
		Eigen::Matrix4d transform_to_local;// the transformation matrix
		ChangeFrame change_frame(world_frame, local_frame, transform_to_local);

		// transform
		Eigen::Vector4d point_world(intersect_points[iSample][0], intersect_points[iSample][1], intersect_points[iSample][2], 1.0);
		Eigen::Vector4d point_local = transform_to_local * point_world;
		interpolate_points.push_back(Eigen::Vector2d(point_local[1], point_local[2]));

	}

	//interpolate all the points : similar to (CrossSectionSurface::setControlCurveByInterpolation)
	//
	// check if first u is 0, if it is, we should remove it for interpolation algorithm not support this...
	if (u.size() > 0)
	{
		if (u[0] < 1e-10)
		{
			u.erase(u.begin());
			interpolate_points.erase(interpolate_points.begin());
		}
	}
	if (u.size() == 0)
	{
		cout << "ERROR!! CushionSurface::setFixedHalf1Boundary interpolation u size =0" << endl;
	}
	// -------------------------------------------------------------------------------------------------------------------------------------------------


	cross_section_surface_.setHalf1BoundaryControlCurveByInterpolation(interpolate_points, u);

	// at last, we adjust the key cross section control points accordingly
	for (auto& iCrossSection : key_cross_sections_)
	{
		Eigen::Vector2d half1_end_point;
		cross_section_surface_.getValue(iCrossSection.u_, 0.0, half1_end_point);
		iCrossSection.cross_section_.setHalf1EndPoint(half1_end_point);
	}

}

void CushionSurface::setFixedBoundary(const Mesh& _original_boundary, int _num_interpolate_point, double _pull_back_connector)
{
	// pull back connector according to target translation length
	Mesh pulled_back_boundary = _original_boundary;
	if (_pull_back_connector > 0.0)
	{
		PullBackConnector* pull_it_back = new PullBackConnector;
		// adjust translate_length here
		connector_pull_back_length_ = pull_it_back->TranslateAlongZ(_pull_back_connector, *this, _original_boundary, _num_interpolate_point);
		if (connector_pull_back_length_ > 0.0)
		{
			pull_it_back->TranslateMeshAlongZ(pulled_back_boundary, connector_pull_back_length_);
		}
		else
		{
			connector_pull_back_length_ = 0.0;
		}
		delete pull_it_back;
	}

	setFixedBoundary(pulled_back_boundary, _num_interpolate_point, false);

}

void CushionSurface::FindConnectorAndSampleFrameIntersectPoints(const std::vector<Mesh::Point>& _connector_polygon, int _reference_num_sample, 
	std::vector<int>& _sample_frames, std::vector<Mesh::Point>& _intersect_point) const
{
	// NOTE this sampling not creating exactly _reference_num_sample of points...
	int interval_size = discrete_cushion_.size() / _reference_num_sample;
	if (discrete_cushion_.size() <= _reference_num_sample)
	{
		interval_size = 1;
	}

	// get the intersect points and sample frame index
	_intersect_point.clear();
	_intersect_point.reserve(_reference_num_sample + 1);// a rough estimation
	_sample_frames.clear();
	_sample_frames.reserve(_reference_num_sample + 1);
	for (int iSample = 0; iSample < discrete_cushion_.size(); iSample += interval_size)
	{
		_sample_frames.push_back(iSample);

		//
		// compute intersection point of frame normal plane with boundary polyline
		// 
		Mesh::Point intersect_point;
		Mesh::Point plane_normal(discrete_cushion_[iSample].rotated_tangent_(0),
			discrete_cushion_[iSample].rotated_tangent_(1),
			discrete_cushion_[iSample].rotated_tangent_(2));
		Mesh::Point rail_point(discrete_cushion_[iSample].frame_value_(0),
			discrete_cushion_[iSample].frame_value_(1),
			discrete_cushion_[iSample].frame_value_(2));
		PolylinePlaneIntersection find_point;
		find_point.ComputeInterectionPoints_Polygon_ClosestOne(_connector_polygon, plane_normal, rail_point, intersect_point);

		_intersect_point.push_back(intersect_point);
	}
}

void CushionSurface::getHalf2Boundary(std::vector<Mesh::Point>& _boundary) const
{
	if (discrete_cushion_.size() == 0)
	{
		cout << "WARNING!!!!!!!!!!!!!!CushionSurface::getHalf2Boundary no discrete cushion construncted" << endl;
	}
	_boundary.resize(discrete_cushion_.size());
	for (int iSample = 0; iSample < discrete_cushion_.size(); ++iSample)
	{
		_boundary[iSample][0] = discrete_cushion_[iSample].discrete_cross_section_.back()[0];
		_boundary[iSample][1] = discrete_cushion_[iSample].discrete_cross_section_.back()[1];
		_boundary[iSample][2] = discrete_cushion_[iSample].discrete_cross_section_.back()[2];
	}
}

void CushionSurface::setKeyCrossSections(const std::vector<double>& _u, int _degree)
{
	if (cushion_ == nullptr)
	{
		cout << "[Warning from CushionSurface::setKeyCrossSections] cushion mesh not set" << endl;
		return;
	}

	key_cross_sections_.resize(_u.size());


	for (int iU = 0; iU < _u.size(); ++iU)
	{
		double u = _u[iU];

		// get cross section polylines

		std::vector<Mesh::Point> cross_section_polyline;

		Mesh::Point rail_tangent;
		Eigen::VectorXd rail_tangent_eigen;
		Mesh::Point rail_u;
		Eigen::VectorXd rail_u_eigen;
		rail_.getNormalPlane(u, rail_tangent_eigen, rail_u_eigen);
		rail_tangent[0] = rail_tangent_eigen[0];
		rail_tangent[1] = rail_tangent_eigen[1];
		rail_tangent[2] = rail_tangent_eigen[2];
		rail_u[0] = rail_u_eigen[0];
		rail_u[1] = rail_u_eigen[1];
		rail_u[2] = rail_u_eigen[2];

		// intersect the rail normal plane with the cushion
		FindCrossSectionPolyline(*cushion_, rail_tangent, rail_u, rail_u, cross_section_polyline);

		// fit the cross section polyline

		//get local frame
		Eigen::Matrix4d local_frame;
		Eigen::Vector3d frame_normal;
		getRailNormal_AxisZ(rail_tangent_eigen, frame_normal);
		CurveInfo get_rail_local_frame;
		get_rail_local_frame.getFrameGivenNormal_3DCurve(u, rail_.getRail(), frame_normal, local_frame);

		CrossSection cross_section;
		cross_section.FitOrigianl3DPolyline(cross_section_polyline, local_frame, _degree);

		key_cross_sections_[iU].u_ = u;
		key_cross_sections_[iU].cross_section_ = std::move(cross_section);
		key_cross_sections_[iU].cross_section_polyline_ = std::move(cross_section_polyline);
	}
}

void CushionSurface::setKeyCrossSections_NewFrame(const std::vector<int>& _i_frame, int _degree)
{
	if (cushion_ == nullptr)
	{
		cout << "[Warning from CushionSurface::setKeyCrossSections_NewFrame] cushion mesh not set" << endl;
		return;
	}

	key_cross_sections_.resize(_i_frame.size());


	for (int iFrame = 0; iFrame < _i_frame.size(); ++iFrame)
	{
		// get cross section polylines by cutting it with frame normal plane

		std::vector<Mesh::Point> cross_section_polyline;
		getOriginalCrossSection(_i_frame[iFrame], cross_section_polyline);
		//Mesh::Point frame_tangent(discrete_cushion_[_i_frame[iFrame]].rotated_tangent_(0), 
		//						  discrete_cushion_[_i_frame[iFrame]].rotated_tangent_(1), 
		//						  discrete_cushion_[_i_frame[iFrame]].rotated_tangent_(2));
		//Mesh::Point frame_value(discrete_cushion_[_i_frame[iFrame]].frame_value_(0),
		//						discrete_cushion_[_i_frame[iFrame]].frame_value_(1), 
		//						discrete_cushion_[_i_frame[iFrame]].frame_value_(2));
		//
		//// intersect the rail normal plane with the cushion
		//FindCrossSectionPolyline(*cushion_, frame_tangent, frame_value, frame_value, cross_section_polyline);

		// fit the cross section polyline

		//get local frame
		Eigen::Matrix4d local_frame;
		getNewFrames(_i_frame[iFrame], local_frame);

		CrossSection cross_section;
		cross_section.FitOrigianl3DPolyline(cross_section_polyline, local_frame, _degree);
		
		key_cross_sections_[iFrame].u_ = discrete_cushion_[_i_frame[iFrame]].u_;
		key_cross_sections_[iFrame].cross_section_ = std::move(cross_section);
		key_cross_sections_[iFrame].cross_section_polyline_ = std::move(cross_section_polyline);
	}
}

void CushionSurface::setKeyCrossSections_NewFrame(const std::vector<double>& _u, int _degree)
{
	std::vector<int> i_frame;
	ClosestFrameIdxToU(_u, i_frame);

	setKeyCrossSections_NewFrame(i_frame, _degree);

}

void CushionSurface::setKeyCrossSections_NewFrame_ChangableDegree(const std::vector<double>& _u, int _degree_half1, int _degree_half2)
{
	// transfer the u value to its closest frame
	std::vector<int> i_frame;
	ClosestFrameIdxToU(_u, i_frame);

	if (cushion_ == nullptr)
	{
		cout << "[Warning from CushionSurface::setKeyCrossSections_NewFrame_ChangableDegree] cushion mesh not set" << endl;
		return;
	}

	key_cross_sections_.resize(i_frame.size());

	for (int iFrame = 0; iFrame < i_frame.size(); ++iFrame)
	{
		// get cross section polylines by cutting it with frame normal plane
		std::vector<Mesh::Point> cross_section_polyline;
		getOriginalCrossSection(i_frame[iFrame], cross_section_polyline);

		// fit the cross section polyline
		
		//get local frame
		Eigen::Matrix4d local_frame;
		getNewFrames(i_frame[iFrame], local_frame);

		CrossSection cross_section;
		cross_section.FitOrigianl3DPolyline_ChangeableDegree(cross_section_polyline, local_frame, _degree_half1, _degree_half2);

		key_cross_sections_[iFrame].u_ = discrete_cushion_[i_frame[iFrame]].u_;
		key_cross_sections_[iFrame].cross_section_ = std::move(cross_section);
		key_cross_sections_[iFrame].cross_section_polyline_ = std::move(cross_section_polyline);
	}
}

void CushionSurface::setKeyCrossSections_NewFrame_ChangableDegree_FindHalf2BeginIndex(const std::vector<double>& _u, const std::vector<Mesh::Point>& _half2_endpoints, int _degree_half1, int _degree_half2)
{
	// transfer the u value to its closest frame
	std::vector<int> i_frame;
	ClosestFrameIdxToU(_u, i_frame);

	if (cushion_ == nullptr)
	{
		cout << "[Warning from CushionSurface::setKeyCrossSections_NewFrame_ChangableDegree_FindHalf2BeginIndex] cushion mesh not set" << endl;
		return;
	}

	key_cross_sections_.resize(i_frame.size());

	for (int iFrame = 0; iFrame < i_frame.size(); ++iFrame)
	{
		// get cross section polylines by cutting it with frame normal plane
		std::vector<Mesh::Point> cross_section_polyline;
		getOriginalCrossSection(i_frame[iFrame], cross_section_polyline);
		
		// re-orient cross_section_polyline, from half1 to half2
		// 
		// the intersect point of the _half2_endpoint curve with the frame normal plane
		// since there might be multiple intersect points, I think we can get the right one by finding the closest one to the frame value
		Mesh::Point frame_tangent(discrete_cushion_[i_frame[iFrame]].rotated_tangent_(0),
			discrete_cushion_[i_frame[iFrame]].rotated_tangent_(1),
			discrete_cushion_[i_frame[iFrame]].rotated_tangent_(2));
		Mesh::Point frame_value(discrete_cushion_[i_frame[iFrame]].frame_value_(0),
			discrete_cushion_[i_frame[iFrame]].frame_value_(1),
			discrete_cushion_[i_frame[iFrame]].frame_value_(2));
		Mesh::Point correspond_half2_endpoint;
		PolylinePlaneIntersection get_correspond_half2_endpoint;
		get_correspond_half2_endpoint.ComputeInterectionPoints_Polygon_ClosestOne(_half2_endpoints, frame_tangent, frame_value, correspond_half2_endpoint);
		if ((cross_section_polyline.back() - correspond_half2_endpoint).norm() > (cross_section_polyline[0] - correspond_half2_endpoint).norm())
		{
			// re-orient cross_section_polyline, from half1 to half2
			std::reverse(cross_section_polyline.begin(), cross_section_polyline.end());
		}

		// find the vertex in polyline that's closet to the frame value
		double min_distance = FLT_MAX;
		std::vector<Mesh::Point>::const_iterator min_iter = cross_section_polyline.cbegin();
		if (cross_section_polyline.size() > 0)
		{
			for (auto iterVertex = cross_section_polyline.cbegin(); iterVertex != cross_section_polyline.end(); ++iterVertex)
			{
				double distance = (*iterVertex - frame_value).sqrnorm();
				if (distance < min_distance)
				{
					min_distance = distance;
					min_iter = iterVertex;
				}
			}
		}
		int half2_begin_idx = std::distance(cross_section_polyline.cbegin(), min_iter);

		// fit the cross section polyline

		//get local frame
		Eigen::Matrix4d local_frame;
		getNewFrames(i_frame[iFrame], local_frame);

		CrossSection cross_section;
		// set half2 begin index as cross_section_polyline.size()/2
		cross_section.FitOrigianl3DPolyline_ChangeableDegree(cross_section_polyline, half2_begin_idx, local_frame, _degree_half1, _degree_half2);
		//cross_section.FitOrigianl3DPolyline(cross_section_polyline, half2_begin_idx, local_frame, _degree_half1);

		//if (iFrame >= 9 && iFrame <=12)
		//{
		//	cout << "u=" << discrete_cushion_[i_frame[iFrame]].u_ << ", at frame" << iFrame << endl;
		//}

		key_cross_sections_[iFrame].u_ = discrete_cushion_[i_frame[iFrame]].u_;
		key_cross_sections_[iFrame].cross_section_ = std::move(cross_section);
		key_cross_sections_[iFrame].cross_section_polyline_ = std::move(cross_section_polyline);
		//key_cross_sections_[iFrame].cross_section_polyline_ = std::vector<Mesh::Point>{ cross_section_polyline.back() };
	}
}

void CushionSurface::setInitialCustomKeyCrossSections(const std::vector<int>& _i_frame, const Mesh& _original_cushion, const std::vector<Mesh::Point>& _original_rail_polyline, int _degree)
{
	key_cross_sections_.resize(_i_frame.size());
	for (int iFrame = 0; iFrame < _i_frame.size(); ++iFrame)
	{
		std::vector<Mesh::Point> cross_section_polyline;
		int center_idx = 0;
		getOriginalCrossSection_OffAxis(_i_frame[iFrame], _original_cushion, _original_rail_polyline, cross_section_polyline, center_idx);

		// fit stretched cross section(i.e. only fix boundary while other control points to origin)
		
		//get local frame
		Eigen::Matrix4d local_frame;
		getNewFrames(_i_frame[iFrame], local_frame);

		CrossSection cross_section;
		cross_section.Fit3DPolyline_FixBoundary(cross_section_polyline, center_idx, local_frame, _degree);

		key_cross_sections_[iFrame].u_ = discrete_cushion_[_i_frame[iFrame]].u_;
		key_cross_sections_[iFrame].cross_section_ = std::move(cross_section);
		key_cross_sections_[iFrame].cross_section_polyline_ = std::move(cross_section_polyline);
		//key_cross_sections_[iFrame].rail_point_ = correspond_original_rail_point;
	}
}

void CushionSurface::setInitialCustomKeyCrossSections(const std::vector<double>& _u, int _degree)
{
	std::vector<int> i_frame;
	ClosestFrameIdxToU(_u, i_frame);

	std::vector<Mesh::Point> original_rail_polyline;
	SamplePolylineFromCurve get_original_rail_poylgon;
	get_original_rail_poylgon.getPolylineAsPointVector3D<Mesh::Point>(original_rail_.getRail(), 10000, original_rail_polyline, original_rail_.getDomain());

	// using the fitted surface mesh from the original cushion scan as reference 
	setInitialCustomKeyCrossSections(i_frame, original_fitted_surface_, original_rail_polyline, _degree);
}

void CushionSurface::setInitialCustomKeyCrossSections_ChangableDegree(const std::vector<double>& _u, int _degree_half1, int _degree_half2)
{
	std::vector<int> i_frame;
	ClosestFrameIdxToU(_u, i_frame);

	std::vector<Mesh::Point> original_rail_polyline;
	SamplePolylineFromCurve get_original_rail_poylgon;
	get_original_rail_poylgon.getPolylineAsPointVector3D<Mesh::Point>(original_rail_.getRail(), 10000, original_rail_polyline, original_rail_.getDomain());

	key_cross_sections_.resize(i_frame.size());
	for (int iFrame = 0; iFrame < i_frame.size(); ++iFrame)
	{
		std::vector<Mesh::Point> cross_section_polyline;
		int center_idx = 0;
		getOriginalCrossSection_OffAxis(i_frame[iFrame], original_fitted_surface_, original_rail_polyline, cross_section_polyline, center_idx);

		// fit stretched cross section(i.e. only fix boundary while other control points to origin)

		//get local frame
		Eigen::Matrix4d local_frame;
		getNewFrames(i_frame[iFrame], local_frame);

		CrossSection cross_section;
		cross_section.Fit3DPolyline_FixBoundary_ChangeableDegree(cross_section_polyline, center_idx, local_frame, _degree_half1, _degree_half2);

		key_cross_sections_[iFrame].u_ = discrete_cushion_[i_frame[iFrame]].u_;
		key_cross_sections_[iFrame].cross_section_ = std::move(cross_section);
		key_cross_sections_[iFrame].cross_section_polyline_ = std::move(cross_section_polyline);
		//key_cross_sections_[iFrame].rail_point_ = correspond_original_rail_point;
	}

}

void CushionSurface::changeKeyCrossSections(const std::vector<CrossSection>& _key_cross_sections)
{
	int num_section = _key_cross_sections.size();
	if (num_section == key_cross_sections_.size())
	{
		for (int iSection = 0; iSection < num_section; ++iSection)
		{
			key_cross_sections_[iSection].cross_section_ = _key_cross_sections[iSection];
		}
	}
	else
	{
		cout << "[WARNING CushionSurface::changeKeyCrossSections] can't modify cushion because key cross section number not match" << endl;
	}
}

void CushionSurface::setKeyCrossSectionFitAnotherCushion(const Mesh& _cushion, const std::vector<Mesh::Point>& _rail_polyline, int _degree)
{
	FindKeyCrossSectionFrameId();
	for (int iFrame = 0; iFrame < key_cross_sections_.size(); ++iFrame)
	{
		int frame_id = key_cross_sections_[iFrame].frame_id_;
		std::vector<Mesh::Point> cross_section_polyline;
		int center_idx = 0;
		getOriginalCrossSection_OffAxis(frame_id, _cushion, _rail_polyline, cross_section_polyline, center_idx);

		//get local frame
		Eigen::Matrix4d local_frame;
		getNewFrames(frame_id, local_frame);

		CrossSection cross_section;
		cross_section.FitOrigianl3DPolyline_ChangeableDegree(cross_section_polyline, center_idx, local_frame, _degree, _degree);

		key_cross_sections_[iFrame].u_ = discrete_cushion_[frame_id].u_;
		key_cross_sections_[iFrame].cross_section_ = std::move(cross_section);
		key_cross_sections_[iFrame].cross_section_polyline_ = std::move(cross_section_polyline);
	}
}

void CushionSurface::setKeyCrossSection(const std::vector<CrossSection>& _key_cross_sections)
{
	if (_key_cross_sections.size() != key_cross_sections_.size())
	{
		cout << "WARNING!!!!!!!!!!!!!!!!!!!!!!CushionSurface::setKeyCrossSection key cross section size not match" << endl;
	}
	for (int iKey = 0; iKey < _key_cross_sections.size(); ++iKey)
	{
		key_cross_sections_[iKey].cross_section_ = _key_cross_sections[iKey];
	}
}

void CushionSurface::ConstructCrossSectionSurface(int _degree)
{
	int num_section = key_cross_sections_.size();
	std::vector<const CrossSection*> cross_sections;
	cross_sections.reserve(num_section);
	std::vector<double> u(num_section);
	for (int iSection = 0; iSection < num_section; ++iSection)
	{
		cross_sections.push_back(&(key_cross_sections_[iSection].cross_section_));
		u[iSection] = key_cross_sections_[iSection].u_;
	}
	cross_section_surface_.setControlCurveByInterpolation(cross_sections, u, _degree, rail_.getDomain());
	cross_section_control_curve_degree_ = _degree;
}


void CushionSurface::getValue_FaceNormal(double _u, double _v, Mesh::Point& _value) const
{
	Eigen::Vector3d rail_value;
	Eigen::Vector3d frame_tangent;
	Eigen::Vector3d frame_normal;

	getRailFrame_FromFace(_u, rail_value, frame_tangent, frame_normal);
	Eigen::Vector3d frame_binormal = frame_tangent.cross(frame_normal);

	Eigen::Vector3d cross_section_value;
	cross_section_surface_.getValue_Degree2(_u, _v, cross_section_value);

	// surface formula
	Eigen::Vector3d value =
		rail_value +
		frame_tangent * cross_section_value[0] +
		frame_normal * cross_section_value[1] +
		frame_binormal * cross_section_value[2];

	_value[0] = value[0];
	_value[1] = value[1];
	_value[2] = value[2];
}

void CushionSurface::getValue(double _u, double _v, const Eigen::Vector3d& _rail_normal, Mesh::Point& _value) const
{
	Eigen::Matrix4d local_frame;
	CurveInfo get_rail_local_frame;
	get_rail_local_frame.getFrameGivenNormal_3DCurve(_u, rail_.getRail(), _rail_normal, local_frame);

	Eigen::Vector3d cross_section_value;
	cross_section_surface_.getValue_Degree2(_u, _v, cross_section_value);

	Eigen::Vector4d cross_section_value_4d(Eigen::Vector4d::Ones());
	cross_section_value_4d[0] = cross_section_value[0];
	cross_section_value_4d[1] = cross_section_value[1];
	cross_section_value_4d[2] = cross_section_value[2];
	Eigen::Vector3d value = (local_frame * cross_section_value_4d).block(0, 0, 3, 1);
	_value[0] = value[0];
	_value[1] = value[1];
	_value[2] = value[2];
}

void CushionSurface::getValue(double _u, double _v, 
	const Eigen::Vector3d& _rail_value, const Eigen::Vector3d& _frame_tangent, const Eigen::Vector3d& _frame_normal,
	Eigen::Vector3d& _value) const
{
	Eigen::Vector3d frame_binormal = _frame_tangent.cross(_frame_normal);

	getValue(_u, _v, _rail_value, _frame_tangent, _frame_normal, frame_binormal, _value);
}

void CushionSurface::getValue(double _u, double _v, 
	const Eigen::Vector3d& _rail_value, const Eigen::Vector3d& _frame_tangent, const Eigen::Vector3d& _frame_normal, const Eigen::Vector3d& _frame_binormal,
	Eigen::Vector3d& _value) const
{
	Eigen::Vector2d cross_section_value;
	cross_section_surface_.getValue(_u, _v, cross_section_value);

	// surface formula
	_value =
		_rail_value +
		_frame_normal * cross_section_value[0] +
		_frame_binormal * cross_section_value[1];
}

void CushionSurface::getRailValue(int _frame, Eigen::Vector3d& _value) const
{
	_value = discrete_cushion_[_frame].frame_value_;
}

void CushionSurface::UpdateDerivativeInfo()
{
	cross_section_surface_.UpdateDerivativeInfo();
}

void CushionSurface::getDerivative_DiscreteCushionAndCrossSectionNormal_KeyControlPoints(
	std::vector<std::vector<std::vector<std::array<Eigen::Matrix3Xd, 2>>>>& _half1_rest, 
	std::vector<std::vector<std::vector<std::array<Eigen::Matrix3Xd, 2>>>>& _half2_rest,
	std::vector<std::vector<Eigen::Matrix3Xd>>& _half1_point1_length,
	std::vector<std::vector<Eigen::Matrix3Xd>>& _tangent_angle, 
	std::vector<std::vector<Eigen::Matrix3Xd>>& _half2_point1_length, 
	std::vector<std::vector<std::array<Eigen::Matrix3Xd, 2>>>& _point0, 
	std::vector<std::vector<std::vector<std::array<Eigen::Matrix3Xd, 2>>>>& _normal_half1_rest, 
	std::vector<std::vector<std::vector<std::array<Eigen::Matrix3Xd, 2>>>>& _normal_half2_rest, 
	std::vector<std::vector<Eigen::Matrix3Xd>>& _normal_half1_point1_length, 
	std::vector<std::vector<Eigen::Matrix3Xd>>& _normal_tangent_angle, 
	std::vector<std::vector<Eigen::Matrix3Xd>>& _normal_half2_point1_length, 
	std::vector<std::vector<std::array<Eigen::Matrix3Xd, 2>>>& _normal_point0) const
{
	std::vector<std::vector<std::vector<std::array<Eigen::Matrix2Xd, 2>>>> half1_rest;
	std::vector<std::vector<std::vector<std::array<Eigen::Matrix2Xd, 2>>>> half2_rest;
	std::vector<std::vector<Eigen::Matrix2Xd>> half1_point1_length;
	std::vector<std::vector<Eigen::Matrix2Xd>> tangent_angle;
	std::vector<std::vector<Eigen::Matrix2Xd>> half2_point1_length;
	std::vector<std::vector<std::array<Eigen::Matrix2Xd, 2>>> point0;
	std::vector<std::vector<std::vector<std::array<Eigen::Matrix2Xd, 2>>>> normal_half1_rest;
	std::vector<std::vector<std::vector<std::array<Eigen::Matrix2Xd, 2>>>> normal_half2_rest;
	std::vector<std::vector<Eigen::Matrix2Xd>> normal_half1_point1_length;
	std::vector<std::vector<Eigen::Matrix2Xd>> normal_tangent_angle;
	std::vector<std::vector<Eigen::Matrix2Xd>> normal_half2_point1_length;
	std::vector<std::vector<std::array<Eigen::Matrix2Xd, 2>>> normal_point0;
	cross_section_surface_.findDerivative_CrossSectionSurfaceAndNormal_KeyControlPoints(
		half1_rest, half2_rest, half1_point1_length, tangent_angle, half2_point1_length, point0, 
		normal_half1_rest, normal_half2_rest, normal_half1_point1_length, normal_tangent_angle, normal_half2_point1_length, normal_point0);

	_half1_rest.resize(half1_rest.size(), std::vector<std::vector<std::array<Eigen::Matrix3Xd, 2>>>(u_values_.size(), std::vector<std::array<Eigen::Matrix3Xd, 2>>(v_values_.size())));
	_half2_rest.resize(half2_rest.size(), std::vector<std::vector<std::array<Eigen::Matrix3Xd, 2>>>(u_values_.size(), std::vector<std::array<Eigen::Matrix3Xd, 2>>(v_values_.size())));
	_half1_point1_length.resize(u_values_.size(), std::vector<Eigen::Matrix3Xd>(v_values_.size()));
	_tangent_angle.resize(u_values_.size(), std::vector<Eigen::Matrix3Xd>(v_values_.size()));
	_half2_point1_length.resize(u_values_.size(), std::vector<Eigen::Matrix3Xd>(v_values_.size()));
	_point0.resize(u_values_.size(), std::vector<std::array<Eigen::Matrix3Xd, 2>>(v_values_.size()));

	_normal_half1_rest.resize(half1_rest.size(), std::vector<std::vector<std::array<Eigen::Matrix3Xd, 2>>>(u_values_.size(), std::vector<std::array<Eigen::Matrix3Xd, 2>>(v_values_.size())));
	_normal_half2_rest.resize(half2_rest.size(), std::vector<std::vector<std::array<Eigen::Matrix3Xd, 2>>>(u_values_.size(), std::vector<std::array<Eigen::Matrix3Xd, 2>>(v_values_.size())));
	_normal_half1_point1_length.resize(u_values_.size(), std::vector<Eigen::Matrix3Xd>(v_values_.size()));
	_normal_tangent_angle.resize(u_values_.size(), std::vector<Eigen::Matrix3Xd>(v_values_.size()));
	_normal_half2_point1_length.resize(u_values_.size(), std::vector<Eigen::Matrix3Xd>(v_values_.size()));
	_normal_point0.resize(u_values_.size(), std::vector<std::array<Eigen::Matrix3Xd, 2>>(v_values_.size()));

	for (int iU = 0; iU < u_values_.size(); ++iU)
	{
		for (int iV = 0; iV < v_values_.size(); ++iV)
		{
			for (int iRest1 = 0; iRest1 < _half1_rest.size(); ++iRest1)
			{
				for (int iCoordinate = 0; iCoordinate < 2; ++iCoordinate)
				{
					_half1_rest[iRest1][iU][iV][iCoordinate] = discrete_cushion_[iU].new_frame_normal_ * half1_rest[iRest1][iU][iV][iCoordinate].row(0) +
						discrete_cushion_[iU].new_frame_binormal_ * half1_rest[iRest1][iU][iV][iCoordinate].row(1);
					_normal_half1_rest[iRest1][iU][iV][iCoordinate] = discrete_cushion_[iU].new_frame_normal_ * normal_half1_rest[iRest1][iU][iV][iCoordinate].row(0) +
						discrete_cushion_[iU].new_frame_binormal_ * normal_half1_rest[iRest1][iU][iV][iCoordinate].row(1);
				}
			}

			for (int iRest2 = 0; iRest2 < _half2_rest.size(); ++iRest2)
			{
				for (int iCoordinate = 0; iCoordinate < 2; ++iCoordinate)
				{
					_half2_rest[iRest2][iU][iV][iCoordinate] = discrete_cushion_[iU].new_frame_normal_ * half2_rest[iRest2][iU][iV][iCoordinate].row(0) +
						discrete_cushion_[iU].new_frame_binormal_ * half2_rest[iRest2][iU][iV][iCoordinate].row(1);
					_normal_half2_rest[iRest2][iU][iV][iCoordinate] = discrete_cushion_[iU].new_frame_normal_ * normal_half2_rest[iRest2][iU][iV][iCoordinate].row(0) +
						discrete_cushion_[iU].new_frame_binormal_ * normal_half2_rest[iRest2][iU][iV][iCoordinate].row(1);
				}
			}

			for (int iCoordinate = 0; iCoordinate < 2; ++iCoordinate)
			{
				_point0[iU][iV][iCoordinate] = discrete_cushion_[iU].new_frame_normal_ * point0[iU][iV][iCoordinate].row(0) +
					discrete_cushion_[iU].new_frame_binormal_ * point0[iU][iV][iCoordinate].row(1);
				_normal_point0[iU][iV][iCoordinate] = discrete_cushion_[iU].new_frame_normal_ * normal_point0[iU][iV][iCoordinate].row(0) +
					discrete_cushion_[iU].new_frame_binormal_ * normal_point0[iU][iV][iCoordinate].row(1);
			}

			_half1_point1_length[iU][iV] = discrete_cushion_[iU].new_frame_normal_ * half1_point1_length[iU][iV].row(0) +
				discrete_cushion_[iU].new_frame_binormal_ * half1_point1_length[iU][iV].row(1);
			_tangent_angle[iU][iV] = discrete_cushion_[iU].new_frame_normal_ * tangent_angle[iU][iV].row(0) +
				discrete_cushion_[iU].new_frame_binormal_ * tangent_angle[iU][iV].row(1);
			_half2_point1_length[iU][iV] = discrete_cushion_[iU].new_frame_normal_ * half2_point1_length[iU][iV].row(0) +
				discrete_cushion_[iU].new_frame_binormal_ * half2_point1_length[iU][iV].row(1);

			_normal_half1_point1_length[iU][iV] = discrete_cushion_[iU].new_frame_normal_ * normal_half1_point1_length[iU][iV].row(0) +
				discrete_cushion_[iU].new_frame_binormal_ * normal_half1_point1_length[iU][iV].row(1);
			_normal_tangent_angle[iU][iV] = discrete_cushion_[iU].new_frame_normal_ * normal_tangent_angle[iU][iV].row(0) +
				discrete_cushion_[iU].new_frame_binormal_ * normal_tangent_angle[iU][iV].row(1);
			_normal_half2_point1_length[iU][iV] = discrete_cushion_[iU].new_frame_normal_ * normal_half2_point1_length[iU][iV].row(0) +
				discrete_cushion_[iU].new_frame_binormal_ * normal_half2_point1_length[iU][iV].row(1);

		}

	}


}

void CushionSurface::getKeyCrossSectionsPolyline(std::vector<std::vector<Mesh::Point>>& _cross_section_polyline) const
{
	int num_key_sections = key_cross_sections_.size();
	_cross_section_polyline.resize(num_key_sections);
	for (int iSection = 0; iSection < num_key_sections; ++iSection)
	{
		_cross_section_polyline[iSection] = key_cross_sections_[iSection].cross_section_polyline_;
	}
}

void CushionSurface::getKeyCrossSectionsRailPoint(std::vector<Mesh::Point>& _rail_point) const
{
	int num_key_sections = key_cross_sections_.size();
	_rail_point.resize(num_key_sections);
	for (int iSection = 0; iSection < num_key_sections; ++iSection)
	{
		//_rail_point[iSection] = key_cross_sections_[iSection].rail_point_;
	}
}

void CushionSurface::getKeyCrossSections(std::vector<CrossSection>& _key_cross_sections) const
{
	int num_key_sections = key_cross_sections_.size();
	_key_cross_sections.resize(num_key_sections);
	for (int iSection = 0; iSection < num_key_sections; ++iSection)
	{
		_key_cross_sections[iSection] = (key_cross_sections_[iSection].cross_section_);
	}
}

void CushionSurface::getKeyCrossSectionsU(std::vector<double>& _key_cross_sections_u) const
{
	int num_key_sections = key_cross_sections_.size();
	_key_cross_sections_u.resize(num_key_sections);
	for (int iSection = 0; iSection < num_key_sections; ++iSection)
	{
		_key_cross_sections_u[iSection] = (key_cross_sections_[iSection].u_);
	}
}

void CushionSurface::getKeyCrossSectionsMesh_3D(std::vector<Mesh>& _curves) const
{
	int num_sample_v = 100;
	int num_key = key_cross_sections_.size();
	// find frame idx from the u value
	std::vector<double> u(num_key);
	for (int iKey = 0; iKey < num_key; ++iKey)
	{
		u[iKey] = key_cross_sections_[iKey].u_;
	}
	std::vector<int> frame;
	ClosestFrameIdxToU(u, frame);

	assert(frame.size() == num_key);
	_curves.resize(num_key);
	for (int iKey = 0; iKey < num_key; ++iKey)
	{
		SubMesh convert_to_mesh;
		convert_to_mesh.PointCloudMeshFromPoints(discrete_cushion_[frame[iKey]].discrete_cross_section_, _curves[iKey]);
	}
}

void CushionSurface::getKeyCrossSection3DBeizerPoints(Mesh& _points) const
{
	int num_key = key_cross_sections_.size();
	// find frame idx from the u value
	std::vector<double> u(num_key);
	for (int iKey = 0; iKey < num_key; ++iKey)
	{
		u[iKey] = key_cross_sections_[iKey].u_;
	}
	std::vector<int> frame;
	ClosestFrameIdxToU(u, frame);

	for (int iKey = 0; iKey < num_key; ++iKey)
	{
		Mesh bezier_points2d;
		Mesh bezier_curve;//not used
		key_cross_sections_[iKey].cross_section_.getBezierCurve(bezier_curve, bezier_points2d);


		for (auto& iVertex : bezier_points2d.vertices())
		{
			Eigen::Vector3d point3d = discrete_cushion_[frame[iKey]].frame_value_ +
				discrete_cushion_[frame[iKey]].new_frame_normal_ * bezier_points2d.point(iVertex)[1] +
				discrete_cushion_[frame[iKey]].new_frame_binormal_ * bezier_points2d.point(iVertex)[2];
			_points.add_vertex(Mesh::Point(point3d[0], point3d[1], point3d[2]));
		}
	}
}

int CushionSurface::getKeyCrossSectionsNumber() const
{
	return key_cross_sections_.size();
}

const CrossSectionSurface& CushionSurface::getCrossSectionSurface() const
{
	return cross_section_surface_;
}

void CushionSurface::getRailNormal_FromFace(const Mesh::Point& _rail_tangent, const Mesh::Point& _rail_value, Eigen::Vector3d& _normal) const
{
	if (head_tree_ == nullptr)
	{
		cout << "[Warning from CushionSurface::getRailNormal_FromFace] head AABB Tree not set" << endl;
		return;
	}

	uint face_id = 0;
	const Mesh* head = head_tree_->getPrimitives();
	head_tree_->LineMeshIntersection(_rail_value, Mesh::Point(0.0, 0.0, -1.0), Mesh::Point(), face_id);
	Mesh::Point face_normal(0, 0, 0);

	// face normal average with vertex normal
	for (Mesh::VertexHandle iVertex : head->fv_range(head->face_handle(face_id)))
	{
		face_normal += head->normal(iVertex);
	}
	face_normal /= 3.0;

	Mesh::Point normal = face_normal.cross(_rail_tangent).normalized();
	_normal = Eigen::Map<Eigen::Vector3d>(normal.data());// convert type to Eigen

}

void CushionSurface::getRailFrame_FromFace(double _u, Eigen::Vector3d& _rail_value, Eigen::Vector3d& _rail_tangent, Eigen::Vector3d& _frame_normal) const
{
	Eigen::VectorXd rail_tangent;
	Eigen::VectorXd rail_value;
	rail_.getNormalPlane(_u, rail_tangent, rail_value);

	// get local frame
	_rail_value = rail_value;
	_rail_tangent = rail_tangent;
	_rail_tangent.normalize();
	Mesh::Point rail_tangent_openmesh(rail_tangent[0], rail_tangent[1], rail_tangent[2]);
	Mesh::Point rail_value_openmesh(rail_value[0], rail_value[1], rail_value[2]);
	getRailNormal_FromFace(rail_tangent_openmesh, rail_value_openmesh, _frame_normal);
}

void CushionSurface::getRailNormal_AxisZ(const Eigen::Vector3d& _rail_tangent, Eigen::Vector3d& _normal) const
{
	//TODO: determine normal direction! notice cross product sequence makes a difference!
	_normal = Eigen::Vector3d::UnitZ().cross(_rail_tangent).normalized();
}

void CushionSurface::getRailFrame_AxisZ(double _u, Eigen::Vector3d& _rail_value, Eigen::Vector3d& _rail_tangent, Eigen::Vector3d& _frame_normal) const
{
	Eigen::VectorXd rail_tangent;
	Eigen::VectorXd rail_value;
	rail_.getNormalPlane(_u, rail_tangent, rail_value);

	// get local frame
	_rail_value = rail_value;
	_rail_tangent = rail_tangent;

	_rail_tangent.z() = 0.0;
	_rail_tangent.normalize();

	getRailNormal_AxisZ(rail_tangent, _frame_normal);
}

void CushionSurface::getFrameGivenTangent_AxisZ(double _u, const Eigen::Vector3d& _frame_tangent, Eigen::Vector3d& _rail_value, Eigen::Vector3d& _frame_normal, Eigen::Vector3d& _frame_binormal) const
{
	Eigen::VectorXd rail_value;
	rail_.getValue(_u, rail_value);
	_rail_value = rail_value;
	getRailNormal_AxisZ(_frame_tangent, _frame_normal);
	_frame_binormal = _frame_tangent.cross(_frame_normal);

}

void CushionSurface::getRailFrameWithSmoothedFaceNormals(int _num_sample_u,
	std::vector < Eigen::Vector3d>& _rail_value,
	std::vector < Eigen::Vector3d>& _rail_tangent,
	std::vector < Eigen::Vector3d>& _frame_normal,
	std::vector < Eigen::Vector3d>& _frame_binormal,
	int _neighbor) const
{
	_rail_value.resize(_num_sample_u);
	_rail_tangent.resize(_num_sample_u);
	_frame_normal.resize(_num_sample_u);
	_frame_binormal.resize(_num_sample_u);

	const std::array<double, 2>& domain_u = rail_.getDomain();

	for (int iU = 0; iU < _num_sample_u; ++iU)
	{
		double u = (double)iU / (double)_num_sample_u * (domain_u[1] - domain_u[0]) + domain_u[0];
		getRailFrame_FromFace(u, _rail_value[iU], _rail_tangent[iU], _frame_normal[iU]);
	}

	// average normals
	std::vector < Eigen::Vector3d>& frame_normal(_frame_normal);// a copy
	DivisionWithReminder get_index;
	for (int iU = 0; iU < _num_sample_u; ++iU)
	{
		for (int iNeighbor = 1; iNeighbor <= _neighbor; ++iNeighbor)
		{
			_frame_normal[iU] += (frame_normal[get_index.getReminder(iU - iNeighbor, _num_sample_u)]);
			_frame_normal[iU] += (frame_normal[(iU + iNeighbor) % _num_sample_u]);
		}
		_frame_normal[iU] /= (2.0 * _neighbor + 1.0);

		// project to tangent plane
		_frame_normal[iU] = _frame_normal[iU] - _frame_normal[iU].dot(_rail_tangent[iU]) * _rail_tangent[iU];
	}

	// get binormals
	for (int iU = 0; iU < _num_sample_u; ++iU)
	{
		_frame_binormal[iU] = _rail_tangent[iU].cross(_frame_normal[iU]);
	}
}

const Rail& CushionSurface::getRail() const
{
	return rail_;
}

const Rail& CushionSurface::getOriginalRail() const
{
	return original_rail_;

}

double CushionSurface::getConnectorPullBackLength() const
{
	return connector_pull_back_length_;
}

int CushionSurface::getHalf1Degree() const
{
	int degree = 0;
	if (key_cross_sections_.size() > 0)
	{
		degree = key_cross_sections_[0].cross_section_.getHalf1Degree();
	}
	else
	{
		cout << "WARNING!!!!!!!!!!!!!!CushionSurface::getHalf1Degree() no key cross section" << endl;
	}
	return degree;
}

int CushionSurface::getHalf2Degree() const
{
	int degree = 0;
	if (key_cross_sections_.size() > 0)
	{
		degree = key_cross_sections_[0].cross_section_.getHalf2Degree();
	}
	else
	{
		cout << "WARNING!!!!!!!!!!!!!!CushionSurface::getHalf2Degree() no key cross section" << endl;
	}
	return degree;
}

void CushionSurface::FindSurfaceNormal(std::vector<std::vector<Eigen::Vector3d>>& _node_normals) const
{

	int num_u = discrete_cushion_.size();
	assert(num_u > 0);
	int num_v = discrete_cushion_[0].discrete_cross_section_.size();
	_node_normals.resize(num_u, std::vector<Eigen::Vector3d>(num_v));
	int num_v_interval = num_v - 1;
	int half_num_v = num_v_interval / 2;

	for (int iU = 0; iU < num_u; ++iU)
	{
		// at u
		// get the 2d cross section bezier curve
		// then get the 2d derivative curve
		double u = discrete_cushion_[iU].u_;
		BezierCurve half1;
		BezierCurve half2;
		cross_section_surface_.getBezierCurve_Half1(u, half1);
		cross_section_surface_.getBezierCurve_Half2(u, half2);
		BezierCurve derivative_half1;
		BezierCurve derivative_half2;
		half1.getDerivativeCurve(derivative_half1);
		half2.getDerivativeCurve(derivative_half2);

		// get the derivative value at v
		// 
		// half1
		for (int iV = 0; iV < half_num_v; ++iV)
		{
			double v = (double)iV / (double)num_v_interval * 2;// by default, v in [0,2]
			assert(v < 1.0);

			Eigen::VectorXd derivative;// 2d vector, in the local frame, need to transfer to global
			derivative_half1.getValue(v, derivative);

			// convert to global
			// also rotate derivative 90 degrees to get the normal: (x,y) -> (-y,x)
			_node_normals[iU][iV] = discrete_cushion_[iU].new_frame_normal_ * (- derivative(1)) + discrete_cushion_[iU].new_frame_binormal_ * derivative(0);
			_node_normals[iU][iV].normalize();
		}
		// half2
		for (int iV = half_num_v; iV <= num_v_interval; ++iV)
		{
			double v = (double)iV / (double)num_v_interval * 2;// by default, v in [0,2]
			assert(v >= 1.0);

			Eigen::VectorXd derivative;// 2d vector, in the local frame, need to transfer to global
			derivative_half2.getValue(v - 1.0, derivative);

			// convert to global
			// also rotate derivative 90 degrees to get the normal: (x,y) -> (-y,x)
			_node_normals[iU][iV] = discrete_cushion_[iU].new_frame_normal_ * (-derivative(1)) + discrete_cushion_[iU].new_frame_binormal_ * derivative(0);
			_node_normals[iU][iV].normalize();
		}
	}


}

void CushionSurface::WriteSampleCurves(std::string _filepath, std::string _filename)
{
	// cushion surface obj
	OpenMesh::IO::write_mesh(surface_, _filepath + _filename + "_cushion_surface.obj");

	// rail
	std::ofstream write_rail;
	write_rail.open(_filepath + _filename + "_rail.txt");
	if (write_rail.is_open())
	{
		int rail_sample_num = 700;
		write_rail << rail_sample_num << endl;
		for (int iSample = 0; iSample <= rail_sample_num; ++iSample)// notice we actually sample rail_sample_num+1 points, which should be like this to match CurveToMesh program
		{
			const std::array<double, 2>& interval = rail_.getDomain();
			float para = (float)iSample / rail_sample_num;
			para *= (interval[1] - interval[0]);
			para += interval[0];
			Eigen::VectorXd value;
			rail_.getValue(para, value);
			write_rail << value(0) << " " << value(1) << " " << value(2) << endl;
		}
		write_rail.close();
	}

	FindKeyCrossSectionFrameId();
	// key cross sections
	if (1)
	{
		for (int iKey = 0; iKey < key_cross_sections_.size(); ++iKey)
		{
			std::ofstream write_key_cross_section;
			write_key_cross_section.open(_filepath + _filename + "_KeyCrossSection" + std::to_string(iKey) + ".txt");
			if (write_key_cross_section.is_open())
			{
				int frame_id = key_cross_sections_[iKey].frame_id_;
				write_key_cross_section << discrete_cushion_[frame_id].discrete_cross_section_.size() * 2 << endl;
				// go from begining to end, then back to begining in a loop
				for (int iV = 0; iV < discrete_cushion_[frame_id].discrete_cross_section_.size(); ++iV)
				{
					write_key_cross_section << discrete_cushion_[frame_id].discrete_cross_section_[iV](0) << " " <<
						discrete_cushion_[frame_id].discrete_cross_section_[iV](1) << " " <<
						discrete_cushion_[frame_id].discrete_cross_section_[iV](2) << endl;
				}
				for (int iV = discrete_cushion_[frame_id].discrete_cross_section_.size() - 1; iV >= 0; --iV)
				{
					write_key_cross_section << discrete_cushion_[frame_id].discrete_cross_section_[iV](0) << " " <<
						discrete_cushion_[frame_id].discrete_cross_section_[iV](1) << " " <<
						discrete_cushion_[frame_id].discrete_cross_section_[iV](2) << endl;
				}

				write_key_cross_section.close();
			}
		}
	}
	
	Mesh axis_x_mesh;
	Mesh axis_y_mesh;
	Mesh axis_z_mesh;
	Mesh origin_mesh;
	OpenMesh::IO::read_mesh(axis_x_mesh, "C:\\Users\\15539\\Desktop\\MaskDesign\\Data\\Render\\axis_x.obj");
	OpenMesh::IO::read_mesh(axis_y_mesh, "C:\\Users\\15539\\Desktop\\MaskDesign\\Data\\Render\\axis_y.obj");
	OpenMesh::IO::read_mesh(axis_z_mesh, "C:\\Users\\15539\\Desktop\\MaskDesign\\Data\\Render\\axis_z.obj");
	OpenMesh::IO::read_mesh(origin_mesh, "C:\\Users\\15539\\Desktop\\MaskDesign\\Data\\Render\\origin0.5.obj");
	// key frames
	if (0)
	{
		for (int iKey = 0; iKey < key_cross_sections_.size(); ++iKey)
		{
			int frame_id = key_cross_sections_[iKey].frame_id_;
			Mesh axis_x_mesh_transformed = axis_x_mesh;
			Mesh axis_y_mesh_transformed = axis_y_mesh;
			Mesh axis_z_mesh_transformed = axis_z_mesh;
			Mesh origin_mesh_transformed = origin_mesh;


			Eigen::Matrix3d basis(Eigen::Matrix3d::Identity());
			basis.block(0, 0, 3, 1) = discrete_cushion_[frame_id].rotated_tangent_;
			basis.block(0, 1, 3, 1) = discrete_cushion_[frame_id].new_frame_normal_;
			basis.block(0, 2, 3, 1) = discrete_cushion_[frame_id].new_frame_binormal_;
			TransformMesh transform_it0(axis_x_mesh_transformed);
			transform_it0.RotateTranslate(basis, discrete_cushion_[frame_id].frame_value_);
			TransformMesh transform_it1(axis_y_mesh_transformed);
			transform_it1.RotateTranslate(basis, discrete_cushion_[frame_id].frame_value_);
			TransformMesh transform_it2(axis_z_mesh_transformed);
			transform_it2.RotateTranslate(basis, discrete_cushion_[frame_id].frame_value_);
			TransformMesh transform_it3(origin_mesh_transformed);
			transform_it3.RotateTranslate(basis, discrete_cushion_[frame_id].frame_value_);


			OpenMesh::IO::write_mesh(axis_x_mesh_transformed, _filepath + _filename + "_KeyFrameX" + std::to_string(iKey) + ".obj");
			OpenMesh::IO::write_mesh(axis_y_mesh_transformed, _filepath + _filename + "_KeyFrameY" + std::to_string(iKey) + ".obj");
			OpenMesh::IO::write_mesh(axis_z_mesh_transformed, _filepath + _filename + "_KeyFrameZ" + std::to_string(iKey) + ".obj");
			OpenMesh::IO::write_mesh(origin_mesh_transformed, _filepath + _filename + "_KeyOrigin" + std::to_string(iKey) + ".obj");
		}
	}

	// control points and control polygon
	if (0)
	{
		// control points
		Mesh cushion_rail_control_mesh;//control points
		rail_.getControlPointMesh(cushion_rail_control_mesh);
		for (int iPoint = 0; iPoint < cushion_rail_control_mesh.n_vertices() - 1; ++iPoint)
		{
			Mesh origin_mesh_transformed = origin_mesh;
			TransformMesh transform_it3(origin_mesh_transformed);
			transform_it3.RotateTranslate(Eigen::Matrix3d::Identity(),
				Eigen::Map<Eigen::Vector3d>(cushion_rail_control_mesh.point(cushion_rail_control_mesh.vertex_handle(iPoint)).data()));
			OpenMesh::IO::write_mesh(origin_mesh_transformed, _filepath + _filename + "_ControlPoints" + std::to_string(iPoint) + ".obj");

		}
		// control polygon
		std::ofstream write_control_polygon;
		write_control_polygon.open(_filepath + _filename + "_ControlPolygon.txt");
		if (write_control_polygon.is_open())
		{
			int sample_num_each_segment = 70;
			int total_sample_num = sample_num_each_segment * (cushion_rail_control_mesh.n_vertices() - 1);
			write_control_polygon << total_sample_num - 1 << endl;
			for (int iPoint = 0; iPoint < cushion_rail_control_mesh.n_vertices() - 1; ++iPoint)
			{
				Mesh::Point current = cushion_rail_control_mesh.point(cushion_rail_control_mesh.vertex_handle(iPoint));
				Mesh::Point next = cushion_rail_control_mesh.point(cushion_rail_control_mesh.vertex_handle(iPoint + 1));

				for (int iSample = 0; iSample < sample_num_each_segment; ++iSample)
				{
					float para = (float)iSample / (sample_num_each_segment - 1);
					Mesh::Point sample_point = (1.0 - para) * current + para * next;
					write_control_polygon << sample_point[0] << " " << sample_point[1] << " " << sample_point[2] << endl;

				}
			}
			write_control_polygon.close();
		}
	}

	// more frame normal samples
	if (0)
	{
		
		int sample_num = 100;
		computeSampleFrames(20, sample_num);
		for (int iFrame = 0; iFrame < sample_num; ++iFrame)
		{
			Mesh axis_y_mesh_transformed = axis_y_mesh;
			Mesh origin_mesh_transformed = origin_mesh;

			Eigen::Matrix3d basis(Eigen::Matrix3d::Identity());
			basis.block(0, 0, 3, 1) = discrete_cushion_[iFrame].rotated_tangent_;
			basis.block(0, 1, 3, 1) = discrete_cushion_[iFrame].new_frame_normal_;
			basis.block(0, 2, 3, 1) = discrete_cushion_[iFrame].new_frame_binormal_;
			TransformMesh transform_it1(axis_y_mesh_transformed);
			transform_it1.RotateTranslate(basis, discrete_cushion_[iFrame].frame_value_);
			TransformMesh transform_it3(origin_mesh_transformed);
			transform_it3.RotateTranslate(basis, discrete_cushion_[iFrame].frame_value_);
			OpenMesh::IO::write_mesh(axis_y_mesh_transformed, _filepath + _filename + "_FrameSamplesNormal" + std::to_string(iFrame) + ".obj");
			OpenMesh::IO::write_mesh(origin_mesh_transformed, _filepath + _filename + "_FrameSamplesOrigin" + std::to_string(iFrame) + ".obj");

		}
	}
}

void CushionSurface::SaveResults(std::string _filepath, std::string _filename) const
{
	// rail curve
	std::string rail_file = _filepath + _filename + "_rail.txt";
	rail_.SaveResults(rail_file);

	// key cross sections control points
	std::vector<std::string> key_cross_section_files(key_cross_sections_.size());
	std::vector<double> key_cross_section_u(key_cross_sections_.size());
	for (int iKey = 0; iKey < key_cross_sections_.size(); ++iKey)
	{
		key_cross_section_files[iKey] = _filepath + _filename + "_KeyCrossSection" + std::to_string(iKey) + ".txt";
		key_cross_sections_[iKey].cross_section_.saveToFile(key_cross_section_files[iKey]);
		key_cross_section_u[iKey] = key_cross_sections_[iKey].u_;
	}

	int sample_u = 300;
	int sample_v = 50;
	if (discrete_cushion_.size() > 0)
	{
		sample_u = discrete_cushion_.size();
		if (discrete_cushion_[0].discrete_cross_section_.size() > 1)
		{
			sample_v = discrete_cushion_[0].discrete_cross_section_.size() - 1;
		}
	}

	//// arrange data in xml file
	//XmlParser save_data;
	//save_data.saveCushionSurface(_filepath + _filename + ".xml", rail_file, key_cross_section_files, key_cross_section_u, 
	//	sample_u, sample_v, cross_section_control_curve_degree_, connector_pull_back_length_);
}

void CushionSurface::ReadResults(std::string _xml_file, int _generalized_cushion)
{
	// read data
	std::string rail_file;
	std::vector<std::string> key_cross_section_files;
	std::vector<double> key_cross_section_u;
	int sample_u = 0;
	int sample_v = 0;
	//XmlParser read_data;
	//read_data.readCushionSurface(_xml_file, rail_file, key_cross_section_files, key_cross_section_u,
	//	sample_u, sample_v, cross_section_control_curve_degree_, connector_pull_back_length_);
	// construct cushion with the data
	rail_.ReadResults(rail_file);

	key_cross_sections_.resize(key_cross_section_files.size());
	assert(key_cross_section_u.size() == key_cross_section_files.size());
	for (int iKey = 0; iKey < key_cross_section_files.size(); ++iKey)
	{
		key_cross_sections_[iKey].cross_section_.loadFromFile(key_cross_section_files[iKey]);
		key_cross_sections_[iKey].u_ = key_cross_section_u[iKey];
	}

	ConstructCrossSectionSurface(cross_section_control_curve_degree_);

	double frame_relaxation_radius = 20.0;
	if (_generalized_cushion == 1)
	{
		frame_relaxation_radius = 1.0;
	}
	computeSampleFrames(frame_relaxation_radius, sample_u);

	//// TODO read the connector curve data 
	//FilePath get_path;
	//std::string data_path = get_path.getExePath_ParentFolder2() + "\\data";
	//setFixedBoundary(data_path + "\\outside_boundary.txt", 101, true);
	

	computeDiscreteCrossSections_NewFrames(sample_v);
	constructSurfaceMesh();

}

void CushionSurface::InitializeAsGenericCushion(std::string _data_file)
{
	// data
	int sample_u = 300;
	int sample_v = 50;
	std::string rail_file = _data_file + "/GenericCushion_rail.txt";
	std::vector<std::string> key_cross_section_files(17);
	for (int iKey = 0; iKey < key_cross_section_files.size(); ++iKey)
	{
		key_cross_section_files[iKey] = _data_file + "/GenericCushion_KeyCrossSection" + std::to_string(iKey) + ".txt";
	}
	std::vector<double> key_cross_section_u{ 0.0266667,0.0633333,0.12,0.18,0.23,0.32,0.403333,0.463333,0.506667,0.536667,0.576667,0.663333,0.723333,0.803333,0.866667,0.92,0.98 };
	cross_section_control_curve_degree_ = 3;
	connector_pull_back_length_ = 0;
	Mesh extracted_centerline_curve;
	OpenMesh::IO::read_mesh(extracted_centerline_curve, _data_file + "/extracted_rail_polygon.obj");
	std::vector<Eigen::VectorXd> extracted_polygon;
	SubMesh get_polygon;
	get_polygon.PointCloudMeshToEigenVectorX(extracted_centerline_curve, extracted_polygon);

	// construct cushion with the data

	setCushionPolyline(extracted_polygon);
	rail_.ReadResults(rail_file);
	setOriginalRail(rail_);
	double frame_relaxation_radius = 20.0;
	computeSampleFrames(frame_relaxation_radius, sample_u);

	key_cross_sections_.resize(key_cross_section_files.size());
	assert(key_cross_section_u.size() == key_cross_section_files.size());
	for (int iKey = 0; iKey < key_cross_section_files.size(); ++iKey)
	{
		key_cross_sections_[iKey].cross_section_.loadFromFile(key_cross_section_files[iKey]);
		key_cross_sections_[iKey].u_ = key_cross_section_u[iKey];
	}
	ConstructCrossSectionSurface(cross_section_control_curve_degree_);

	setFixedBoundary(_data_file + "/outside_boundary.txt", 100, true);
	computeDiscreteCrossSections_NewFrames(sample_v);
	constructSurfaceMesh();
	setCurrentMeshAsOriginalCushionMesh();
}

void CushionSurface::constructStandardSurfaceWireFrame(int _num_sample_u, int _num_sample_v)
{
	computeStandardFrames(_num_sample_u);
	computeDiscreteCrossSections_StandardFrames(_num_sample_v);
}

//void CushionSurface::constructRelaxedSurfaceWireFrame(float _radius, int _num_sample_u, int _num_sample_v)
//{
//	computeFrames(_radius, _num_sample_u);
//	computeDiscreteCrossSections_NewFrames(_num_sample_v);
//}


void CushionSurface::ConstructSelfIntersectPartMesh(std::vector<Mesh>& _part_meshes) const
{
	std::vector<std::pair<int, int>> intersect_part;
	CheckIntersection(discrete_cushion_, intersect_part);

	_part_meshes.resize(intersect_part.size());
	for (int iPart = 0; iPart < intersect_part.size(); ++iPart)
	{
		int begin_idx = intersect_part[iPart].first;
		int end_idx = intersect_part[iPart].second;
		int total_num_frame = discrete_cushion_.size();
		if (end_idx < begin_idx)//periodic case
		{
			end_idx += total_num_frame;
		}
		assert(begin_idx < total_num_frame);

		int num_frames = end_idx - begin_idx + 1;
		std::vector<std::vector<Mesh::Point>> uv_points;
		uv_points.reserve(num_frames);
		for (int iFrame = 0; iFrame < num_frames; ++iFrame)
		{
			std::vector<Eigen::Vector3d> v_points_eigen = discrete_cushion_[(begin_idx + iFrame) % total_num_frame].discrete_cross_section_;
			std::vector<Mesh::Point> v_points;//convert to EigenVector to Mesh::Point
			v_points.reserve(v_points_eigen.size());
			for (auto& point : v_points_eigen)
			{
				v_points.push_back(Mesh::Point(point(0), point(1), point(2)));
			}
			uv_points.push_back(v_points);
		}

		SubMesh construct_mesh;
		construct_mesh.MeshFromUVSurface(uv_points, _part_meshes[iPart]);
	}


}

void CushionSurface::FrameRelaxation(int _interval_side_expand, bool _extra_iteration)
{
	std::vector<std::pair<int, int>> intersect_part;
	CheckIntersection(discrete_cushion_, intersect_part);

	FrameRelax relax_it;
	relax_it.Relax(discrete_cushion_, intersect_part, _interval_side_expand, _extra_iteration);


}

void CushionSurface::constructVUPoints(std::vector<std::vector<Mesh::Point>>& _vu_points) const
{

	//constructSurfaceWireFrame(_num_sample_u, _num_sample_v);
	int num_sample_u = discrete_cushion_.size();
	if (num_sample_u > 0)
	{
		int num_sample_v = discrete_cushion_[0].discrete_cross_section_.size();

		_vu_points.resize(num_sample_v);
		for (int iV = 0; iV < num_sample_v; ++iV)
		{
			std::vector<Mesh::Point> u_points(num_sample_u);
			for (int iU = 0; iU < num_sample_u; ++iU)
			{
				Eigen::Vector3d point_eigen = discrete_cushion_[iU].discrete_cross_section_[iV];
				u_points[iU] = Mesh::Point(point_eigen(0), point_eigen(1), point_eigen(2));
			}
			_vu_points[iV] = u_points;
		}
	}
}

void CushionSurface::constructUVPoints(std::vector<std::vector<Mesh::Point>>& _uv_points) const
{
	int num_sample_u = discrete_cushion_.size();
	if (num_sample_u > 0)
	{
		int num_sample_v = discrete_cushion_[0].discrete_cross_section_.size();

		_uv_points.resize(num_sample_u);
		for (int iU = 0; iU < num_sample_u; ++iU)
		{
			std::vector<Mesh::Point> v_points(num_sample_v);
			for (int iV = 0; iV < num_sample_v; ++iV)
			{
				Eigen::Vector3d point_eigen = discrete_cushion_[iU].discrete_cross_section_[iV];
				v_points[iV] = Mesh::Point(point_eigen(0), point_eigen(1), point_eigen(2));
			}
			_uv_points[iU]= v_points;
		}
	}
}

void CushionSurface::constructSurfaceMesh(/*int _num_sample_u, int _num_sample_v*/)
{
	surface_.clear();

	std::vector<std::vector<Mesh::Point>> uv_points;
	constructUVPoints(uv_points);

	SubMesh construct_mesh;
	construct_mesh.MeshFromUVSurface_PeriodU(uv_points, surface_);
}

void CushionSurface::constructSmoothedSurfaceMesh(int _num_sample_u, int _num_sample_v, int _smoothness)
{
	surface_.clear();

	std::vector < Eigen::Vector3d> rail_value;
	std::vector < Eigen::Vector3d> rail_tangent;
	std::vector < Eigen::Vector3d> frame_normal;
	std::vector < Eigen::Vector3d> frame_binormal;
	getRailFrameWithSmoothedFaceNormals(_num_sample_u, rail_value, rail_tangent, frame_normal, frame_binormal, _smoothness);

	const std::array<double, 2>& domain_u = rail_.getDomain();
	std::vector<std::vector<Mesh::Point>> uv_points;
	uv_points.reserve(_num_sample_u + 1);// +1 for forming a closed surface
	std::vector<Mesh::Point> v_points(_num_sample_v);

	for (int iU = 0; iU < _num_sample_u; ++iU)
	{
		double u = (double)iU / (double)_num_sample_u * (domain_u[1] - domain_u[0]) + domain_u[0];
		for (int iV = 0; iV < _num_sample_v; ++iV)
		{
			double v = (double)iV / (double)_num_sample_v * 2;// by default, v in [0,2]

			Eigen::Vector3d cross_section_value;
			cross_section_surface_.getValue_Degree2(u, v, cross_section_value);

			// surface formula
			Eigen::Vector3d value =
				rail_value[iU] +
				rail_tangent[iU] * cross_section_value[0] +
				frame_normal[iU] * cross_section_value[1] +
				frame_binormal[iU] * cross_section_value[2];
			v_points[iV][0] = value[0];
			v_points[iV][1] = value[1];
			v_points[iV][2] = value[2];
		}
		uv_points.push_back(v_points);
	}
	SubMesh construct_mesh;
	construct_mesh.MeshFromUVSurface_PeriodU(uv_points, surface_);
}

const Mesh& CushionSurface::getSurfaceMesh() const
{
	return surface_;
}

void CushionSurface::setCurrentMeshAsOriginalCushionMesh()
{
	original_fitted_surface_ = surface_;
}

void CushionSurface::setOriginalRail(const Rail& _rail)
{
	original_rail_ = _rail;
}

void CushionSurface::setOriginalCushionMesh(const Mesh& _cushion_surface)
{
	original_fitted_surface_ = _cushion_surface;
}

const std::vector<DiscreteCushionPiece>* CushionSurface::getDiscreteCushion() const
{
	return &discrete_cushion_;
}

void CushionSurface::FlattenInnerHalf()
{
	for (auto& iKeyCrossSection : key_cross_sections_)
	{
		iKeyCrossSection.cross_section_.FlattenHalf1_ParallelToAxisX();
	}
}

double CushionSurface::ComputeInnerBoundaryLength() const
{
	double length = 0.0;
	for (int iSection = 0; iSection < discrete_cushion_.size() - 1; ++iSection)
	{
		double segment_length = (discrete_cushion_[iSection].discrete_cross_section_[0] - discrete_cushion_[iSection + 1].discrete_cross_section_[0]).norm();
		length += segment_length;
	}
	return length;// unit is mm
}

void CushionSurface::computeStandardFrames(int _num_sample_u)
{
	discrete_cushion_.resize(_num_sample_u);
	const std::array<double, 2>& domain_u = rail_.getDomain();

	for (int iU = 0; iU < _num_sample_u; ++iU)
	{
		discrete_cushion_[iU].u_ = (double)iU / (double)_num_sample_u * (domain_u[1] - domain_u[0]) + domain_u[0];

		getRailFrame_AxisZ(discrete_cushion_[iU].u_, 
			discrete_cushion_[iU].frame_value_, discrete_cushion_[iU].frame_tangent_, discrete_cushion_[iU].frame_normal_);

		discrete_cushion_[iU].frame_binormal_ = discrete_cushion_[iU].frame_tangent_.cross(discrete_cushion_[iU].frame_normal_);

		discrete_cushion_[iU].rotation_ = Eigen::Quaterniond(1.0, 0.0, 0.0, 0.0);
		discrete_cushion_[iU].rotated_tangent_ = discrete_cushion_[iU].frame_tangent_;

	}
}

void CushionSurface::computeNewFrames_AxisZ()
{
	for (auto& iPiece : discrete_cushion_)
	{
		getRailNormal_AxisZ(iPiece.rotated_tangent_, iPiece.new_frame_normal_);
		iPiece.new_frame_binormal_ = iPiece.rotated_tangent_.cross(iPiece.new_frame_normal_);
	}
}


void CushionSurface::computeDiscreteCrossSections_StandardFrames(int _num_sample_v)
{
	for (auto& iFrame : discrete_cushion_)
	{
		iFrame.discrete_cross_section_.resize(_num_sample_v);

		for (int iV = 0; iV < _num_sample_v; ++iV)
		{
			double v = (double)iV / (double)_num_sample_v * 2;// by default, v in [0,2]

			getValue(iFrame.u_, v, iFrame.frame_value_, iFrame.frame_tangent_, iFrame.frame_normal_, iFrame.frame_binormal_,
				iFrame.discrete_cross_section_[iV]);
		}
	}

}

void CushionSurface::UpdateFrameValue()
{
	for (auto& iCrossSection : discrete_cushion_)
	{
		Eigen::VectorXd current_rail_value;
		rail_.getValue(iCrossSection.u_, current_rail_value);

		iCrossSection.frame_value_ = current_rail_value;
	}
}

void CushionSurface::computeDiscreteCrossSections_NewFrames(int _num_sample_v)
{
	// update cross section surface values
	u_values_.resize(discrete_cushion_.size());
	v_values_.resize(_num_sample_v + 1);
	for (int iU = 0; iU < discrete_cushion_.size(); ++iU)
	{
		u_values_[iU] = discrete_cushion_[iU].u_;
	}

	for (int iV = 0; iV <= _num_sample_v; ++iV)
	{
		v_values_[iV] = (double)iV / (double)_num_sample_v * 2;// by default, v in [0,2]
	}

	cross_section_surface_.ComputeSampledValues(u_values_, v_values_);

	// update discrete_cushion_
	for (int iU = 0; iU < discrete_cushion_.size(); ++iU)
	{
		discrete_cushion_[iU].discrete_cross_section_.resize(_num_sample_v + 1);

		for (int iV = 0; iV <= _num_sample_v; ++iV)
		{
			Eigen::Vector2d cross_section_value;

			cross_section_surface_.getValue(iU, iV, cross_section_value);

			discrete_cushion_[iU].discrete_cross_section_[iV] = 	// surface formula
				discrete_cushion_[iU].frame_value_ +
				discrete_cushion_[iU].new_frame_normal_ * cross_section_value[0] +
				discrete_cushion_[iU].new_frame_binormal_ * cross_section_value[1];

		}
	}
}

void CushionSurface::FindCrossSectionForAllFrames(std::vector<CrossSection>& _cross_sections) const
{
	_cross_sections.resize(discrete_cushion_.size());
	for (int iFrame = 0; iFrame < discrete_cushion_.size(); ++iFrame)
	{
		FindCrossSectionAtOneFrame(_cross_sections[iFrame], iFrame);
	}
}

void CushionSurface::CheckIntersection(const std::vector<DiscreteCushionPiece>& _cross_sections, std::vector<std::pair<int, int>>& _intersect_parts) const
{
	SelfIntersectionTest find_intersection;
	find_intersection.TestIntersection(_cross_sections, _intersect_parts);
}

void CushionSurface::ClosestFrameIdxToU(const std::vector<double>& _u, std::vector<int>& _i_frame) const
{
	_i_frame.clear();
	_i_frame.reserve(_u.size());
	int iFrame = 0;
	for (int iU = 0; iU < _u.size(); ++iU)
	{
		int counter = 0;
		while (iFrame < discrete_cushion_.size() - 1 && discrete_cushion_[iFrame].u_ < _u[iU])
		{
			++iFrame;
			++counter;
		}
		if (counter > 0)// prevent duplicate frame index
		{
			_i_frame.push_back(iFrame);
		}
	}
}

double CushionSurface::getLocalLinearParameter(double _u, double _u0, double _u1) const
{
	return (_u - _u0) / (_u1 - _u0);
}

void CushionSurface::FindCrossSectionPolyline(const Mesh& _cushion, const Mesh::Point& _frame_tangent, const Mesh::Point& _frame_value, const Mesh::Point& _indicator, 
	std::vector<Mesh::Point>& _result_polyline) const
{
	// intersect the rail normal plane with the cushion
	MeshPlaneIntersection find_section_polyline;
	find_section_polyline.setData(&_cushion, _frame_tangent, _frame_value);
	find_section_polyline.FindIntersectPolyline_AllParts();
	find_section_polyline.getIntersectPoints_ClosetToAPoint(_indicator, _result_polyline);
}

int CushionSurface::AddCenterToCrossSectionPolyline(std::vector<Mesh::Point>& _cross_section, const Mesh::Point& _center) const
{
	// find the two adjacent vertex in polyline that's closet to the center point
	double min_distance = FLT_MAX;
	std::vector<Mesh::Point>::const_iterator min_iter = _cross_section.cbegin();
	if (_cross_section.size() > 0)
	{
		std::vector<Mesh::Point>::const_iterator iter_end = _cross_section.end() - 1;
		for (auto iterVertex = _cross_section.cbegin(); iterVertex != iter_end; ++iterVertex)
		{
			double distance = (*iterVertex - _center).sqrnorm() + (*(iterVertex + 1) - _center).sqrnorm();
			if (distance < min_distance)
			{
				min_distance = distance;
				min_iter = iterVertex;
			}
		}
	}
	//cout << "point at min_iter + 1\n" << *(min_iter + 1) << endl;
	//cout << "std::distance(_cross_section.cbegin(), min_iter + 1)=" << std::distance(_cross_section.cbegin(), min_iter + 1) << endl;
	std::vector<Mesh::Point>::const_iterator center_iter = _cross_section.insert(min_iter + 1, _center);
	//cout << "std::distance(_cross_section.cbegin(), center_iter)\n" << std::distance(_cross_section.cbegin(), center_iter) << endl;
	//cout << "point at center_iter\n" << *center_iter << endl;
	//cout << "while center=\n" << _center << endl;

	return std::distance(_cross_section.cbegin(), center_iter);
}

void CushionSurface::getNewFrames(int _dicrete_cross_section_idx, Eigen::Matrix4d& _new_frame) const
{
	_new_frame = Eigen::Matrix4d::Identity();
	_new_frame.block(0, 0, 3, 1) = discrete_cushion_[_dicrete_cross_section_idx].rotated_tangent_;
	_new_frame.block(0, 1, 3, 1) = discrete_cushion_[_dicrete_cross_section_idx].new_frame_normal_;
	_new_frame.block(0, 2, 3, 1) = discrete_cushion_[_dicrete_cross_section_idx].new_frame_binormal_;
	_new_frame.block(0, 3, 3, 1) = discrete_cushion_[_dicrete_cross_section_idx].frame_value_;
}

void CushionSurface::AdjustCrossSectionPolylineToPointingOutward(std::vector<Mesh::Point>& _cross_section) const
{
	if (_cross_section.empty() == false)
	{
		//Mesh::Point direction = _cross_section.back() - _cross_section[0];
		//// we regard cases:
		//// [direction.dot(_cross_section[0]) < 0 && direction.dot(_cross_section.back()) < 0] equivalent to inward order
		//// [direction.dot(_cross_section[0]) > 0 && direction.dot(_cross_section.back()) > 0] equivalent to outward order
		//// otherwise if [direction.dot(_cross_section[0]) * direction.dot(_cross_section.back())] < 0, is an unexpected abnormal case, which this trick cannot handle
		//if (direction.dot(_cross_section[0]) < 0 && direction.dot(_cross_section.back()) < 0)
		//{
		//	// inward order
		//	// reverse the polyline
		//	std::reverse(_cross_section.begin(), _cross_section.end());
		//}
		//else if ((direction.dot(_cross_section[0]) * direction.dot(_cross_section.back())) <= 0)
		//{
		//	cout << "[WARNING FROM CushionSurface::AdjustCrossSectionPolylineToPointingOutward] unable to detect correct order of the cross section polyline" << endl;
		//}

		if (_cross_section.back().norm() < _cross_section[0].norm())
		{
			// inward order
			// reverse the polyline
			std::reverse(_cross_section.begin(), _cross_section.end());
		}

	}
}

void CushionSurface::ProjectPolylineToFace(const std::vector<Eigen::VectorXd>& _polyline, std::vector<Eigen::VectorXd>& _projected_polyline) const
{
	int num_data = _polyline.size();
	_projected_polyline.clear();
	_projected_polyline.reserve(num_data);
	for (int iPoint = 0; iPoint < num_data; ++iPoint)
	{
		Eigen::VectorXd point_eigen = _polyline[iPoint];
		Mesh::Point point(point_eigen(0), point_eigen(1), point_eigen(2));
		Mesh::Point projected;
		uint face_id = 0; // no use, only for input of the function
		bool is_intersect = head_tree_->LineMeshIntersection(point, Mesh::Point(0.0, 0.0, -1.0), projected, face_id);
		if (!is_intersect)
		{
			cout << "//WARNING// polyline projection fail, possibly bad mask face alignment" << endl;
		}
		_projected_polyline.push_back(Eigen::Map<Eigen::VectorXd>(projected.data(), 3));
	}
}

void CushionSurface::getOriginalCrossSection(int _frame, std::vector<Mesh::Point>& _cross_section_polyline) const
{
	// get cross section polylines by cutting it with frame normal plane

	Mesh::Point frame_tangent(discrete_cushion_[_frame].rotated_tangent_(0),
		discrete_cushion_[_frame].rotated_tangent_(1),
		discrete_cushion_[_frame].rotated_tangent_(2));
	Mesh::Point frame_value(discrete_cushion_[_frame].frame_value_(0),
		discrete_cushion_[_frame].frame_value_(1),
		discrete_cushion_[_frame].frame_value_(2));

	// intersect the rail normal plane with the cushion
	FindCrossSectionPolyline(*cushion_, frame_tangent, frame_value, frame_value, _cross_section_polyline);
}

void CushionSurface::getOriginalCrossSection_OffAxis(int _frame, const Mesh& _original_cushion, const std::vector<Mesh::Point>& _original_rail_polyline,
	std::vector<Mesh::Point>& _cross_section_polyline, int& _center_index) const
{
	Mesh::Point frame_tangent(discrete_cushion_[_frame].rotated_tangent_(0),
		discrete_cushion_[_frame].rotated_tangent_(1),
		discrete_cushion_[_frame].rotated_tangent_(2));
	Mesh::Point frame_value(discrete_cushion_[_frame].frame_value_(0),
		discrete_cushion_[_frame].frame_value_(1),
		discrete_cushion_[_frame].frame_value_(2));

	// the intersect point of the original rail with the frame normal plane
	// since there might be multiple intersect points, I think we can get the right one by finding the closest one to the frame value
	Mesh::Point correspond_original_rail_point;
	PolylinePlaneIntersection get_correspond_rail_point;
	get_correspond_rail_point.ComputeInterectionPoints_Polygon_ClosestOne(_original_rail_polyline, frame_tangent, frame_value, correspond_original_rail_point);

	// the intersect cross_section polyline of the original cushion with the frame normal plane
	FindCrossSectionPolyline(_original_cushion, frame_tangent, frame_value, correspond_original_rail_point, _cross_section_polyline);
	AdjustCrossSectionPolylineToPointingOutward(_cross_section_polyline);

	_center_index = AddCenterToCrossSectionPolyline(_cross_section_polyline, correspond_original_rail_point);

}

void CushionSurface::FindCrossSectionAtOneFrame(CrossSection& _cross_section, int _frame) const
{
	const double& u = discrete_cushion_[_frame].u_;

	std::vector<Eigen::Vector2d> half1_rest;
	double half1_point1_length;
	double tangent_angle;
	double half2_point1_length;
	std::vector<Eigen::Vector2d> half2_rest;
	cross_section_surface_.getControlParameters(u, half1_rest, half1_point1_length, tangent_angle, half2_point1_length, half2_rest);

	Eigen::Vector2d point0;
	cross_section_surface_.getPoint0(u, point0);

	_cross_section.setControlPoints(half1_rest, half1_point1_length, tangent_angle, half2_point1_length, half2_rest);
	_cross_section.setPoint0(point0);
}

void CushionSurface::FindKeyCrossSectionFrameId()
{
	int num_key = key_cross_sections_.size();
	// find frame idx from the u value
	std::vector<double> u(num_key);
	for (int iKey = 0; iKey < num_key; ++iKey)
	{
		u[iKey] = key_cross_sections_[iKey].u_;
	}
	std::vector<int> frame;
	ClosestFrameIdxToU(u, frame);

	for (int iKey = 0; iKey < num_key; ++iKey)
	{
		key_cross_sections_[iKey].frame_id_ = frame[iKey];
	}
}

void CushionSurface::setFixedBoundary(const std::string& _original_boundary, int _num_interpolate_point, bool _do_connector_pull_back)
{
	// read original boundary from a Rhino modelled file
	std::vector<Mesh::Point> outside_boundary;
	ReadRhino read_boundary;
	read_boundary.Read3DPoints_Brackets(_original_boundary, outside_boundary, 1);
	Mesh original_boundary;
	SubMesh convert_typer;
	convert_typer.PointCloudMeshFromPoints(outside_boundary, original_boundary);

	setFixedBoundary(original_boundary, _num_interpolate_point, _do_connector_pull_back);
}

void CushionSurface::setFixedBoundary(const Mesh& _boundary, int _num_interpolate_point, bool _do_connector_pull_back)
{
	Mesh pulled_back_boundary = _boundary;
	if (_do_connector_pull_back)
	{
		// pull back connector according to target translation length
		if (connector_pull_back_length_ > 0.0)
		{
			PullBackConnector* pull_it_back = new PullBackConnector;
			pull_it_back->TranslateMeshAlongZ(pulled_back_boundary, connector_pull_back_length_);
			delete pull_it_back;
		}
	}

	// NOTE this sampling not creating exactly _num_interpolate_point of points......
	// choose frames to interpolate, according to designated interpolation point number
	//int interval_size = discrete_cushion_.size() / _num_interpolate_point;
	//if (discrete_cushion_.size() <= _num_interpolate_point)
	//{
	//	interval_size = 1;
	//}

	// convert boundary data type to vector
	std::vector<Mesh::Point> pulled_back_boundary_points;
	SubMesh convert_type;
	convert_type.PointCloudMeshToPoints(pulled_back_boundary, pulled_back_boundary_points);

	std::vector<int> sample_frames;
	std::vector<Mesh::Point> intersect_points;
	FindConnectorAndSampleFrameIntersectPoints(pulled_back_boundary_points, _num_interpolate_point, sample_frames, intersect_points);


	// set boundary

	// get the interpolation points and corresponding u values

	std::vector<Eigen::VectorXd> interpolate_points;
	std::vector<double> u;
	interpolate_points.reserve(sample_frames.size());
	u.reserve(sample_frames.size());
	for (int iSample = 0; iSample < sample_frames.size(); ++iSample)
	{
		u.push_back(discrete_cushion_[sample_frames[iSample]].u_);

		//
		// 1. compute intersection point of frame normal plane with boundary polyline
		// 
		//Mesh::Point intersect_point;
		//Mesh::Point plane_normal(discrete_cushion_[iSample].rotated_tangent_(0),
		//						discrete_cushion_[iSample].rotated_tangent_(1),
		//						discrete_cushion_[iSample].rotated_tangent_(2));
		//Mesh::Point rail_point(discrete_cushion_[iSample].frame_value_(0),
		//						discrete_cushion_[iSample].frame_value_(1),
		//						discrete_cushion_[iSample].frame_value_(2));
		//PolylinePlaneIntersection find_point;
		//find_point.ComputeInterectionPoints_Polygon_ClosestOne(original_boundary, plane_normal, rail_point, intersect_point);



		// 
		// 2. get the point in 2d local frame : similar to (CrossSection::compute2DLocalCoordinateOf3DPolyline)
		//
		// find transformation matrix
		Eigen::Matrix4d local_frame;
		getNewFrames(sample_frames[iSample], local_frame);
		Eigen::Matrix4d world_frame(Eigen::Matrix4d::Identity());
		Eigen::Matrix4d transform_to_local;// the transformation matrix
		ChangeFrame change_frame(world_frame, local_frame, transform_to_local);

		// transform
		Eigen::Vector4d point_world(intersect_points[iSample][0], intersect_points[iSample][1], intersect_points[iSample][2], 1.0);
		Eigen::Vector4d point_local = transform_to_local * point_world;
		interpolate_points.push_back(Eigen::Vector2d(point_local[1], point_local[2]));

	}

	//interpolate all the points : similar to (CrossSectionSurface::setControlCurveByInterpolation)
	//
	// check if first u is 0, if it is, we should remove it for interpolation algorithm not support this...
	if (u.size() > 0)
	{
		if (u[0] < 1e-10)
		{
			u.erase(u.begin());
			interpolate_points.erase(interpolate_points.begin());
		}
	}
	if (u.size() == 0)
	{
		cout << "ERROR!! CushionSurface::setFixedBoundary interpolation u size =0" << endl;
	}
	cross_section_surface_.setBoundaryControlCurveByInterpolation(interpolate_points, u);

	// at last, we just the key cross section control points accordingly
	for (auto& iCrossSection : key_cross_sections_)
	{
		Eigen::Vector2d half2_end_point;
		cross_section_surface_.getValue(iCrossSection.u_, 2.0, half2_end_point);
		iCrossSection.cross_section_.setHalf2EndPoint(half2_end_point);
	}

}

Eigen::VectorXd CushionSurface::SurfaceSamplePointFunction(const Eigen::VectorXd& _key_parameters)
{

	// get key cross-section from x
	std::vector<CrossSection> key_cross_sections(key_cross_sections_.size());
	for (int iKey = 0; iKey < key_cross_sections.size(); ++iKey)
	{
		key_cross_sections[iKey] = key_cross_sections_[iKey].cross_section_;
	}
	int half1_degree = key_cross_sections_[0].cross_section_.getHalf1Degree();
	int half2_degree = key_cross_sections_[0].cross_section_.getHalf2Degree();
	int dim_each_cross_section = 2 * (half1_degree - 1) + 2 * (half2_degree - 2) + 5;

	for (int iSection = 0; iSection < key_cross_sections_.size(); ++iSection)
	{
		std::vector<Eigen::Vector2d> half1_rest(half1_degree - 1);
		double half1_point1_length = 0.0;
		double tangent_angle = 0.0;
		double half2_point1_length = 0.0;
		std::vector<Eigen::Vector2d> half2_rest(half2_degree - 2);//except end point

		int base_index = iSection * dim_each_cross_section;
		half1_point1_length = _key_parameters[base_index];
		tangent_angle = _key_parameters[base_index + 1];
		half2_point1_length = _key_parameters[base_index + 2];

		int start_index_half1 = base_index + 3;
		for (int iRest1 = 0; iRest1 < half1_rest.size(); ++iRest1)
		{
			half1_rest[iRest1].x() = _key_parameters[start_index_half1 + 2 * iRest1];
			half1_rest[iRest1].y() = _key_parameters[start_index_half1 + 2 * iRest1 + 1];
		}
		int start_index_half2 = start_index_half1 + 2 * half1_rest.size();
		for (int iRest2 = 0; iRest2 < half2_rest.size(); ++iRest2)
		{
			half2_rest[iRest2].x() = _key_parameters[start_index_half2 + 2 * iRest2];
			half2_rest[iRest2].y() = _key_parameters[start_index_half2 + 2 * iRest2 + 1];
		}
		int start_index_point0 = start_index_half2 + 2 * half2_rest.size();
		Eigen::Vector2d point0(_key_parameters[start_index_point0], _key_parameters[start_index_point0 + 1]);
		key_cross_sections[iSection].setControlPoints_SameDegree(half1_rest, half1_point1_length, tangent_angle, half2_point1_length, half2_rest, point0);
	}

	changeKeyCrossSections(key_cross_sections);

	ConstructCrossSectionSurface(3);

	computeDiscreteCrossSections_NewFrames(v_values_.size() - 1);// TODO this can be speed up by avoid repeatedly constructing bezier curve at the same frame

	//UpdateDerivativeInfo();

	Eigen::VectorXd sample_points(u_values_.size() * v_values_.size() * 3);
	for (int iU = 0; iU < u_values_.size(); ++iU)
	{
		for (int iV = 0; iV < v_values_.size(); ++iV)
		{
			int point_id = iU * v_values_.size() + iV;
			sample_points(3 * point_id + 0) = discrete_cushion_[iU].discrete_cross_section_[iV](0);
			sample_points(3 * point_id + 1) = discrete_cushion_[iU].discrete_cross_section_[iV](1);
			sample_points(3 * point_id + 2) = discrete_cushion_[iU].discrete_cross_section_[iV](2);
		}
	}

	return sample_points;
}

Eigen::VectorXd CushionSurface::SurfaceNormalSamplePointFunction(const Eigen::VectorXd& _key_parameters)
{
	// get key cross-section from x
	std::vector<CrossSection> key_cross_sections(key_cross_sections_.size());
	for (int iKey = 0; iKey < key_cross_sections.size(); ++iKey)
	{
		key_cross_sections[iKey] = key_cross_sections_[iKey].cross_section_;
	}
	int half1_degree = key_cross_sections_[0].cross_section_.getHalf1Degree();
	int half2_degree = key_cross_sections_[0].cross_section_.getHalf2Degree();
	int dim_each_cross_section = 2 * (half1_degree - 1) + 2 * (half2_degree - 2) + 5;

	for (int iSection = 0; iSection < key_cross_sections_.size(); ++iSection)
	{
		std::vector<Eigen::Vector2d> half1_rest(half1_degree - 1);
		double half1_point1_length = 0.0;
		double tangent_angle = 0.0;
		double half2_point1_length = 0.0;
		std::vector<Eigen::Vector2d> half2_rest(half2_degree - 2);//except end point

		int base_index = iSection * dim_each_cross_section;
		half1_point1_length = _key_parameters[base_index];
		tangent_angle = _key_parameters[base_index + 1];
		half2_point1_length = _key_parameters[base_index + 2];

		int start_index_half1 = base_index + 3;
		for (int iRest1 = 0; iRest1 < half1_rest.size(); ++iRest1)
		{
			half1_rest[iRest1].x() = _key_parameters[start_index_half1 + 2 * iRest1];
			half1_rest[iRest1].y() = _key_parameters[start_index_half1 + 2 * iRest1 + 1];
		}
		int start_index_half2 = start_index_half1 + 2 * half1_rest.size();
		for (int iRest2 = 0; iRest2 < half2_rest.size(); ++iRest2)
		{
			half2_rest[iRest2].x() = _key_parameters[start_index_half2 + 2 * iRest2];
			half2_rest[iRest2].y() = _key_parameters[start_index_half2 + 2 * iRest2 + 1];
		}
		int start_index_point0 = start_index_half2 + 2 * half2_rest.size();
		Eigen::Vector2d point0(_key_parameters[start_index_point0], _key_parameters[start_index_point0 + 1]);
		key_cross_sections[iSection].setControlPoints_SameDegree(half1_rest, half1_point1_length, tangent_angle, half2_point1_length, half2_rest, point0);

	}

	changeKeyCrossSections(key_cross_sections);

	ConstructCrossSectionSurface(3);

	computeDiscreteCrossSections_NewFrames(v_values_.size() - 1);// TODO this can be speed up by avoid repeatedly constructing bezier curve at the same frame

	//UpdateDerivativeInfo();

	std::vector<std::vector<Eigen::Vector3d>> normal_global;
	FindSurfaceNormal(normal_global);
	Eigen::VectorXd sample_points(u_values_.size() * v_values_.size() * 3);
	for (int iU = 0; iU < u_values_.size(); ++iU)
	{
		for (int iV = 0; iV < v_values_.size(); ++iV)
		{
			int point_id = iU * v_values_.size() + iV;
			sample_points(3 * point_id + 0) = discrete_cushion_[iU].frame_value_(0) + normal_global[iU][iV](0);
			sample_points(3 * point_id + 1) = discrete_cushion_[iU].frame_value_(1) + normal_global[iU][iV](1);
			sample_points(3 * point_id + 2) = discrete_cushion_[iU].frame_value_(2) + normal_global[iU][iV](2);

		}
	}

	return sample_points;

}
