#pragma once

#include "CrossSection.h"
#include "Rail.h"
#include "CrossSectionSurface.h"
#include "DiscreteCushionPiece.h"
#include "../MeshProcessing/AabbTree.h"
#include "../MeshViewer/MeshDefinition.h"
#include <Eigen/Core>

/*
cushion surface is modeled as a generalized sweep surface: 
s(u,v) = r(u) + M(u) * c(u,v); v in [0,2]
r(u) is the rail; M(u) is the rotation matrix;
at each u, c(u,v) is the cross section; along the u, interpolated from some sample key cross sections

process of fitting a generic cushion :
1. set rail (i.e. r(u)):
	1) set a polyline extracted from the cushion (setCushionPolyline)
	2) fit a rail from the polyline (setRailByFittingCushionPolyline)
2. find a frame (i.e. M(u)): computeSampleFrames(); this function will find a proper relaxed frame
3. construct cross section surface:
	1) set key cross section parameters on the rail (setKeyCrossSections_NewFrame(); using the relaxed frame here)
	2) fit a cross section surface from the key cross sections (ConstructCrossSectionSurface())
4. get a surface wireframe for creating a mesh object
	0) resample the frame if you want a frame different from the sampled frames(in 2.) (ResampleFrames())
	1) computeDiscreteCrossSections_NewFrames
	2) constructSurfaceMesh
*/
/*
* orientation convention:
* when looking at cushion in the wearing position, (along z axis)
* rail curve is counterclock wise; normal is pointing outward( UnitZ x Tangent )
* cross section orients from inside to outside
*/


class CushionSurface
{
public:
	CushionSurface();
	~CushionSurface();

	void setHeadTree(const AabbTree* _head_tree);
	void setCushionMesh(const Mesh* _cushion);
	void freeCushionMesh();

	void setCushionPolyline(const std::vector<Eigen::VectorXd>& _cushion_polyline);
	void setRailByFittingCushionPolyline(int _segments, int _degree, int _smoothness);
	void setRailByFittingCushionPolyline(const std::vector<int>& _parts_idx, int _degree, int _smoothness);// manually define a non uniform segment 
	void changeRailByFittingProjectedCushionPolylineOnFace(int _segments, int _degree, int _smoothness);//project direction is unit z
	void changeRailByFittingProjectedCushionPolylineOnFace(const std::vector<int>& _parts_idx, int _degree, int _smoothness);//project direction is unit z
	void setRailByFitting(const std::vector<Eigen::VectorXd>& _polyline, int _segments, int _degree, int _smoothness);
	void setRailByFitting(const std::vector<Eigen::VectorXd>& _polyline, const std::vector<int>& _parts_idx, int _degree, int _smoothness);
	void setRail(const BSpline& _rail);
	void setCrossSectionSurface(
		const BSpline& _half1_rest, const BSpline& _half1_point1_length,
		const BSpline& _tangent,
		const BSpline& _half2_point1_length, const BSpline& _half2_rest);
	void setLineSegmentAsCrossSectionSurface(double _line_radius);// this is for visualization of the relaxed frame

	// the standard pipeline to fit a cushion(can be original or other cushion for recover the parametric representation)
	void FitCushion(const Mesh* _cushion_surface, 
		const std::vector<Eigen::VectorXd>& _rail_polygon, const std::vector<int>& _polygon_segment_parts_idx, const std::vector<double>& _key_cross_section_sample,
		int _sample_frames, int _sample_v, const Mesh& _connector_boundary,
		int _half1_degree, int _half2_degree);

	// number of sample = number of u sample in the final result mesh
	void computeSampleFrames(float _radius, int _num_sample_u = 1000);
	void computeStandardSampleFrames(int _num_sample_u = 1000);// not doing frame relaxation
	void computeInterpolatedFrameAtU(double _u, int _frame0_idx, int _frame1_idx, DiscreteCushionPiece& _cushion_piece) const;// input u in [u0,u1]; here u0 is frame0 parameter and u1 is frame1 parameter
	void ResampleFrames(int _num_sample_u = 1300);// after sample frames and done the relaxation, resample it by interpolation the tangent rotation

	// interpolate the boundary, so that surface keeps the boundary
	// i.e. find the boundary 2d spline formed by control points in CrossSectionSurface
	// rail and frame must been set already
	// also translate connector back if it's too close to the human face, according to target translation length:_pull_back_connector
	void setFixedBoundary(const std::string& _original_boundary, int _num_interpolate_point = 100, double _pull_back_connector = 0.0);
	void setFixedHalf1Boundary(const std::string& _original_boundary, int _num_interpolate_point = 50);// TODO: copy from setFixedBoundary, need to combine the code
	void setFixedBoundary(const Mesh& _original_boundary, int _num_interpolate_point = 101, double _pull_back_connector = 0.0);
	// utility function for setFixedBoundary
	// input sample number, output sample frame index and intersect points
	// NOTE for convenience, the real sample number maybe slightly different from the input number
	void FindConnectorAndSampleFrameIntersectPoints(const std::vector<Mesh::Point>& _connector_polygon,
		int _reference_num_sample, std::vector<int>& _sample_frames, std::vector<Mesh::Point>& _intersect_point) const;
	void getHalf2Boundary(std::vector<Mesh::Point>& _boundary) const;// sample points of the boundary curve

	// this is for fitting the original cushion
	// 1. fit the rail first 
	// 2. at given u, intersect the rail normal plane with the cushion, get a 2D polyline
	// 3. fit the polyline, get the cross section (two joined 2D Bezier curve) 
	void setKeyCrossSections(const std::vector<double>& _u, int _degree = 3);// initial implementation(using standard frames)
	void setKeyCrossSections_NewFrame(const std::vector<int>& _i_frame, int _degree = 3);// after sample a surface wire frame, set key cross sections by frame index
	void setKeyCrossSections_NewFrame(const std::vector<double>& _u, int _degree = 3);// after sample a surface wire frame, set key cross sections by finding closest frame index from the u parameters
	void setKeyCrossSections_NewFrame_ChangableDegree(const std::vector<double>& _u, int _degree_half1, int _degree_half2);
	void setKeyCrossSections_NewFrame_ChangableDegree_FindHalf2BeginIndex(const std::vector<double>& _u, const std::vector<Mesh::Point>& _half2_endpoints, int _degree_half1, int _degree_half2);
	// this is for finding key cross sections for a general rail 
	void setInitialCustomKeyCrossSections(const std::vector<int>& _i_frame, const Mesh& _original_cushion, const std::vector<Mesh::Point>& _original_rail_polyline, int _degree = 3);
	void setInitialCustomKeyCrossSections(const std::vector<double>& _u, int _degree = 3);
	void setInitialCustomKeyCrossSections_ChangableDegree(const std::vector<double>& _u, int _degree_half1, int _degree_half2);
	void changeKeyCrossSections(const std::vector<CrossSection>& _key_cross_sections);// not change the u, only change the cross section
	// Major difference is this does not translate the connect point of two halfs to origin  
	void setKeyCrossSectionFitAnotherCushion(const Mesh& _cushion, const std::vector<Mesh::Point>& _rail_polyline, int _degree = 2);// TODO: mostly copy from setInitialCustomKeyCrossSections function
	void setKeyCrossSection(const std::vector<CrossSection>& _key_cross_sections);

	//void StretchKeyCrossSection_BoundaryFixed();
	

	// interpolate the Bezier control points of key cross sections in 2D plane
	// degree is the interpolation curve degree
	// not compute the boundary control points, it should be fixed as the original boundary
	void ConstructCrossSectionSurface(int _degree);

	// standard frames (take frame tangent normal as rail tangent, or say cross section in the curve normal plane)
	void constructStandardSurfaceWireFrame(int _num_sample_u = 1000, int _num_sample_v = 50);
	// use circle with radius c as bound, test self intersection of the standard frame, then do a relaxation to get a no self intersected frame
	//void constructRelaxedSurfaceWireFrame(float _radius, const std::vector<double>& _u, int _cross_section_degree = 3, int _cross_section_surface_degree = 3, int _num_sample_u = 1000, int _num_sample_v = 50);

	void ConstructSelfIntersectPartMesh(std::vector<Mesh>& _part_meshes) const;
	void FrameRelaxation(int _interval_side_expand = 2, bool _extra_iteration = false);

	void UpdateFrameValue();
	void computeDiscreteCrossSections_NewFrames(int _num_sample_v = 50);// after computing the frames, compute sampled the cross section values in each frame
	void FindCrossSectionForAllFrames(std::vector<CrossSection>& _cross_sections) const;

	// constructSurfaceWireFrame beforehand
	void constructVUPoints(std::vector<std::vector<Mesh::Point>>& _vu_points) const;// this representation is good for silicone printing
	void constructUVPoints(std::vector<std::vector<Mesh::Point>>& _uv_points) const;
	void constructSurfaceMesh(/*int _num_sample_u = 1000, int _num_sample_v = 50*/);//construct mesh from the sampled surface wireframe 
	void constructSmoothedSurfaceMesh(int _num_sample_u = 1000, int _num_sample_v = 50, int _smoothness = 6);// (older ideas)using face normal info for defining rail normal  
	const Mesh& getSurfaceMesh() const;
	void setCurrentMeshAsOriginalCushionMesh();
	void setOriginalRail(const Rail& _rail);// original rail is used in setInitialCustomKeyCrossSections functions
	void setOriginalCushionMesh(const Mesh& _cushion_surface);// original cushion surface is used in setInitialCustomKeyCrossSections functions

	const std::vector<DiscreteCushionPiece>* getDiscreteCushion() const;

	void getValue_FaceNormal(double _u, double _v, Mesh::Point& _value) const;// TODO: delete/update this previous version
	void getValue(double _u, double _v, const Eigen::Vector3d& _rail_normal, Mesh::Point& _value) const;// TODO: delete/update this previous version
	void getValue(double _u, double _v,
		const Eigen::Vector3d& _rail_value, const Eigen::Vector3d& _frame_tangent, const Eigen::Vector3d& _frame_normal,
		Eigen::Vector3d& _value) const;// get value with the local frame given
	void getValue(double _u, double _v,
		const Eigen::Vector3d& _rail_value, const Eigen::Vector3d& _frame_tangent, const Eigen::Vector3d& _frame_normal, const Eigen::Vector3d& _frame_binormal,
		Eigen::Vector3d& _value) const;// get value with the local frame given
	void getRailValue(int _frame, Eigen::Vector3d& _value) const;

	// derivative

	void UpdateDerivativeInfo();
	// The descrete points r(u) + M2(u)c1(u,v) + M3(u)c2(u,v)
	// and cross-section normal: r(u) - M2(u)c2'(u,v) + M3(u)c1'(u,v)
	// derivative w.r.t. key control point parameters
	// pre-request: call computeDiscreteCrossSections_NewFrames() and UpdateDerivativeInfo()
	// note _half2_rest not include the end point
	void getDerivative_DiscreteCushionAndCrossSectionNormal_KeyControlPoints(
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
		std::vector<std::vector<std::array<Eigen::Matrix3Xd, 2>>>& _normal_point0) const;

	// this function is basically a copy of CrossSectionOptimizeInfo::ChangeCushionFromX 
	Eigen::VectorXd SurfaceSamplePointFunction(const Eigen::VectorXd& _key_parameters);
	Eigen::VectorXd SurfaceNormalSamplePointFunction(const Eigen::VectorXd& _key_parameters);


	// for visualization
	void getKeyCrossSectionsPolyline(std::vector<std::vector<Mesh::Point>>& _cross_section_polyline) const;
	void getKeyCrossSectionsRailPoint(std::vector<Mesh::Point>& _rail_point) const;// this function is currently unavailable
	void getKeyCrossSections(std::vector<CrossSection>& _key_cross_sections) const;
	void getKeyCrossSectionsU(std::vector<double>& _key_cross_sections_u) const;
	void getKeyCrossSectionsMesh_3D(std::vector<Mesh>& _curves) const;
	void getKeyCrossSection3DBeizerPoints(Mesh& _points) const;
	int getKeyCrossSectionsNumber() const;
	const CrossSectionSurface& getCrossSectionSurface() const;


	// normal is the cross product of rail tangent and the normal of the face
	// other options -- 
	// option1: average tangent of the curve which is the intersection of a square on the normal plane of rail and the face 
	// option2: cross product of rail tangent and axis z
	void getRailNormal_FromFace(const Mesh::Point& _rail_tangent, const Mesh::Point& _rail_value, Eigen::Vector3d& _normal) const;
	void getRailFrame_FromFace(double _u, Eigen::Vector3d& _rail_value, Eigen::Vector3d& _rail_tangent, Eigen::Vector3d& _frame_normal) const;
	void getRailNormal_AxisZ(const Eigen::Vector3d& _rail_tangent, Eigen::Vector3d& _normal) const;
	// new version: project tangent to XOY plane!!
	void getRailFrame_AxisZ(double _u, Eigen::Vector3d& _rail_value, Eigen::Vector3d& _rail_tangent, Eigen::Vector3d& _frame_normal) const;
	void getFrameGivenTangent_AxisZ(double _u, const Eigen::Vector3d& _frame_tangent, Eigen::Vector3d& _rail_value, Eigen::Vector3d& _frame_normal, Eigen::Vector3d& _frame_binormal) const;
	void getRailFrameWithSmoothedFaceNormals(int _num_sample_u,
		std::vector < Eigen::Vector3d>& _rail_value, 
		std::vector < Eigen::Vector3d>& _rail_tangent,
		std::vector < Eigen::Vector3d>& _frame_normal,
		std::vector < Eigen::Vector3d>& _frame_binormal,
		int _neighbor = 4) const;
	const Rail& getRail() const;
	const Rail& getOriginalRail() const;

	double getConnectorPullBackLength() const;
	int getHalf1Degree() const;
	int getHalf2Degree() const;

	// output normal data arrangement correspond to vertices in the discrete_cushion_
	void FindSurfaceNormal(std::vector<std::vector<Eigen::Vector3d>>& _node_normals) const;
	// _filename should not include suffix
	// not const because have to find key cross section frame id
	void WriteSampleCurves(std::string _filepath, std::string _filename);
	void SaveResults(std::string _filepath, std::string _filename) const;// TODO use XmlParser to read and save cushion raw data
	void ReadResults(std::string _xml_file, int _generalized_cushion = 0);
	void InitializeAsGenericCushion(std::string _data_file);

public:
	// some useful functions (like for fabrications)

	//project direction is unit z
	void ProjectPolylineToFace(const std::vector<Eigen::VectorXd>& _polyline, std::vector<Eigen::VectorXd>& _projected_polyline) const;

	// special design for 3d printing	
	// flatten the inner half of the key cross sections
	void FlattenInnerHalf();

	double ComputeInnerBoundaryLength() const;// unit is mm

private:
	void computeStandardFrames(int _num_sample_u);// (take frame tangent normal as rail tangent, or say cross section in the curve normal plane)
	void computeNewFrames_AxisZ();// (use newly defined frame tangent to compute new frame normal as cross product of tangent and unit z vector)
	void computeDiscreteCrossSections_StandardFrames(int _num_sample_v);// after computing the frames, compute sampled the cross section values in each frame
	////void recomputeDiscreteCrossSections();// rotate the cross sections according to the rotated frame
	void CheckIntersection(const std::vector<DiscreteCushionPiece>& _cross_sections, std::vector<std::pair<int, int>>& _intersect_parts) const;

	// find closest frame index for each u value
	// require the input u in increasing order
	void ClosestFrameIdxToU(const std::vector<double>& _u, std::vector<int>& _i_frame) const;

	// return t, such that, t*u1 + (1-t)*u0 = u 
	// i.e. t=0 <=> u=u0; t=1 <=> u=u1
	double getLocalLinearParameter(double _u, double _u0, double _u1) const;

	// for fitting key cross sections
	void FindCrossSectionPolyline(const Mesh& _cushion, const Mesh::Point& _frame_tangent, const Mesh::Point& _frame_value, const Mesh::Point& _indicator, 
		std::vector<Mesh::Point>& _result_polyline) const;
	// return the center index that is newly added in the polyline
	int AddCenterToCrossSectionPolyline(std::vector<Mesh::Point>& _cross_section, const Mesh::Point& _center) const;
	void getNewFrames(int _dicrete_cross_section_idx, Eigen::Matrix4d& _new_frame) const;
	// use a trick to get the polyline in consistent outward pointing order
	void AdjustCrossSectionPolylineToPointingOutward(std::vector<Mesh::Point>& _cross_section) const;



	// setKeyCrossSection helper functions
	void getOriginalCrossSection(int _frame, std::vector<Mesh::Point>& _cross_section_polyline) const;
	void getOriginalCrossSection_OffAxis(int _frame, const Mesh& _original_cushion, const std::vector<Mesh::Point>& _original_rail_polyline,
		std::vector<Mesh::Point>& _cross_section_polyline, int& _center_index) const;

	void FindCrossSectionAtOneFrame(CrossSection& _cross_section, int _frame) const;

	void FindKeyCrossSectionFrameId();

	// set connector_pull_back_length_ well before this function call
	void setFixedBoundary(const std::string& _original_boundary, int _num_interpolate_point, bool _do_connector_pull_back);
	void setFixedBoundary(const Mesh& _boundary, int _num_interpolate_point, bool _do_connector_pull_back);


public:
	struct KeyCrossSection
	{
		double u_ = 0.0;
		CrossSection cross_section_;
		
		std::vector<Mesh::Point> cross_section_polyline_;// sampled from the original cushion mesh, used to fit a cross section 

		int frame_id_ = 0;// this shuold be updated before use, via FindKeyCrossSectionFrameId();
		// this is for visualization only
		//Mesh::Point rail_point_;// the rail point on the original cushion correspond to the cross section;
	};

private:
	const AabbTree* head_tree_;// the aabb tree structure of the head model; for providing local frame normal information
	const Mesh* cushion_;// the original cushion for reconstruction(fit key cross sections)
	std::vector<Eigen::VectorXd> cushion_polyline_;//the extracted polyline on the original cushion(as an initial rail)

	Rail original_rail_;//fitting rail of the original cushion 
	Mesh original_fitted_surface_;//fitting of the original cushion 

	Rail rail_;
	std::vector<KeyCrossSection> key_cross_sections_;
	CrossSectionSurface cross_section_surface_;

	std::vector<DiscreteCushionPiece> discrete_cushion_;// sampled surface wireframe (information of M(u) is here)
	// u and v values in the sampled discrete cushion
	// will be updated in computeDiscreteCrossSections_NewFrames
	std::vector<double> u_values_;
	std::vector<double> v_values_;
	Mesh surface_;

	// additional info
	double connector_pull_back_length_;
	int cross_section_control_curve_degree_;
};

