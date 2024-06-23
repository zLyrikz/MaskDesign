#pragma once
#include <Eigen/Core>
#include "../MeshViewer/MeshDefinition.h"

class BSpline;
/* modeled as a two 2D Bezier curve G1 continuously joined at origin
// XOY local frame is the YOZ plane in the world coordinate
*
* for both half1 and half2:
*	first conrtol point is origin by default;
*	the second conrtol point defined by the length to the origin, (the two points of the two halfs share same direction, for G1 continuous)
*
*/
class CrossSection
{
public:
	CrossSection();
	CrossSection(const CrossSection& _other);// have to define copy constructor because there is pointer as members
	CrossSection& operator=(const CrossSection& _other);
	~CrossSection();

	// fit a cross section polyline represented in 3d space (it's actually in the 2d normal plane of the frame)
	// 1. get local 2d coordinate of the polyline (polyline in the plane spanned by the second and third basis of the frame)
	// 2. fit it with the 2d Bezier curve
	// 3. do a translation, s.t. the connection point of the 2 halfs to origin
	// input the polyline and the local frame of the polyline
	// this fits the polyline with both halfs are Bezier curve
	// notice we reverse the polyline if necessary to a correct order 
	void FitOrigianl3DPolyline(const std::vector<Mesh::Point>& _polyline, const Eigen::Matrix4d& _local_frame, int _degree = 3);
	void FitOrigianl3DPolyline(const std::vector<Mesh::Point>& _polyline, const int _half2_begin_idx, const Eigen::Matrix4d& _local_frame, int _degree = 3);
	void FitOrigianl3DPolyline_ChangeableDegree(const std::vector<Mesh::Point>& _polyline, const Eigen::Matrix4d& _local_frame, int _degree_half1, int _degree_half2);
	void FitOrigianl3DPolyline_ChangeableDegree(const std::vector<Mesh::Point>& _polyline, const int _half2_begin_idx, const Eigen::Matrix4d& _local_frame, int _degree_half1, int _degree_half2);
	// notice a orthogonal transformation plus translation of the data points won't affect the least square fitting result of the bspline!
	// notice the polyline need to already have a correct order 
	void Fit3DPolyline_FixBoundary(const std::vector<Mesh::Point>& _polyline, const int _half2_begin_idx, const Eigen::Matrix4d& _local_frame, int _degree = 3);
	void Fit3DPolyline_FixBoundary_ChangeableDegree(const std::vector<Mesh::Point>& _polyline, const int _half2_begin_idx, const Eigen::Matrix4d& _local_frame, int _degree_half1, int _degree_half2);

	void setHalf1EndPoint(const Eigen::Vector2d& _end_point);
	void setHalf2EndPoint(const Eigen::Vector2d& _end_point);
	void setTangentAngle(double _angle);
	void setHalf2Point1Length(double _length);
	void setControlPoints(const std::vector<Eigen::Vector2d>& _half1_rest, double _half1_point1_length,
		double _tangent_angle,
		double _half2_point1_length, const std::vector<Eigen::Vector2d>& _half2_rest);
	void setPoint0(const Eigen::Vector2d& _point0);
	void setControlPoints_SameDegree(const std::vector<Eigen::Vector2d>& _half1_rest, double _half1_point1_length,
		double _tangent_angle,
		double _half2_point1_length, const std::vector<Eigen::Vector2d>& _half2_rest_except_end, const Eigen::Vector2d& _connect_point);
	// without setting the boundary point (end point of half2)
	// note, have to keep the degree same with the current degree
	void setControlPoints_SameDegree(const std::vector<Eigen::Vector2d>& _half1_rest, double _half1_point1_length,
		double _tangent_angle,
		double _half2_point1_length, const std::vector<Eigen::Vector2d>& _half2_rest_except_end);
	// this if for degree3, without setting the boundary point(_half2_point3)
	void setControlPoints_SameDegree(const std::vector<Eigen::Vector2d>& _half1_rest, double _half1_point1_length,
		double _tangent_angle,
		double _half2_point1_length, const Eigen::Vector2d& _half2_point2);

	// make the curve mesh x=0, y=x, z=y
	void getCurveMeshIn3D(Mesh& _curve, const std::vector<Eigen::VectorXd>& _polyline_2d) const;
	void getHalf1EndPoint(Eigen::Vector2d& _end_point) const;
	void getHalf2EndPoint(Eigen::Vector2d& _end_point) const;
	// these get functions do not give point0, find point0 from getPoint0
	void getControlPoints(
		std::vector<Eigen::Vector2d>& _half1_rest, double& _half1_point1_length,
		Eigen::Vector2d& _tangent, 
		double& _half2_point1_length, std::vector<Eigen::Vector2d>& _half2_rest) const;
	void getControlPoints(
		std::vector<Eigen::Vector2d>& _half1_rest, double& _half1_point1_length,
		double& _tangent_angle,
		double& _half2_point1_length, std::vector<Eigen::Vector2d>& _half2_rest) const;
	void getControlPoints(
		const std::vector<Eigen::Vector2d>*& _half1_rest, const double*& _half1_point1_length,
		const double*& _tangent_angle,
		const double*& _half2_point1_length, const std::vector<Eigen::Vector2d>*& _half2_rest) const;
	void getControlPoints(
		const std::vector<Eigen::Vector2d>*& _half1_rest, const double*& _half1_point1_length,
		const double*& _tangent_angle,
		const double*& _half2_point1_length, const std::vector<Eigen::Vector2d>*& _half2_rest, const Eigen::Vector2d*& _point0) const;
	// this if for degree3, without getting the boundary point(_half2_point3)
	void getControlPoints(
		const std::vector<Eigen::Vector2d>*& _half1_rest, const double*& _half1_point1_length,
		const double*& _tangent_angle,
		const double*& _half2_point1_length, const Eigen::Vector2d*& _half2_point2) const;

	double getTangentAngle() const;
	double getConnectPointX() const;
	double getHalf1Point2X() const;


	int getHalf1Degree() const;
	int getHalf2Degree() const;


	// save and read result
	void saveToFile(const std::string& _file) const;
	void loadFromFile(const std::string& _file);
public:
	// update the information before calling these functions
	void UpdateMoreInfo();

	// for evaluation
	// 
	// all the rest should be called after UpdateMoreInfo() called once
	bool IsConvex_Half1Point1() const;
	bool IsConvex_Half1Point2() const;
	bool IsConvex_Half2Point1() const;
	bool IsConvex_Half2Point2() const;
	double AngleAtHalf1Point1() const;// range in [0, pi]
	double AngleAtHalf1Point2() const;
	double AngleAtHalf2Point1() const;
	double AngleAtHalf2Point2() const;
	double AngleBetweenHalf1Vector01And23() const;	
	double GetHalf1Vector01Length() const;		
	double GetHalf1Vector12Length() const;		// TODO pre compute norm, and just return the value instead of computing the norm
	double GetHalf1Vector23Length() const;
	double GetHalf2Vector01Length() const;
	double GetHalf2Vector12Length() const;
	double GetHalf2Vector23Length() const;
	double GetHalf1EndNorm() const;
	const Eigen::Vector2d* GetHalf1Point1() const;
	const Eigen::Vector2d* GetHalf1Point2() const;
	const Eigen::Vector2d* GetHalf2Point1() const;
	const Eigen::Vector2d* GetHalf2Point2() const;

	double CrossProductHalf1Vector10AndVector12() const;
	double CrossProductHalf2Vector12AndVector10() const;
	const Eigen::Vector2d* GetHalf1Vector01() const;
	const Eigen::Vector2d* GetHalf1Vector12() const;
	const Eigen::Vector2d* GetHalf2Vector01() const;
	const Eigen::Vector2d* GetHalf2Vector12() const;
	const Eigen::Vector2d* GetPoint0() const;
	const Eigen::Vector2d* GetTangent() const;

public:
	// for visualization
	void getBezierCurve(Mesh& _bezier_curve, Mesh& _bezier_points) const;

public:
	// special design for 3d printing
	
	// flatten the inner half of the cross section
	void FlattenHalf1();
	// project the end point on the line of two mid points
	void FlattenEndPoint_Half1();
	// project rest control points to the line parallel to axis x and go through first mid point
	void FlattenHalf1_ParallelToAxisX();

private:
	void compute2DLocalCoordinateOf3DPolyline(const std::vector<Mesh::Point>& _3d_polyline, const Eigen::Matrix4d& _local_frame,
		std::vector<Eigen::VectorXd>& _2d_polyline) const;

	// suppose x values are monotonic, we separate two halfs by x negative and positive
	// half1 x <=0; half2 x > 0
	// also make the polyline contains half1 in the front part, and half2 in the back part
	// i.e. reverse the polyline if half2 is in the front part
	// at last, return index of the half2 first point in the polyline
	int SeparatePolyineTo2Halfs(std::vector<Eigen::VectorXd>& _2d_polyline) const;

	void ComputeBezierPointsFromBSpline(const BSpline& _b_spline);

	void Fit2DPolyline(const std::vector<Eigen::VectorXd>& _polyline_2d, const int _half2_begin_idx, int _degree);
	void Fit2DPolyline_ChangebaleDegree(const std::vector<Eigen::VectorXd>& _polyline_2d, const int _half2_begin_idx,
		const Eigen::Vector2d& _connect_point, int _degree_half1, int _degree_half2);

	void TranslateConnectPointToOrigin();
	void TranslateConnectPointToOrigin_FixHalf2EndPoint();

	// used for testing convex
	// note order of two vector matters 
	bool IsAngleSmallerThan180(const Eigen::Vector2d& _vector1, const Eigen::Vector2d& _vector2) const;
	double CrossProduct2D(const Eigen::Vector2d& _vector1, const Eigen::Vector2d& _vector2) const;
	double AngleAtTwoVector(const Eigen::Vector2d& _vector1, const Eigen::Vector2d& _vector2) const;// range in [0,pi]
public:
	/*******************************************************************************************************/
	// for visulization
	//// 

	// the original polyline that's fitted (before doing any translation or modification on bezier point)
	// and the fitting spline
	Mesh polyline_2d_;
	Mesh spline_;

	// the final bezier points and curve
	Mesh bezier_;
	Mesh bezier_points_;
	/*******************************************************************************************************/

private:
	// control points 
	// the origin is regarded as the first control point for both halfs
	// note for the current BMC-F2 testing example, half1 is the inside half, half2 is the outside
	std::vector<Eigen::Vector2d> half1_rest_;// size = degree - 1
	double half1_point1_length_ = 0.0;// > 0
	Eigen::Vector2d tangent_;//// note, this may not be updated accordingly!  by default, this is the direction of half2 point1; normalized
	double tangent_angle_ = 0.0;// represent the angle between tangent and unit x; angle range in [-pi, pi]
	double half2_point1_length_ = 0.0;// > 0
	std::vector<Eigen::Vector2d> half2_rest_;// size = degree - 1
	Eigen::Vector2d connect_point_;// the bezier control point that connects two halfs; it may not be zero if not translated to origin yet

	int half1_degree_ = 0;
	int half2_degree_ = 0;

	// the following info are computed from function call: UpdateMoreInfo();
	// record some data for conviently computing angle at control points
	// dynamicly allocate memory if these are needed
	// note they are not updated if control points are changed
	// coding note! using newed pointer must provide copy constructor and assignment operator
	Eigen::Vector2d* half1_point1_;
	Eigen::Vector2d* half1_vector01_;
	Eigen::Vector2d* half2_point1_;
	Eigen::Vector2d* half2_vector01_;
	Eigen::Vector2d* half1_vector12_;// the vector of (point2 - point1)
	Eigen::Vector2d* half2_vector12_;// the vector of (point2 - point1)
	Eigen::Vector2d* half1_vector23_;// the vector of (point3 - point2) for degree 3 curve
	Eigen::Vector2d* half2_vector23_;// the vector of (point3 - point2) for degree 3 curve
	double half1_end_norm_ = 0.0;
};

