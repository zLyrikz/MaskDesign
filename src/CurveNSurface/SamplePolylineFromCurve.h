#pragma once
#include <array>
#include "../Alglib/interpolation.h"
#include "../MeshViewer/MeshDefinition.h"
#include "../CurveNSurface/BSpline.h"

class SamplePolylineFromCurve
{
public:
	SamplePolylineFromCurve();
	~SamplePolylineFromCurve();

	// 3d curve, domain = [0,1] 
	void getPolylineAsPointCloud3D(const std::array<alglib::spline1dinterpolant, 3>& _curve, int _sample_num, Mesh& _polyline);

	// 3d curve, set domain manually
	void getPolylineAsPointCloud3D(const BSpline& _curve, int _sample_num, Mesh& _polyline, const std::array<double, 2>& _interval = std::array<double, 2>{0, 1});
	// TODO: to be tested
	// Point3D like Mesh::Point or CGAL::Point, and it has constructor like Point3D(x,y,z)
	template<class Point3D>
	void getPolylineAsPointVector3D(const BSpline& _curve, int _sample_num, std::vector<Point3D>& _polyline, const std::array<double, 2>& _interval = std::array<double, 2>{0, 1});
	
	// 2d curve, map to 3d, x=0,y=x,z=y
	void getPolylineAsPointCloud2D(const BSpline& _curve, int _sample_num, Mesh& _polyline, const std::array<double,2>& _interval= std::array<double, 2>{0,1});

	// 1d curve, set domain manually; mesh z=0, x=curve parameter, y=curve value
	void getPolylineAsPointCloud1D(const BSpline& _curve, int _sample_num, Mesh& _polyline, const std::array<double, 2>& _interval = std::array<double, 2>{0, 1});
	void getDerivatePolylineAsPointCloud1D(const BSpline& _curve, int _sample_num, Mesh& _polyline, const std::array<double, 2>& _interval = std::array<double, 2>{0, 1});
};

template<class Point3D>
inline void SamplePolylineFromCurve::getPolylineAsPointVector3D(const BSpline& _curve, int _sample_num, std::vector<Point3D>& _polyline, const std::array<double, 2>& _interval)
{
	//// get the mesh of the curve
	_polyline.reserve(_sample_num);
	for (int iSample = 0; iSample < _sample_num; ++iSample)
	{
		float para = (float)iSample / (_sample_num - 1);
		para *= (_interval[1] - _interval[0]);
		para += _interval[0];
		Eigen::VectorXd value;
		_curve.getValue(para, value);
		_polyline.push_back(std::move(Point3D(value(0), value(1), value(2))));
	}
}
