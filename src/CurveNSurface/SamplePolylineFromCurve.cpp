#include "SamplePolylineFromCurve.h"
#include "../MeshProcessing/SubMesh.h"

SamplePolylineFromCurve::SamplePolylineFromCurve()
{
}

SamplePolylineFromCurve::~SamplePolylineFromCurve()
{
}

void SamplePolylineFromCurve::getPolylineAsPointCloud3D(const std::array<alglib::spline1dinterpolant, 3>& _curve, int _sample_num, Mesh& _polyline)
{
	//// get the mesh of the curve
	std::vector<std::array<double, 3>> vertices(_sample_num);
	std::vector<std::array<int, 3>> triangles;
	for (int iSample = 0; iSample < _sample_num; ++iSample)
	{
		float para = (float)iSample / (_sample_num - 1);
		vertices[iSample][0] = alglib::spline1dcalc(_curve[0], para);
		vertices[iSample][1] = alglib::spline1dcalc(_curve[1], para);
		vertices[iSample][2] = alglib::spline1dcalc(_curve[2], para);
	}
	SubMesh construct_mesh;

	construct_mesh.MeshFromVerticesAndTriangles(_polyline, vertices, triangles);
}

void SamplePolylineFromCurve::getPolylineAsPointCloud3D(const BSpline& _curve, int _sample_num, Mesh& _polyline, const std::array<double, 2>& _interval)
{
	//// get the mesh of the curve
	std::vector<std::array<double, 3>> vertices(_sample_num);
	std::vector<std::array<int, 3>> triangles;
	for (int iSample = 0; iSample < _sample_num; ++iSample)
	{
		float para = (float)iSample / (_sample_num - 1);
		para *= (_interval[1] - _interval[0]);
		para += _interval[0];
		Eigen::VectorXd value;
		_curve.getValue(para, value);
		vertices[iSample][0] = value(0);
		vertices[iSample][1] = value(1);
		vertices[iSample][2] = value(2);
	}
	SubMesh construct_mesh;

	construct_mesh.MeshFromVerticesAndTriangles(_polyline, vertices, triangles);
}



void SamplePolylineFromCurve::getPolylineAsPointCloud2D(const BSpline& _curve, int _sample_num, Mesh& _polyline, const std::array<double, 2>& _interval)
{
	//// get the mesh of the curve
	std::vector<std::array<double, 3>> vertices(_sample_num);
	std::vector<std::array<int, 3>> triangles;
	for (int iSample = 0; iSample < _sample_num; ++iSample)
	{
		float para = (float)iSample / (_sample_num - 1);
		para *= (_interval[1] - _interval[0]);
		para += _interval[0];
		Eigen::VectorXd value;
		_curve.getValue(para, value);
		vertices[iSample][0] = 0.0;
		vertices[iSample][1] = value(0);
		vertices[iSample][2] = value(1);
	}
	SubMesh construct_mesh;

	construct_mesh.MeshFromVerticesAndTriangles(_polyline, vertices, triangles);
}

void SamplePolylineFromCurve::getPolylineAsPointCloud1D(const BSpline& _curve, int _sample_num, Mesh& _polyline, const std::array<double, 2>& _interval)
{
	//// get the mesh of the curve
	std::vector<std::array<double, 3>> vertices(_sample_num);
	std::vector<std::array<int, 3>> triangles;
	for (int iSample = 0; iSample < _sample_num; ++iSample)
	{
		float para = (float)iSample / (_sample_num - 1);
		para *= (_interval[1] - _interval[0]);
		para += _interval[0];
		Eigen::VectorXd value;
		_curve.getValue(para, value);
		vertices[iSample][0] = para;
		vertices[iSample][1] = value(0);
		vertices[iSample][2] = 0.0;
	}
	SubMesh construct_mesh;

	construct_mesh.MeshFromVerticesAndTriangles(_polyline, vertices, triangles);
}

void SamplePolylineFromCurve::getDerivatePolylineAsPointCloud1D(const BSpline& _curve, int _sample_num, Mesh& _polyline, const std::array<double, 2>& _interval)
{
	//// get the mesh of the curve
	std::vector<std::array<double, 3>> vertices(_sample_num);
	std::vector<std::array<int, 3>> triangles;
	for (int iSample = 0; iSample < _sample_num; ++iSample)
	{
		float para = (float)iSample / (_sample_num - 1);
		para *= (_interval[1] - _interval[0]);
		para += _interval[0];
		Eigen::VectorXd value;
		_curve.getDerivative(para, value);
		vertices[iSample][0] = para;
		vertices[iSample][1] = value(0);
		vertices[iSample][2] = 0.0;
	}
	SubMesh construct_mesh;

	construct_mesh.MeshFromVerticesAndTriangles(_polyline, vertices, triangles);
}

