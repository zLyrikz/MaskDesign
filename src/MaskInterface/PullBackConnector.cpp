#include "PullBackConnector.h"
#include "../MeshProcessing/SubMesh.h"
#include <iostream>
using std::cout;
using std::endl;
PullBackConnector::PullBackConnector()
{
}

PullBackConnector::~PullBackConnector()
{
}

double PullBackConnector::TranslateAlongZ(double _target_length, const CushionSurface& _cushion, 
	const Mesh& _connector_curve, int _curve_sample) const
{
	// convert boundary data type to vector
	std::vector<Mesh::Point> connector_curve;
	SubMesh convert_type;
	convert_type.PointCloudMeshToPoints(_connector_curve, connector_curve);

	double final_translation_length = 0.0;// sum of all translation
	double translate_length = 0.0;// translation in each iteration
	double initial_distance = 0.0;
	for (int iIteration = 0; iIteration < 10; ++iIteration)
	{
		// find distance
		std::vector<int> sample_frames;
		std::vector<Mesh::Point> intersect_points;
		_cushion.FindConnectorAndSampleFrameIntersectPoints(connector_curve, _curve_sample, sample_frames, intersect_points);
		// distance = point.z - rail.z;
		double distance = FindMinZDistanceOfRailCurveToPoints(sample_frames, intersect_points, _cushion);
		if (iIteration == 0)
		{
			initial_distance = distance;
		}

		// find translation
		translate_length = _target_length - distance;
		if (abs(translate_length) < 1e-5)
		{
			// no need more iteration if translate_length already small
			break;
		}

		// translate
		TranslatePointsAlongZ(connector_curve, translate_length);
		final_translation_length += translate_length;
	}


	// find final distance
	std::vector<int> sample_frames;
	std::vector<Mesh::Point> intersect_points;
	_cushion.FindConnectorAndSampleFrameIntersectPoints(connector_curve, _curve_sample, sample_frames, intersect_points);
	double final_distance = FindMinZDistanceOfRailCurveToPoints(sample_frames, intersect_points, _cushion);

	cout << "target min distance=" << _target_length 
		<< ", initial distance=" << initial_distance 
		<< ", translate length=" << final_translation_length 
		<< ", min distance if after translation=" << final_distance << endl;

	return final_translation_length;
}

void PullBackConnector::TranslateMeshAlongZ(Mesh& _mesh, double _translation_length) const
{
	for (auto& iVertex : _mesh.vertices())
	{
		Mesh::Point new_point = _mesh.point(iVertex);
		new_point[2] += _translation_length;
		_mesh.set_point(iVertex, new_point);
	}
}

void PullBackConnector::TranslatePointsAlongZ(std::vector<Mesh::Point>& _points, double _translate_length) const
{
	for (auto& iPoint : _points)
	{
		iPoint[2] += _translate_length;
	}
}

double PullBackConnector::FindMinZDistanceOfRailCurveToPoints(const std::vector<int>& _sample_frames, const std::vector<Mesh::Point>& _points, const CushionSurface& _cushion) const
{
	double min_distance = FLT_MAX;
	for (int iFrame = 0; iFrame < _sample_frames.size(); ++iFrame)
	{
		Eigen::Vector3d rail_value;
		_cushion.getRailValue(_sample_frames[iFrame], rail_value);

		double distance = _points[iFrame][2] - rail_value[2];
		if (distance < min_distance)
		{
			min_distance = distance;
		}
	}
	return min_distance;
}
