#include "CurveDistanceToFace.h"
#include "../CurveNSurface/SamplePolylineFromCurve.h"
CurveDistanceToFace::CurveDistanceToFace()
{
}

CurveDistanceToFace::~CurveDistanceToFace()
{
}

void CurveDistanceToFace::setHead(const tmd::TriangleMeshDistance* _head)
{
	head_ = _head;
}

double CurveDistanceToFace::Distance(const Rail& _rail) const
{
	std::vector<Mesh::Point> rail_points;
	SamplePolylineFromCurve getSample;
	getSample.getPolylineAsPointVector3D<Mesh::Point>(_rail.getRail(), 300, rail_points);

	return Distance(rail_points);
}

double CurveDistanceToFace::Distance(const std::vector<Mesh::Point>& _curve) const
{
	double average_distance = 0.0;
	for (int iPoint = 0; iPoint < _curve.size(); ++iPoint)
	{
		average_distance += head_->signed_distance({ _curve[iPoint][0], _curve[iPoint][1] ,_curve[iPoint][2] }).distance;
	}
	average_distance /= double(_curve.size());
	return average_distance;
}

double CurveDistanceToFace::MinDistance(const std::vector<Mesh::Point>& _curve) const
{
	double min_distance = FLT_MAX;
	for (int iPoint = 0; iPoint < _curve.size(); ++iPoint)
	{
		double distance = head_->signed_distance({ _curve[iPoint][0], _curve[iPoint][1] ,_curve[iPoint][2] }).distance;
		if (distance < min_distance)
		{
			min_distance = distance;
		}
	}
	return min_distance;
}
