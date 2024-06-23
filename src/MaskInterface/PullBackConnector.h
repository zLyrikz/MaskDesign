#pragma once
#include "CushionSurface.h"
class PullBackConnector
{
public:
	PullBackConnector();
	~PullBackConnector();

	// find the translation length of the connector along a given direction
	// target length is the minimal connector to rail curve length that we want to reach
	double TranslateAlongZ(double _target_length, const CushionSurface& _cushion, const Mesh& _connector_curve, int _curve_sample) const;
	void TranslateMeshAlongZ(Mesh& _mesh, double _translation_length) const;
	void TranslatePointsAlongZ(std::vector<Mesh::Point>& _points, double _translate_length) const;

private:
	double FindMinZDistanceOfRailCurveToPoints(const std::vector<int>& _sample_frames, const std::vector<Mesh::Point>& _points, const CushionSurface& _cushion) const;

};

