#pragma once
#include "Rail.h"
#include "../MeshProcessing/TriangleMeshDistance.h"

class CurveDistanceToFace
{
public:
	CurveDistanceToFace();
	~CurveDistanceToFace();

	void setHead(const tmd::TriangleMeshDistance* _head);
	double Distance(const Rail& _rail) const;
	double Distance(const std::vector<Mesh::Point>& _curve) const;
	double MinDistance(const std::vector<Mesh::Point>& _curve) const;

private:
	const tmd::TriangleMeshDistance* head_ = nullptr;
};

