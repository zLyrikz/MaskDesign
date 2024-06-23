#pragma once
#include "../MeshViewer/MeshDefinition.h"
#include "../ANN/ANN.h"

// notice copy constructor is shallow copy now
class FindNearestPoint
{
public:
	FindNearestPoint();
	~FindNearestPoint();

	void setPoints(const std::vector<Mesh::Point>& _points);
	int AnnInPoints(const Mesh::Point& _point);
		
private:
	ANNkd_tree* points_tree_ = nullptr;
};

