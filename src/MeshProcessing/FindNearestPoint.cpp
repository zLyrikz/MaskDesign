#include "FindNearestPoint.h"
#include <iostream>
using std::cout;
using std::endl;
FindNearestPoint::FindNearestPoint()
{
}

FindNearestPoint::~FindNearestPoint()
{
	if (points_tree_ != nullptr)
	{
		delete points_tree_;
		points_tree_ = nullptr;
	}
}

void FindNearestPoint::setPoints(const std::vector<Mesh::Point>& _points)
{
	unsigned num_point = _points.size();
	ANNpointArray dataPts = annAllocPts(num_point, 3);
	{
		int count = 0;
		for (const Mesh::Point& iPoint : _points)
		{
			dataPts[count][0] = iPoint[0]; dataPts[count][1] = iPoint[1]; dataPts[count][2] = iPoint[2];
			++count;
		}
	}

	if (points_tree_ != nullptr)
	{
		delete points_tree_;
		points_tree_ = nullptr;
	}
	points_tree_ = new ANNkd_tree(dataPts, num_point, 3);// ANN lib related question: why not working if not new a pointer?
}

int FindNearestPoint::AnnInPoints(const Mesh::Point& _point)
{
	int id_found = -1;

	if (points_tree_ != nullptr)
	{
		ANNpoint tp = annAllocPt(3); tp[0] = _point[0]; tp[1] = _point[1]; tp[2] = _point[2];
		ANNidxArray nnIdx = new ANNidx[1]; ANNdistArray dists = new ANNdist[1];

		points_tree_->annkSearch(tp, 1, nnIdx, dists);

		id_found = nnIdx[0];

		delete[] nnIdx;
		delete[] dists;
	}
	else
	{
		cout << "[WARNING FROM FindNearestPoint::AnnInPoints] points tree not set" << endl;
	}

	return id_found;
}
