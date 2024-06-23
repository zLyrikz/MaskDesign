#pragma once
#include "../MeshViewer/MeshDefinition.h"

class MeshInfo
{
public:
	MeshInfo(const Mesh* _mesh);

	void getBoundingBox(Mesh::Point& _bbMin, Mesh::Point& _bbMax, bool _add_epsilon = false) const;

	double getFaceMaxCoordinate(int _face_id, int _coordinate) const;// coordinate = 0, 1, or 2
	double getFaceMinCoordinate(int _face_id, int _coordinate) const;

	double getAverageEdgeLength() const;// a coarse measure of mesh face size(length)

	void getVerticesAndTriangles(std::vector<std::array<double, 3>>& _vertices, std::vector<std::array<int, 3>>& _triangles) const;

	void getInitialBoundaryHalfedge(std::vector<Mesh::HalfedgeHandle>& _initial_halfedges) const;
private:
	const Mesh* mesh_;
};

