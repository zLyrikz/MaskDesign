#pragma once
#include "../MeshViewer/MeshDefinition.h"
#include <Eigen/Core>

/* a class for construct meshes (various types)
*/

class SubMesh
{
public:
	SubMesh();
	~SubMesh();

	void setData(const Mesh* _whole_mesh, std::vector<int>& _face_area);
	void ConstructSubMesh(); // construct a sub area mesh if given the whole mesh and area vertices
	void getSubMesh(Mesh& _sub_mesh) const;

	// get triangle mesh in libigl Eigen mesh representation
	void Mesh2EigenMatrix(const Mesh& _mesh, Eigen::MatrixXd& _vertices, Eigen::MatrixXi& _faces) const;
	void MeshFromEigenMatrix(Mesh& _mesh, const Eigen::MatrixXd& _vertices, const Eigen::MatrixXi& _faces) const;

	void MeshFromVerticesAndTriangles(Mesh& _result_mesh,
		const std::vector<std::array<double, 3>>& _vertices, const std::vector<std::array<int, 3>>& _triangles) const;
	//void CombineDisjointMeshes() const;

	// input _uv_points[i][j] is point u_i and v_j
	void MeshFromUVSurface(const std::vector<std::vector<Mesh::Point>>& _uv_points, Mesh& _surface_mesh);
	void MeshFromUVSurface_PeriodU(const std::vector<std::vector<Mesh::Point>>& _uv_points, Mesh& _surface_mesh);
	void MeshFromUVSurface_PeriodUV(const std::vector<std::vector<Mesh::Point>>& _uv_points, Mesh& _surface_mesh, bool _flip_normal = false);
	// 4 polyline enclose a watertight cuboid tube shape mesh, 
	// notice input polyline should have same point sequence
	// also, 4 polyline having correct order
	void CuboidTubeMeshFrom4Polylines_Period(const std::vector<Mesh::Point>& _polyline1,
		const std::vector<Mesh::Point>& _polyline2, const std::vector<Mesh::Point>& _polyline3, const std::vector<Mesh::Point>& _polyline4,
		Mesh& _result);
	// this is like last function, just input many polylines (pointers to the polylines)
	void CylinricalMeshFromPolylines_Period(const std::vector<const std::vector<Mesh::Point>*>& _polylines, Mesh& _result);

	void PointCloudMeshFromPoints(const std::vector<Mesh::Point>& _vertices, Mesh& _point_cloud) const;
	void PointCloudMeshFromPoints(const std::vector<Eigen::Vector3d>& _vertices, Mesh& _point_cloud) const;
	void PointCloudMeshToEigenVectorX(const Mesh& _point_cloud, std::vector<Eigen::VectorXd>& _points) const;
	void PointCloudMeshToPoints(const Mesh& _point_cloud, std::vector<Mesh::Point>& _points) const;

	// merge a mesh to main mesh
	// i.e. simply give all vertices and faces to the main mesh
	// the added vertex index in this vector same as the vertex id of merge mesh (handle is the main mesh handle)
	void MergeMesh(Mesh& _main, const Mesh& _merge, std::vector<Mesh::VertexHandle>& _added_vertex) const;

	//void PolylineResample(const std::vector<Mesh::Point>& _original, std::vector<Mesh::Point>& _resampled, int _resample_number) const;

private:
	// for mesh construction from uv surface 
	void MeshAddVertexFromUVPoints(const std::vector<std::vector<Mesh::Point>>& _uv_points, Mesh& _surface_mesh, 
		std::vector<std::vector<Mesh::VertexHandle>>& _uv_handle) const;
	// functions below provide different face normal directions and periodic options
	void MeshAddFaceFromUVPoints(Mesh& _surface_mesh, const std::vector<std::vector<Mesh::VertexHandle>>& _uv_handle) const;
	void MeshAddFaceFromUVPoints_PeriodicU(Mesh& _surface_mesh, const std::vector<std::vector<Mesh::VertexHandle>>& _uv_handle) const;
	void MeshAddFaceFromUVPoints_PeriodicUV(Mesh& _surface_mesh, const std::vector<std::vector<Mesh::VertexHandle>>& _uv_handle) const;
	void MeshAddFaceFromUVPoints_FlipNormal(Mesh& _surface_mesh, const std::vector<std::vector<Mesh::VertexHandle>>& _uv_handle) const;
	void MeshAddFaceFromUVPoints_PeriodicUV_FlipNormal(Mesh& _surface_mesh, const std::vector<std::vector<Mesh::VertexHandle>>& _uv_handle) const;

private:
	const Mesh* whole_mesh_;
	std::vector<int> face_area_; // indices of the faces to form a new mesh

	Mesh* sub_mesh_;
};
