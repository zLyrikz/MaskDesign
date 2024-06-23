#pragma once
#include "../MeshViewer/MeshDefinition.h"
#include <Eigen/Geometry>
#include <stack>


/*
* usage: 
1. initialize object --- AabbTree tree;
2. construct tree -- tree.constructTree(mesh);
3. call
*/


class AabbTree
{
public:
	class Box3
	{
	public:
		Mesh::Point min = Mesh::Point(1e40, 1e40, 1e40);
		Mesh::Point max = Mesh::Point(-1e40, -1e40, -1e40);

		Box3(){};
		Box3(Mesh::Point _min, Mesh::Point _max);
		//~Box3();
		//Box3(const Box3& other);

		bool isEmpty() const;
		bool intersectsBox(Box3& other);

		void expandByPoint(Mesh::Point p);
		void expandByOffset(double offset);
		void expandByFactor(double factor);
		Mesh::Point getCenter() const;
		Mesh::Point getSize() const;
		void setToCenter(Mesh::Point* target);
		void setToSize(Mesh::Point* target);
		void setToHalfSize(Mesh::Point* target);
		bool equals(Box3& other);

		Mesh::Point* getBoundById(unsigned int id);

		bool isInside(Mesh::Point& pt);
	};

	struct AabbNode
	{
		std::shared_ptr<AabbNode> left = nullptr;
		std::shared_ptr<AabbNode> right = nullptr;

		AabbTree* tree = nullptr;

		Box3 bbox;
		uint axis = 2;// split axis
		uint depth = 0;
		float splitPosition = 0.0f;// basically, bboxLeft.max[this->axis] = this->splitPosition; 

		std::vector<uint> primitiveIds = {};

		AabbNode() = default;
		//AabbNode(const AabbNode & other);
		AabbNode(const std::vector<uint>&primitiveIds, const Box3 & bbox, AabbTree * tree, uint depthLeft);
		~AabbNode() = default;

		void construct(const std::vector<uint>&primitiveIds, uint depthLeft);
		bool isALeaf();
		bool isALeafWithPrimitives();

		bool getPlaneBoxIntersection(Mesh::Point* normal, Mesh::Point* vert, Mesh::Point* boxMax);
		bool getTiangleBoxIntersection(Mesh::Point* T, Mesh::Point* boxCenter, Mesh::Point* boxHalfSize);

		// find an optimal split position using an "Adaptive Error-Bounded Heuristic" (Hunt, Mark, Stoll)
		double getSplitPosition(const std::vector<uint>&primitiveIds, std::vector<uint>&out_left, std::vector<uint>&out_right);
		// TODO: notice, some build error happen if uncomment this function
		//float getAdaptivelyResampledSplitPosition(const std::vector<uint>&primitiveIds, std::vector<uint>&out_left, std::vector<uint>&out_right);
		float getCostEstimate(float splitPos, uint nLeft, uint nRight);
		bool hasEnoughBranching(size_t nLeftPrims, size_t nRightPrims, size_t nPrims);
		void filterPrimitives();  // filters only primitives which actually intersect leaf box

		//void applyMatrix(Matrix4 & m);

		//bool useIntrinscs() const
		//{
		//	return tree->useIntrinsics;
		//}
	};

public:
	AabbTree();
	~AabbTree() = default;

	void constructTree(const Mesh * _tri_mesh);// tree for triangle mesh

	// TODO: now it always computes signed distance, which means [_is_signed] is of no use
	// note the sign is a bit inaccurate when the point is close too the mesh
	// note, time not counted!
	double Point2PrimitivesDistance(Mesh::Point& _point, bool _is_signed = true, long int* _time = nullptr) const;// find min distance
	// signed distance
	double PointCloud2PrimitivesMaxDistance(Mesh& _point_cloud) const;// max{ |distance(p,tree)| for all p in point cloud} sign is same as the distance taken
	void ClosestPointOnPrimitives(Mesh::Point& _point, Mesh::Point& _closest_point) const;// find the closest point on the primitives to a query point

	// only for triangle mesh
	// return first intersection if there is more than one (we can prove the first hit AABB leaf is the closest intersect leaf to the ray start point)
	// return true iff intersect; if not intersect, return the intersection point as the start point of the ray
	bool RayMeshIntersection(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, Mesh::Point& _intersect_point) const;
	// same as last function, but also output intersection face id
	bool RayMeshIntersection(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, Mesh::Point& _intersect_point, 
		uint& _face_id) const;
	// same as the function above, only return all the intersection points (also face id info version)
	// all intersect points are push_back to _intersect_point; if no intersect, then no new element added to _intersect_point
	void RayMeshIntersection(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, std::vector<Mesh::Point>& _intersect_point) const;
	void RayMeshIntersection(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, std::vector<Mesh::Point>& _intersect_point, std::vector<uint>& _face_id) const;

	// this is like the ray mesh intersection, only check the both sides of the ray
	// output the closest point on the mesh to the ray start point 
	bool LineMeshIntersection(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, Mesh::Point& _intersect_point, uint& _face_id) const;
	// this only restrict the intersection within the line segment; if more than one intersection, we return the closest one to the line start; include the t value
	// return -1 if no intersection even if for the ray; (here the ray is extending the line ending point along the line segment)
	// return 0 if ray intersect but already extending outside of the line end point
	// return 1 if line segment intersect the mesh
	int LineSegmentMeshIntersection(const Mesh::Point& _line_start, const Mesh::Point& _line_end, Mesh::Point& _intersect_point, uint& _face_id, double& _t) const;
	// this is to check line segment mesh intersection; 
	// but in case intersect out of line segment, return the closest intersect point to the line segment on the line, _t is reference distance to the line start
	// 
	// return -1 if no intersect for the whole line
	// return 0 if line intersect but already extending outside of the line segment
	// return 1 if line segment intersect the mesh
	int LineMeshIntersection_ClosestToASegment(const Mesh::Point& _line_start, const Mesh::Point& _line_end, Mesh::Point& _intersect_point, uint& _face_id, double& _t) const;

	AabbNode* getRootNode() const;
	//void getAllNodes(std::vector<AabbTree::AabbNode>* buffer);// Why does this cause link error?
	void getBoxes(std::vector<std::pair<Mesh::Point, Mesh::Point>>& _nodes_boxes) const;
	void getLeafBoxes(std::vector<std::pair<Mesh::Point, Mesh::Point>>& _leaf_boxes) const;

	uint getMaxDepth() const;
	const Mesh* getPrimitives() const;


public:
	// utility function

	// if point inside box, distance = 0
	// see &5.1.3 http://staff.ustc.edu.cn/~zhuang/hpg/books/Morgan_Kaufmann.Real-Time%20Collision%20Detection(2005).pdf
	double Point2AabbSquareDistance(Mesh::Point& _point, const Box3& _box) const;
	double Point2LeafPrimitivesDistance(const Eigen::Vector3d& _point, const AabbNode* _leaf, bool _is_signed = true) const;
	void ClosestPointOnLeafPrimitives(const Eigen::Vector3d& _point, const AabbNode* _leaf, Mesh::Point& _closest_point) const;

	// see &5.3.3 P179  http://staff.ustc.edu.cn/~zhuang/hpg/books/Morgan_Kaufmann.Real-Time%20Collision%20Detection(2005).pdf
	bool IntersectRayAABB(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, const Box3& _AABB, /*float& _t_min,*/ Mesh::Point& _intersect_point) const;// coding note: a const function cannot call a non-const function
	bool RayLeafPrimitiveIntersection(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, const AabbNode* _leaf, Mesh::Point& _intersect_point) const;
	bool RayLeafPrimitiveIntersection(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, const AabbNode* _leaf, Mesh::Point& _intersect_point,
		uint& _face_id) const;// same as last function, but also output intersection face id
	void RayLeafPrimitiveIntersection(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, const AabbNode* _leaf,
		std::vector<Mesh::Point>& _intersect_point, std::vector<uint>& _face_id) const;
	void RayLeafPrimitiveIntersection(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, const AabbNode* _leaf,
		std::vector<Mesh::Point>& _intersect_point, std::vector<float>& _t) const;
	// the t parameter has property: intersect_point = _ray_start + t*_ray_direction
	void RayLeafPrimitiveIntersection(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, const AabbNode* _leaf,
		std::vector<Mesh::Point>& _intersect_point, std::vector<uint>& _face_id, std::vector<float>& _t) const;
	// return -1 if no intersection even if for the ray; (here the ray is extending the line ending point along the line segment)
	// return 0 if ray intersect but already extending outside of the line end point
	// return 1 if line segment intersect the mesh
	// if ray intersects, output _t as the closest t value
	int LineSegmentLeafPrimitiveIntersection(const Mesh::Point& _line_start, const Mesh::Point& _line_end, const AabbNode* _leaf, Mesh::Point& _intersect_point, uint& _face_id, double& _t) const;

	// see &5.3.6 P190  http://staff.ustc.edu.cn/~zhuang/hpg/books/Morgan_Kaufmann.Real-Time%20Collision%20Detection(2005).pdf
	// Given ray and triangle abc, returns whether ray intersects
	// triangle and if so, also returns the barycentric coordinates (u,v,w)
	// of the intersection point
	// notice, for the floating point error, we allow a bit wider range of the triangle to regard a point inside the triangle
	bool IntersectRayTriangle(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, 
		const Mesh::Point& _triangle_a, const Mesh::Point& _triangle_b, const Mesh::Point& _triangle_c,
		float& _u, float& _v, float& _w, float& _t) const;

private:
	// some core part shared by similar functions
	
	// for ray mesh intersection
	// traverse a given aabb stack until a leaf node is found; not change the _leaf_node value is no leaf is hit
	void TraverseAabbStack(/*input*/ std::stack<std::pair<Mesh::Point*, AabbNode*>>& _aabb_stack, 
		const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction,
		/*output*/AabbNode*& _leaf_node) const;
	void ClearAabbStack(std::stack<std::pair<Mesh::Point*, AabbNode*>>& _aabb_stack) const;
	// output closest point index
	int GetClosestOne(const std::vector<Mesh::Point>& _intersect_points, const Mesh::Point& _reference_point) const;
	// for RayLeafPrimitiveIntersection
	//

	// may contain duplicate points
	void RayMeshIntersection_Raw(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, std::vector<Mesh::Point>& _intersect_point, std::vector<uint>& _face_id) const;

private:
	Box3 bbox;

	uint depth = 0;
	uint MAX_DEPTH = INT_MAX/*20*/;// root depth = 0; leaf depth <= MAX_DEPTH

	//std::vector<Primitive> primitives = {};
	std::shared_ptr<AabbNode> root = nullptr;
	const Mesh* mesh_ = nullptr;
};

