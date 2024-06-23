#include "AabbTree.h"
#include "../MeshProcessing/MeshInfo.h"
#include "../MeshProcessing/PointTriangleDistance.h"
#include "../Utility/VectorOperations.h"

//#include <Windows.h>
#include <iomanip>
#include <iostream>
#include <Eigen/Core>

using std::cout;
using std::endl;

AabbTree::AabbTree()
{
}

void AabbTree::constructTree(const Mesh* _tri_mesh)
{
	mesh_ = _tri_mesh;

	MeshInfo find_bounding_box(_tri_mesh);
	find_bounding_box.getBoundingBox(bbox.min, bbox.max, true);
	//cout << "bounding box:\n min:\n" << bbox.min << endl << "max:\n" << bbox.max << endl;

	// generate index array to primitives
	std::vector<uint> primitiveIds(mesh_->n_faces());
	std::iota(primitiveIds.begin(), primitiveIds.end(), 0);

	// generate root and all that follow 
	
	this->root = std::make_shared<AabbNode>(primitiveIds, bbox, this, MAX_DEPTH);
}

double AabbTree::Point2PrimitivesDistance(Mesh::Point& _point, bool _is_signed, long int* _time) const
{
	// traverse the tree from root node
	// best-first search: always explore the node in stack with the minimum distance
	// a property of this search: 
	// Property: For any node that has a larger distance to the point than optimal, it will not be explored.

	// record computing time on leaf nodes to test efficiency
	if (_time)
	{
		*_time = 0;
	}

	//int num_triangles = 0;
	//int num_leaf = 0;

	double distance = FLT_MAX;
	Eigen::Vector3d point_eigen = Eigen::Map<Eigen::Vector3d>(_point.data());

	double root_square_distance = Point2AabbSquareDistance(_point, bbox);
	std::map<float, AabbNode*> nodes_to_serach = { {root_square_distance, root.get()} };// initial root node

	while (!nodes_to_serach.empty())
	{		
		if (nodes_to_serach.begin()->first >= distance * distance) // rest of the nodes have larger distances  
		{
			//cout << "num_triangles evaluated " << num_triangles << endl;
			//cout << "num_leaf " << num_leaf << endl;
			return distance;
		}
		else
		{
			// take min distance node out
			AabbNode* current_node = nodes_to_serach.begin()->second;
			nodes_to_serach.erase(nodes_to_serach.begin());

			if (current_node->isALeaf())
			{
				// note, time not counted!
				//long int time0 = GetTickCount64();

				double distance2node = Point2LeafPrimitivesDistance(point_eigen, current_node, _is_signed);
				//cout << "distance2node = " << distance2node << endl;

				//long int time1 = GetTickCount64();
				if (_time)
				{
					//*_time += time1 - time0;
				}
				//num_triangles += current_node->primitiveIds.size();
				//num_leaf += 1;

				if (abs(distance2node) < abs(distance))
				{
					distance = distance2node;
				}
			}
			else
			{
				AabbNode* left = current_node->left.get();
				AabbNode* right = current_node->right.get();
				if (left != nullptr)
				{
					double left_square_distance = Point2AabbSquareDistance(_point, left->bbox);
					if (left_square_distance < distance * distance)
					{
						nodes_to_serach.insert({ left_square_distance, left });
					}
				}
				if (right != nullptr)
				{
					double right_square_distance = Point2AabbSquareDistance(_point, right->bbox);
					if (right_square_distance < distance * distance)
					{
						nodes_to_serach.insert({ right_square_distance, right });
					}
				}
			}
		}
	}

	//cout << "num_triangles evaluated " << num_triangles << endl;
	//cout << "num_leaf " << num_leaf << endl;

	return distance;
}

double AabbTree::PointCloud2PrimitivesMaxDistance(Mesh& _point_cloud) const
{
	double result = 0.0;
	for (Mesh::VertexHandle iVertex : _point_cloud.vertices())
	{
		Mesh::Point point = _point_cloud.point(iVertex);
		double current_distance = Point2PrimitivesDistance(point);
		if (abs(current_distance) > abs(result))
		{
			result = current_distance;
		}
	}

	return result;
}

void AabbTree::ClosestPointOnPrimitives(Mesh::Point& _point, Mesh::Point& _closest_point) const
{
	double distance = FLT_MAX;
	Eigen::Vector3d point_eigen = Eigen::Map<Eigen::Vector3d>(_point.data());
	Mesh::Point temp_closest_point;

	double root_square_distance = Point2AabbSquareDistance(_point, bbox);
	std::map<float, AabbNode*> nodes_to_serach = { {root_square_distance, root.get()} };// initial root node


	while (!nodes_to_serach.empty())
	{
		if (nodes_to_serach.begin()->first >= distance * distance) // rest of the nodes have larger distances  
		{
			return;
		}
		else
		{

			// take min distance node out
			AabbNode* current_node = nodes_to_serach.begin()->second;
			nodes_to_serach.erase(nodes_to_serach.begin());

			if (current_node->isALeaf())
			{
				
				ClosestPointOnLeafPrimitives(point_eigen, current_node, temp_closest_point);
				double distance2node = (_point - temp_closest_point).norm();
				//cout << "distance2node = " << distance2node << endl;

				if (distance2node < distance)
				{
					distance = distance2node;
					_closest_point = temp_closest_point;
				}
			}
			else
			{
				AabbNode* left = current_node->left.get();
				AabbNode* right = current_node->right.get();
				if (left != nullptr)
				{
					double left_square_distance = Point2AabbSquareDistance(_point, left->bbox);
					if (left_square_distance < distance * distance)
					{
						nodes_to_serach.insert({ left_square_distance, left });
					}
				}
				if (right != nullptr)
				{
					double right_square_distance = Point2AabbSquareDistance(_point, right->bbox);
					if (right_square_distance < distance * distance)
					{
						nodes_to_serach.insert({ right_square_distance, right });
					}
				}
			}
		}
	}
}

bool AabbTree::RayMeshIntersection(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, Mesh::Point& _intersect_point) const
{
	// initialize output value
	bool is_intersect = false;
	_intersect_point = _ray_start;

	// initialize a stack to begin the search

	// the point in stack is intersction of the node box with the ray, 
	// it'll be null if the intersection hasn't been tested.
	std::stack<std::pair<Mesh::Point*, AabbNode*>> aabb_stack; 
	aabb_stack.push(std::make_pair(nullptr, root.get()));

	do {
		AabbNode* first_hit_leaf = nullptr;
		TraverseAabbStack(aabb_stack, _ray_start, _ray_direction, first_hit_leaf);
		if (first_hit_leaf != nullptr)
		{
			assert(first_hit_leaf->isALeaf());
			is_intersect = RayLeafPrimitiveIntersection(_ray_start, _ray_direction, first_hit_leaf, _intersect_point);

			if (is_intersect)
			{
				//// before return, delete all the point pointer in stack
				ClearAabbStack(aabb_stack);

			}
		}
		else
		{
			assert(aabb_stack.empty() && !is_intersect);
		}
	} while (!is_intersect && !aabb_stack.empty());// didn't find intersection while there is still something in the stack
	
	return is_intersect;
}

bool AabbTree::RayMeshIntersection(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, Mesh::Point& _intersect_point, uint& _face_id) const
{
	// initialize output value
	bool is_intersect = false;
	_intersect_point = _ray_start;

	// initialize a stack to begin the search

	// the point in stack is intersction of the node box with the ray, 
	// it'll be null if the intersection hasn't been tested.
	std::stack<std::pair<Mesh::Point*, AabbNode*>> aabb_stack;
	aabb_stack.push(std::make_pair(nullptr, root.get()));

	do {
		AabbNode* first_hit_leaf = nullptr;
		TraverseAabbStack(aabb_stack, _ray_start, _ray_direction, first_hit_leaf);
		
		if (first_hit_leaf != nullptr)
		{
			assert(first_hit_leaf->isALeaf());

			// see if we will have an intersection
			is_intersect = RayLeafPrimitiveIntersection(_ray_start, _ray_direction, first_hit_leaf, _intersect_point, _face_id);
			if (is_intersect)
			{
				//// before return, delete all the point pointer in stack
				ClearAabbStack(aabb_stack);

			}
		}
		else
		{
			assert(aabb_stack.empty() && !is_intersect);
		}
	} while (!is_intersect && !aabb_stack.empty());// didn't find intersection while there is still something in the stack

	return is_intersect;
}

void AabbTree::RayMeshIntersection(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, std::vector<Mesh::Point>& _intersect_point) const
{
	// face_id not used
	std::vector<uint> face_id;
	std::vector<Mesh::Point> intersect_point;

	RayMeshIntersection_Raw(_ray_start, _ray_direction, intersect_point, face_id);

	// unique duplicate elements
	VectorOperations unique_it;
	std::vector<int> kept_element_id;//not used
	unique_it.UniquePoint(intersect_point, _intersect_point, kept_element_id, 1e-4);
}

void AabbTree::RayMeshIntersection(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, std::vector<Mesh::Point>& _intersect_point, std::vector<uint>& _face_id) const
{
	std::vector<Mesh::Point> intersect_point;
	std::vector<uint> face_id;
	RayMeshIntersection_Raw(_ray_start, _ray_direction, intersect_point, face_id);

	// unique duplicate elements
	VectorOperations unique_it;
	std::vector<int> kept_element_id;
	unique_it.UniquePoint(intersect_point, _intersect_point, kept_element_id, 1e-4);

	for (int i = 0; i < kept_element_id.size(); ++i)
	{
		_face_id.push_back(face_id[kept_element_id[i]]);
	}
}

bool AabbTree::LineMeshIntersection(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, Mesh::Point& _intersect_point, uint& _face_id) const
{
	Mesh::Point point1;
	uint face_id1;
	bool is_intersect1 = RayMeshIntersection(_ray_start, _ray_direction, point1, face_id1);
	Mesh::Point point2;
	uint face_id2;
	bool is_intersect2 = RayMeshIntersection(_ray_start, -_ray_direction, point2, face_id2);

	if (is_intersect1 && is_intersect2)// intersect on both sides
	{
		double distance1 = (point1 - _ray_start).norm();
		double distance2 = (point2 - _ray_start).norm();
		if (distance1 < distance2)
		{
			_face_id = face_id1;
			_intersect_point = point1;
		}
		else
		{
			_face_id = face_id2;
			_intersect_point = point2;
		}
		return true;
	}
	else if (is_intersect1)
	{
		_face_id = face_id1;
		_intersect_point = point1;
		return true;
	}
	else if (is_intersect2)
	{
		_face_id = face_id2;
		_intersect_point = point2;
		return true;
	}
	else
	{
		_intersect_point = _ray_start;
		return false;
	}
}

int AabbTree::LineSegmentMeshIntersection(const Mesh::Point& _line_start, const Mesh::Point& _line_end, Mesh::Point& _intersect_point, uint& _face_id, double& _t) const
{
	// the frame work is the same with the ray mesh intersection functions

	// initialize output value
	int intersect_state = 0;

	// initialize a stack to begin the search
	// the point in stack is intersction of the node box with the ray, 
	// it'll be null if the intersection hasn't been tested.
	std::stack<std::pair<Mesh::Point*, AabbNode*>> aabb_stack;
	aabb_stack.push(std::make_pair(nullptr, root.get()));
	Mesh::Point line_direction = _line_end - _line_start;
	do {
		AabbNode* first_hit_leaf = nullptr;
		TraverseAabbStack(aabb_stack, _line_start, line_direction, first_hit_leaf);

		if (first_hit_leaf != nullptr)
		{
			assert(first_hit_leaf->isALeaf());

			// see if we will have an intersection
			intersect_state = LineSegmentLeafPrimitiveIntersection(_line_start, _line_end, first_hit_leaf, _intersect_point, _face_id, _t);
			if (intersect_state != -1)// find intersection for the ray, will clear the aabb stack
			{
				//if (intersect_state == 1)
				//{
				//	is_intersect = true;
				//}
				//else 
				//{
				//	assert(intersect_state == 0);
				//	
				//}

				// no need to find other aabb, because the first one is the closest one
				//// before return, delete all the point pointer in stack
				ClearAabbStack(aabb_stack);
			}
		}
		else
		{
			assert(aabb_stack.empty() && intersect_state != 1);
		}
	} while (intersect_state != 1 && !aabb_stack.empty());// didn't find intersection while there is still something in the stack

	return intersect_state;

}

int AabbTree::LineMeshIntersection_ClosestToASegment(const Mesh::Point& _line_start, const Mesh::Point& _line_end, Mesh::Point& _intersect_point, uint& _face_id, double& _t) const
{
	int final_state = -1;//initialize

	// check the ray: start->end
	int intersect_state1 = LineSegmentMeshIntersection(_line_start, _line_end, _intersect_point, _face_id, _t);
	if (intersect_state1 == 1)
	{
		final_state = 1;
	}
	else
	{
		// check the ray: end->start
		double t2 = 0.0;
		Mesh::Point intersect_point2;
		uint face_id2 = 0;
		int intersect_state2 = LineSegmentMeshIntersection(_line_end, _line_start, intersect_point2, face_id2, t2);

		if (intersect_state2 == 1)
		{
			// this case is strange, it should only happen due to floating point error
			_intersect_point = std::move(intersect_point2);
			_face_id = face_id2;
			_t = 1 - t2;// this is to make the parameter based on the line_start point;
			final_state = 1;
			cout << "[WARINING] AabbTree::LineMeshIntersection_ClosestToASegment strange case" << endl;
		}
		else if (intersect_state1 == -1 && intersect_state2 == -1)
		{
			final_state = -1;
		}
		else if (intersect_state1 == -1 && intersect_state2 == 0)
		{
			_intersect_point = std::move(intersect_point2);
			_face_id = face_id2;
			_t = 1 - t2;// this is to make the parameter based on the line_start point;
			final_state = 0;
		}
		else if (intersect_state1 == 0 && intersect_state2 == 0)
		{
			if (t2 < _t)// compare which is closer
			{
				_intersect_point = std::move(intersect_point2);
				_face_id = face_id2;
				_t = 1 - t2;// this is to make the parameter based on the line_start point;
			}
			final_state = 0;
		}
		else
		{
			assert(intersect_state1 == 0 && intersect_state2 == -1);

			final_state = 0;
		}
	}

	return final_state;
}

AabbTree::AabbNode* AabbTree::getRootNode() const
{
	return root.get();
}

//void AabbTree::getAllNodes(std::vector<AabbTree::AabbNode>* buffer)
//{
//	std::stack<AabbNode*> nodeStack;
//	nodeStack.push(this->root.get());
//
//	while (!nodeStack.empty()) {
//		AabbNode* item = nodeStack.top();
//		nodeStack.pop();
//
//		buffer->push_back(*item);
//
//		if (item->left) {
//			nodeStack.push(item->left.get());
//		}
//		if (item->right) {
//			nodeStack.push(item->right.get());
//		}
//	}
//}

void AabbTree::getBoxes(std::vector<std::pair<Mesh::Point, Mesh::Point>>& _nodes_boxes) const
{
	std::stack<AabbNode*> nodeStack;
	nodeStack.push(root.get());

	while (!nodeStack.empty())
	{
		AabbNode* top_node = nodeStack.top();
		nodeStack.pop();
		
		_nodes_boxes.push_back(std::pair<Mesh::Point, Mesh::Point >(top_node->bbox.min, top_node->bbox.max));

		if (top_node->left)
		{
			nodeStack.push(top_node->left.get());
		}
		if (top_node->right)
		{
			nodeStack.push(top_node->right.get());
		}
	}
}

void AabbTree::getLeafBoxes(std::vector<std::pair<Mesh::Point, Mesh::Point>>& _leaf_boxes) const
{
	std::stack<AabbNode*> nodeStack;
	nodeStack.push(root.get());

	while (!nodeStack.empty())
	{
		AabbNode* top_node = nodeStack.top();
		nodeStack.pop();

		if (top_node->isALeaf())
		{
			_leaf_boxes.push_back(std::pair<Mesh::Point, Mesh::Point >(top_node->bbox.min, top_node->bbox.max));
		}

		if (top_node->left)
		{
			nodeStack.push(top_node->left.get());
		}
		if (top_node->right)
		{
			nodeStack.push(top_node->right.get());
		}
	}
}

uint AabbTree::getMaxDepth() const
{
	return MAX_DEPTH;
}

const Mesh* AabbTree::getPrimitives() const
{
	return mesh_;
}

double AabbTree::Point2AabbSquareDistance(Mesh::Point& _point, const Box3& _box) const
{
	double distance = 0.0;

	for (int iCoordinate = 0; iCoordinate < 3; ++iCoordinate)
	{
		if (_point[iCoordinate] < _box.min[iCoordinate])
		{
			distance += (_box.min[iCoordinate] - _point[iCoordinate]) * 
						(_box.min[iCoordinate] - _point[iCoordinate]);
		}
		else if (_point[iCoordinate] > _box.max[iCoordinate])
		{
			distance += (_point[iCoordinate] - _box.max[iCoordinate]) *
						(_point[iCoordinate] - _box.max[iCoordinate]);
		}
	}

	return distance;
}

double AabbTree::Point2LeafPrimitivesDistance(const Eigen::Vector3d& _point, const AabbNode* _leaf, bool _is_signed) const
{
	float distance = FLT_MAX;
	PointTriangleDistance find_distance;
	Eigen::Vector3f point = _point.cast<float>();

	// get distance to triangles in this leaf
	for (int iPrimitive = 0; iPrimitive < _leaf->primitiveIds.size(); ++iPrimitive)
	{
		// get vertices of this face
		Mesh::FaceHandle face = mesh_->face_handle(_leaf->primitiveIds[iPrimitive]);
		Mesh::Normal normal = mesh_->normal(face);
		Eigen::Vector3f normal_eigen = Eigen::Map<Eigen::Vector3d>(normal.data()).cast<float>();
		std::vector<Eigen::Vector3f> triangle;
		triangle.reserve(3);
		for (auto iVertexTriangle = mesh_->cfv_ccwiter(face); iVertexTriangle.is_valid(); ++iVertexTriangle)
		{
			Mesh::Point vertex = mesh_->point(*iVertexTriangle);
			Eigen::Vector3f vertex_eigen = Eigen::Map<Eigen::Vector3d>(vertex.data()).cast<float>();
			triangle.push_back(vertex_eigen);
		}

		/////////////////////////////////////////////////////////
		////for debug check if normal is accurate
		//Eigen::Vector3d recompute_normal = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]).normalized().cast<double>();
		//Eigen::Vector3d openmesh_normal = normal_eigen.cast<double>();
		//Eigen::Vector3d difference_normal = openmesh_normal-recompute_normal;
		//if (abs(difference_normal[0]) > 1e-5 || abs(difference_normal[1]) > 1e-5 || abs(difference_normal[2]) > 1e-5)
		//{
		//	cout << "!!!!!!!!!!!!!!!!!!!!!!NORMAL IN CONSISTENT!!!!!!!!!!!!!!!!!!!!!!" << difference_normal << endl;
		//}
		/////////////////////////////////////////////////////////


		float current_distance =
			find_distance.Point2Triangle(point, triangle[0], triangle[1], triangle[2], normal_eigen);


		if (fabs(current_distance) < fabs(distance))
		{
			distance = current_distance;
		}
	}


	return distance;
}

void AabbTree::ClosestPointOnLeafPrimitives(const Eigen::Vector3d& _point, const AabbNode* _leaf, Mesh::Point& _closest_point) const
{
	float distance = FLT_MAX;
	PointTriangleDistance find_distance;
	Eigen::Vector3f point = _point.cast<float>();
	Eigen::Vector3f closest_point;

	// get distance to triangles in this leaf
	for (int iPrimitive = 0; iPrimitive < _leaf->primitiveIds.size(); ++iPrimitive)
	{
		// get vertices of this face
		Mesh::FaceHandle face = mesh_->face_handle(_leaf->primitiveIds[iPrimitive]);
		Mesh::Normal normal = mesh_->normal(face);
		Eigen::Vector3f normal_eigen = Eigen::Map<Eigen::Vector3d>(normal.data()).cast<float>();
		std::vector<Eigen::Vector3f> triangle;
		triangle.reserve(3);
		for (auto iVertexTriangle = mesh_->cfv_ccwiter(face); iVertexTriangle.is_valid(); ++iVertexTriangle)
		{
			Mesh::Point vertex = mesh_->point(*iVertexTriangle);
			Eigen::Vector3f vertex_eigen = Eigen::Map<Eigen::Vector3d>(vertex.data()).cast<float>();
			triangle.push_back(vertex_eigen);
		}


		find_distance.ClosestPointOnTriangle(point, triangle[0], triangle[1], triangle[2], closest_point);
		float current_distance = (point - closest_point).norm();

		if (current_distance < distance)
		{
			distance = current_distance;
			_closest_point[0] = closest_point[0];
			_closest_point[1] = closest_point[1];
			_closest_point[2] = closest_point[2];
		}
	}
}

bool AabbTree::IntersectRayAABB(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, const Box3& _AABB, Mesh::Point& _intersect_point) const
{
	float tmin = 0.0f; // set to -FLT_MAX to get first hit on line
	float tmax = FLT_MAX; // set to max distance ray can travel (for segment)

	// For all three slabs(the space between a pair of parallel planes)
	for (int i = 0; i < 3; i++) {
		if (abs(_ray_direction[i]) < 1e-10) {
			// Ray is parallel to slab. No hit if origin not within slab
			if (_ray_start[i] < _AABB.min[i] || _ray_start[i] > _AABB.max[i]) return false;
		}
		else {
			// Compute intersection t value of ray with near and far plane of slab
			float ood = 1.0f / _ray_direction[i];
			float t1 = (_AABB.min[i] - _ray_start[i]) * ood;
			float t2 = (_AABB.max[i] - _ray_start[i]) * ood;
			// Make t1 be intersection with near plane, t2 with far plane; i.e. make t1 < t2
			if (t1 > t2)
			{
				// swap t1 and t2
				float t1_backup = t1;
				t1 = t2;
				t2 = t1_backup;
			}
			// Compute the intersection of slab intersection intervals
			if (t1 > tmin) tmin = t1;
			if (t2 < tmax) tmax = t2;// notice here, the original code on the book is [if (t2 > tmax)]
			// Exit with no collision as soon as slab intersection becomes empty
			if (tmin > tmax) return false;
		}
	}
	// Ray intersects all 3 slabs. Return point (q) and intersection t value (tmin)
	_intersect_point = _ray_start + _ray_direction * tmin;
	return true;

}

bool AabbTree::RayLeafPrimitiveIntersection(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, 
	const AabbNode* _leaf, Mesh::Point& _intersect_point) const
{
	std::vector<Mesh::Point> intersect_points;
	std::vector<uint> face_id;// not used later; only for input to the function below
	RayLeafPrimitiveIntersection(_ray_start, _ray_direction, _leaf, intersect_points, face_id);

	bool is_intersect = false;// initialize output intersection state
	if (intersect_points.size() > 0)
	{
		is_intersect = true;

		int closest_idx = GetClosestOne(intersect_points, _ray_start);
		_intersect_point = intersect_points[closest_idx];

	}
	else
	{
		_intersect_point = _ray_start;
	}

	return is_intersect;
}

bool AabbTree::RayLeafPrimitiveIntersection(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, const AabbNode* _leaf, 
	Mesh::Point& _intersect_point, uint& _face_id) const
{
	std::vector<Mesh::Point> intersect_points;
	std::vector<uint> face_id;
	RayLeafPrimitiveIntersection(_ray_start, _ray_direction, _leaf, intersect_points, face_id);

	bool is_intersect = false;// initialize output intersection state
	if (intersect_points.size() > 0)
	{
		is_intersect = true;

		int closest_idx = GetClosestOne(intersect_points, _ray_start);
		_intersect_point = intersect_points[closest_idx];
		_face_id = face_id[closest_idx];

	}
	else
	{
		_intersect_point = _ray_start;
	}

	return is_intersect;
}

void AabbTree::RayLeafPrimitiveIntersection(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, const AabbNode* _leaf,
	std::vector<Mesh::Point>& _intersect_point, std::vector<uint>& _face_id) const
{
	std::vector<float> t;// t is not used
	RayLeafPrimitiveIntersection(_ray_start, _ray_direction, _leaf, _intersect_point, _face_id, t);
}

void AabbTree::RayLeafPrimitiveIntersection(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, const AabbNode* _leaf,
	std::vector<Mesh::Point>& _intersect_point, std::vector<float>& _t) const
{
	std::vector<uint> face_id;// face_id is not used
	RayLeafPrimitiveIntersection(_ray_start, _ray_direction, _leaf, _intersect_point, face_id, _t);
}


void AabbTree::RayLeafPrimitiveIntersection(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, const AabbNode* _leaf,
	std::vector<Mesh::Point>& _intersect_point, std::vector<uint>& _face_id, std::vector<float>& _t) const
{
	for (int iTriangle = 0; iTriangle < _leaf->primitiveIds.size(); ++iTriangle)
	{
		// get vertices of current triangle
		uint face_id = _leaf->primitiveIds[iTriangle];
		Mesh::FaceHandle triangle = mesh_->face_handle(face_id);
		std::vector<Mesh::Point> points;
		points.reserve(3);
		for (Mesh::VertexHandle iVertex : mesh_->fv_range(triangle))
		{
			points.push_back(mesh_->point(iVertex));
		}

		// test intersection
		float barycentric_u = 0.f;
		float barycentric_v = 0.f;
		float barycentric_w = 0.f;
		float ray_hit_parameter = 0.f;
		bool is_intersect = IntersectRayTriangle(_ray_start, _ray_direction,
			points[0], points[1], points[2],
			barycentric_u, barycentric_v, barycentric_w, ray_hit_parameter);

		if (is_intersect)
		{
			_face_id.push_back(face_id);
			_intersect_point.push_back(barycentric_u * points[0] + barycentric_v * points[1] + barycentric_w * points[2]);
			_t.push_back(ray_hit_parameter);
		}
	}
}

int AabbTree::LineSegmentLeafPrimitiveIntersection(const Mesh::Point& _line_start, const Mesh::Point& _line_end, const AabbNode* _leaf, Mesh::Point& _intersect_point, uint& _face_id, double& _t) const
{
	Mesh::Point line_direction = _line_end - _line_start;

	std::vector<Mesh::Point> intersect_points;
	std::vector<float> t;
	std::vector<uint> face_id;
	RayLeafPrimitiveIntersection(_line_start, line_direction, _leaf, intersect_points, face_id, t);

	int intersect_state = -1;// initialize output intersection state
	if (intersect_points.size() > 0)
	{
		intersect_state = 0;// the ray intersects the mesh, will check if it's inside the line segment

		// TODO: we can get the closest one by using the t values! smallest t is the one! this suitable for all other similar functions in this class
		//int closest_idx = GetClosestOne(intersect_points, _line_start);
		std::vector<float>::iterator closest_t = std::min_element(t.begin(), t.end());
		int closest_idx = std::distance(t.begin(), closest_t);
		_t = *closest_t;
		if (*closest_t > -1e-5 && *closest_t < 1.0 + 1e-5)// intersect within the line segment
		{
			intersect_state = 1;
		}
		_intersect_point = std::move(intersect_points[closest_idx]);// pass the value anyway
		_face_id = face_id[closest_idx];
	}

	return intersect_state;
}

bool AabbTree::IntersectRayTriangle(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction,
	const Mesh::Point& _triangle_a, const Mesh::Point& _triangle_b, const Mesh::Point& _triangle_c, 
	float& _u, float& _v, float& _w, float& _t) const
{
	//soving linear system:
	// (-ray_direction, ab, ac).dot(t,v,w) = ray_start - a
	// i.e. _ray_start + t* ray_direction - a = v*ab + w* ac (intersection of ray and traingle's supporting plane)
	// intersect point = u*a + v*b + w*c;

	const Eigen::Vector3f a = Eigen::Map<const Eigen::Vector3d>(_triangle_a.data()).cast<float>();
	const Eigen::Vector3f ray_start = Eigen::Map<const Eigen::Vector3d>(_ray_start.data()).cast<float>();
	const Eigen::Vector3f ray_direction = Eigen::Map<const Eigen::Vector3d>(_ray_direction.data()).cast<float>();
	const Eigen::Vector3f ab = Eigen::Map<const Eigen::Vector3d>((_triangle_b - _triangle_a).data()).cast<float>();
	const Eigen::Vector3f ac = Eigen::Map<const Eigen::Vector3d>((_triangle_c - _triangle_a).data()).cast<float>();

	Eigen::Matrix3f coefficient_matrix;
	coefficient_matrix << -ray_direction, ab, ac;

	Eigen::Vector3f t_v_w = coefficient_matrix.colPivHouseholderQr().solve(ray_start - a);
	_t = t_v_w[0];
	_v = t_v_w[1];
	_w = t_v_w[2];
	_u = 1 - _v - _w;


	if (_t < -1e-5) // ray travel away from triangle
	{
		return false;
	}
	if (_u >= -1e-5 && _v >= -1e-5 && _w >= -1e-5)/*(_u >= 0 && _v >= 0 && _w >= 0)*/
	{
		return true;
	}
	else
	{
		return false;
	}
}

void AabbTree::TraverseAabbStack(std::stack<std::pair<Mesh::Point*, AabbNode*>>& _aabb_stack, 
	const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, 
	AabbNode*& _leaf_node) const
{
	while (!_aabb_stack.empty())
	{
		// initialize the intersection information
		bool is_intersect_aabb = true;

		Mesh::Point* intersect_aabb = _aabb_stack.top().first;
		AabbNode* current_node = _aabb_stack.top().second;
		if (intersect_aabb == nullptr)
		{
			// check intersection
			intersect_aabb = new Mesh::Point();
			is_intersect_aabb = IntersectRayAABB(_ray_start, _ray_direction, current_node->bbox, *intersect_aabb);
		}
		_aabb_stack.pop();

		// deal with the current node
		if (is_intersect_aabb)
		{
			if (current_node->isALeaf())
			{
				delete intersect_aabb; // this pointer is not to be used
				intersect_aabb = nullptr;

				_leaf_node = current_node;
				return; // return upon hitting a leaf
			}
			else
			{
				// find out the intersection is on left or right split box
				if ((*intersect_aabb)[current_node->axis] < current_node->splitPosition)
				{
					// on the left box

					if (current_node->right)
					{
						_aabb_stack.push(std::make_pair(nullptr, current_node->right.get()));
					}
					if (current_node->left)
					{
						_aabb_stack.push(std::make_pair(intersect_aabb, current_node->left.get()));
					}
					else
					{
						delete intersect_aabb;
						intersect_aabb = nullptr;
					}
				}
				else
				{
					// on the right box

					if (current_node->left)
					{
						_aabb_stack.push(std::make_pair(nullptr, current_node->left.get()));
					}
					if (current_node->right)
					{
						_aabb_stack.push(std::make_pair(intersect_aabb, current_node->right.get()));
					}
					else
					{
						delete intersect_aabb;
						intersect_aabb = nullptr;
					}
				}
			}
		}
		else
		{
			delete intersect_aabb; // this pointer is not to be used
			intersect_aabb = nullptr;
		}
	}

}

void AabbTree::ClearAabbStack(std::stack<std::pair<Mesh::Point*, AabbNode*>>& _aabb_stack) const
{
	// before return, delete all the point pointer in stack
	while (!_aabb_stack.empty())
	{
		Mesh::Point* intersect_point = _aabb_stack.top().first;
		if (intersect_point)
		{
			delete intersect_point;
			intersect_point = nullptr;
		}
		_aabb_stack.pop();
	}
}

int AabbTree::GetClosestOne(const std::vector<Mesh::Point>& _intersect_points, const Mesh::Point& _reference_point) const
{
	// find the closest one
	int idx = 0;
	double min_distance = FLT_MAX;
	for (int iPoint = 0; iPoint < _intersect_points.size(); ++iPoint)
	{
		double distance = (_intersect_points[iPoint] - _reference_point).sqrnorm();// squared norm is enough
		if (distance < min_distance)
		{
			min_distance = distance;
			idx = iPoint;
		}
	}
	return idx;
}

void AabbTree::RayMeshIntersection_Raw(const Mesh::Point& _ray_start, const Mesh::Point& _ray_direction, std::vector<Mesh::Point>& _intersect_point, std::vector<uint>& _face_id) const
{
	// initialize a stack to begin the search
	// the point in stack is intersction of the node box with the ray, 
	// it'll be null if the intersection hasn't been tested.
	std::stack<std::pair<Mesh::Point*, AabbNode*>> aabb_stack;
	aabb_stack.push(std::make_pair(nullptr, root.get()));

	// find all the intersect points
	do {
		AabbNode* first_hit_leaf = nullptr;
		TraverseAabbStack(aabb_stack, _ray_start, _ray_direction, first_hit_leaf);
		if (first_hit_leaf != nullptr)
		{
			assert(first_hit_leaf->isALeaf());

			std::vector<Mesh::Point> intersect_points;
			std::vector<uint> face_id;
			RayLeafPrimitiveIntersection(_ray_start, _ray_direction, first_hit_leaf, intersect_points, face_id);

			for (int iIntersect = 0; iIntersect < intersect_points.size(); ++iIntersect)
			{
				_intersect_point.push_back(std::move(intersect_points[iIntersect]));
				_face_id.push_back(face_id[iIntersect]);
			}

		}
		else
		{
			assert(aabb_stack.empty());
		}
	} while (!aabb_stack.empty());
}


AabbTree::AabbNode::AabbNode(const std::vector<uint>& primitiveIds, const Box3& bbox, AabbTree* tree, uint depthLeft)
{
	this->tree = tree;
	this->bbox = bbox;
	this->depth = tree->getMaxDepth() - depthLeft;// root depth = 0; leaf depth <= MAX_DEPTH

	this->construct(primitiveIds, depthLeft);// (children nodes are constructed by the parent node (not themselves))
}

void AabbTree::AabbNode::construct(const std::vector<uint>& primitiveIds, uint depthLeft)
{
	if (depthLeft == 0 || primitiveIds.size() <= 2) {
		this->primitiveIds = primitiveIds;
		this->filterPrimitives();// what is this used for?
		return;
	}

	// choosing the longest split axis seems to be faster
	if ((bbox.max[0] - bbox.min[0]) > (bbox.max[1] - bbox.min[1]) && (bbox.max[0] - bbox.min[0]) > (bbox.max[2] - bbox.min[2])) {
		this->axis = 0;
	}
	else if ((bbox.max[1] - bbox.min[1]) > (bbox.max[2] - bbox.min[2])) {
		this->axis = 1;
	}
	else {
		this->axis = 2;
	}

	std::vector<uint> leftPrimitiveIds = {};
	std::vector<uint> rightPrimitiveIds = {};

	if (false) // TODO: try this different split strategy
	{
		
		//this->splitPosition = this->getAdaptivelyResampledSplitPosition(primitiveIds, leftPrimitiveIds, rightPrimitiveIds); // classify primitives		
	}
	else
	{
		this->splitPosition = this->getSplitPosition(primitiveIds, leftPrimitiveIds, rightPrimitiveIds); // classify primitives	
	}

	const bool stopBranching = !this->hasEnoughBranching(leftPrimitiveIds.size(), rightPrimitiveIds.size(), primitiveIds.size());
	if (stopBranching) {
		this->primitiveIds = primitiveIds;
		this->filterPrimitives();
		return;
	}

	if (!leftPrimitiveIds.empty()) {
		Box3 bboxLeft = this->bbox;
		bboxLeft.max[this->axis] = this->splitPosition;
		bboxLeft.expandByFactor(1.01);
		this->left = std::make_shared<AabbNode>(leftPrimitiveIds, bboxLeft, this->tree, depthLeft - 1);

		if (this->left->primitiveIds.empty() && this->left->left == nullptr && this->left->right == nullptr) {
			this->left = nullptr;
		}
	}

	if (!rightPrimitiveIds.empty()) {
		Box3 bboxRight = this->bbox;
		bboxRight.min[this->axis] = this->splitPosition;
		bboxRight.expandByFactor(1.01);
		this->right = std::make_shared<AabbNode>(rightPrimitiveIds, bboxRight, this->tree, depthLeft - 1);

		if (this->right->primitiveIds.empty() && this->right->left == nullptr && this->right->right == nullptr) {
			this->right = nullptr;
		}
	}
}

bool AabbTree::AabbNode::isALeaf()
{
	return (!this->left && !this->right);
}

bool AabbTree::AabbNode::isALeafWithPrimitives()
{
	return this->isALeaf() && !this->primitiveIds.empty();
}

#pragma region Helper macro functions
// ====== Helper macros for vectors (to increase speed) ============

#define CROSS(dest, v1, v2)						\
          dest[0] = v1[1] * v2[2] - v1[2] * v2[1];   \
          dest[1] = v1[2] * v2[0] - v1[0] * v2[2];	\
          dest[2] = v1[0] * v2[1] - v1[1] * v2[0];

#define DOT(v1, v2) (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])

#define SUB(dest, v1, v2)						\
          dest[0] = v1[0] - v2[0];					\
          dest[1] = v1[1] - v2[1];					\
          dest[2] = v1[2] - v2[2];

#define FINDMINMAX(x0, x1, x2, min, max)		\
		  min = max = x0;						\
		  if (x1 < min) min = x1;				\
		  if (x1 > max) max = x1;				\
		  if (x2 < min) min = x2;				\
		  if (x2 > max) max = x2;

// =================== Helper macros for axis tests ======================
// X-tests:

#define AXISTEST_X01(a, b, fa, fb)									   \
		p0 = a * v0[1] - b * v0[2];			       				       \
		p2 = a * v2[1] - b * v2[2];			       					   \
        if (p0 < p2) {min = p0; max = p2;} else {min = p2; max = p0;}  \
		rad = fa * (*boxHalfSize)[1] + fb * (*boxHalfSize)[2];		       \
		if (min > rad || max < -rad) return false;

#define AXISTEST_X2(a, b, fa, fb)									   \
		p0 = a * v0[1] - b * v0[2];									   \
		p1 = a * v1[1] - b * v1[2];		       						   \
        if (p0 < p1) {min = p0; max = p1;} else {min = p1; max = p0;}  \
		rad = fa * (*boxHalfSize)[1] + fb * (*boxHalfSize)[2];			   \
		if (min > rad || max < -rad) return false;

// Y-tests:

#define AXISTEST_Y02(a, b, fa, fb)									   \
		p0 = -a * v0[0] + b * v0[2];				   					   \
		p2 = -a * v2[0] + b * v2[2];				       				   \
		if (p0 < p2) {min = p0; max = p2;} else {min = p2; max = p0;}  \
		rad = fa * (*boxHalfSize)[0] + fb * (*boxHalfSize)[2];		       \
		if (min > rad || max <- rad) return false;

#define AXISTEST_Y1(a, b, fa, fb)							    	   \
		p0 = -a * v0[0] + b * v0[2];		      					       \
		p1 = -a * v1[0] + b * v1[2];					                   \
		if (p0 < p1) {min = p0; max = p1;} else {min = p1; max = p0;}  \
		rad = fa * (*boxHalfSize)[0] + fb * (*boxHalfSize)[2];			   \
		if (min > rad || max < -rad) return false;

// Z-tests:

#define AXISTEST_Z12(a, b, fa, fb)									   \
		p1 = a * v1[0] - b * v1[1];			 			               \
		p2 = a * v2[0] - b * v2[1];			       	                   \
		if (p2 < p1) {min = p2; max = p1;} else {min = p1; max = p2;}  \
		rad = fa * (*boxHalfSize)[0] + fb * (*boxHalfSize)[1];			   \
		if (min > rad || max < -rad) return false;

#define AXISTEST_Z0(a, b, fa, fb)									   \
		p0 = a * v0[0] - b * v0[1];									   \
		p1 = a * v1[0] - b * v1[1];									   \
		if (p0 < p1) {min = p0; max = p1;} else {min = p1; max = p0;}  \
		rad = fa * (*boxHalfSize)[0] + fb * (*boxHalfSize)[1];			   \
		if (min > rad || max < -rad) return false;
#pragma endregion 

bool AabbTree::AabbNode::getPlaneBoxIntersection(Mesh::Point* normal, Mesh::Point* vert, Mesh::Point* boxMax)
{
	int q;
	Mesh::Point vmin, vmax;
	double v;

	for (q = 0; q <= 2; q++) {
		v = (*vert)[q];

		if ((*normal)[q] > 0.0) {
			vmin[q] = -(*boxMax)[q] - v;
			vmax[q] = (*boxMax)[q] - v;
		}
		else {
			vmin[q] = (*boxMax)[q] - v;
			vmax[q] = -(*boxMax)[q] - v;
		}
	}

	if (DOT((*normal), vmin) > 0.0) return false;
	if (DOT((*normal), vmax) >= 0.0) return true;
	return false;
}

bool AabbTree::AabbNode::getTiangleBoxIntersection(Mesh::Point* T, Mesh::Point* boxCenter, Mesh::Point* boxHalfSize)
{
	// this is a copied funtion, I think it returns true if the triangle intersects the box
	// but what the logic behind this method? looks like in the end, it reduces to box and plane of the triangle intersection

	Mesh::Point v0, v1, v2;// local coordinates of triangle vertices
	double min, max, p0, p1, p2, rad, fex, fey, fez;
	Mesh::Point normal, e0, e1, e2;

	SUB(v0, (T[0]), (*boxCenter));
	SUB(v1, (T[1]), (*boxCenter));
	SUB(v2, (T[2]), (*boxCenter));

	// tri edges:
	SUB(e0, v1, v0);
	SUB(e1, v2, v1);
	SUB(e2, v0, v2);

	// 9 axis tests:
	fex = fabsf(e0[0]);
	fey = fabsf(e0[1]);
	fez = fabsf(e0[2]);

	AXISTEST_X01(e0[2], e0[1], fez, fey);
	AXISTEST_Y02(e0[2], e0[0], fez, fex);
	AXISTEST_Z12(e0[1], e0[0], fey, fex);

	fex = fabsf(e1[0]);
	fey = fabsf(e1[1]);
	fez = fabsf(e1[2]);

	AXISTEST_X01(e1[2], e1[1], fez, fey);
	AXISTEST_Y02(e1[2], e1[0], fez, fex);
	AXISTEST_Z0(e1[1], e1[0], fey, fex);

	fex = fabsf(e2[0]);
	fey = fabsf(e2[1]);
	fez = fabsf(e2[2]);

	AXISTEST_X2(e2[2], e2[1], fez, fey);
	AXISTEST_Y1(e2[2], e2[0], fez, fex);
	AXISTEST_Z12(e2[1], e2[0], fey, fex);

	// test for AABB overlap in x, y, and z:
	// test in x:
	FINDMINMAX(v0[0], v1[0], v2[0], min, max);
	if (min > (*boxHalfSize)[0] || max < -(*boxHalfSize)[0]) return false;

	// test in y:
	FINDMINMAX(v0[1], v1[1], v2[1], min, max);
	if (min > (*boxHalfSize)[1] || max < -(*boxHalfSize)[1]) return false;

	// test in z:
	FINDMINMAX(v0[2], v1[2], v2[2], min, max);
	if (min > (*boxHalfSize)[2] || max < -(*boxHalfSize)[2]) return false;

	// test if the box intersects the triangle plane  dot(normal, x) + d = 0

	CROSS(normal, e0, e1);

	if (!getPlaneBoxIntersection(&normal, &v0, boxHalfSize)) return false;

	return true;
}

double AabbTree::AabbNode::getSplitPosition(const std::vector<uint>& primitiveIds, std::vector<uint>& out_left, std::vector<uint>& out_right)
{
	const double tot_BoxArea = 2.0 * (this->bbox.max[0] - this->bbox.min[0]) +
		2.0 * (this->bbox.max[1] - this->bbox.min[1]) +
		2.0 * (this->bbox.max[2] - this->bbox.min[2]);
	const uint CUTS = 4;// cut sample number
	const uint N_primitives = primitiveIds.size();
	uint i, j;

	// === Stage 1: Initial sampling of C_L(x) and C_R(x) ========================================================

	// split interval limits:
	const double a = bbox.min[axis];
	const double b = bbox.max[axis];

	// splits:
	std::vector<double> bCutPos = std::vector<double>(CUTS + 2);
	for (i = 0; i <= CUTS + 1; i++) bCutPos[i] = a * (1.0 - ((double)i / (double)(CUTS + 1))) + b * ((double)i / (double)(CUTS + 1));

	// set C_L(x) = 0, C_R(x) = 0
	auto C_L = std::vector<uint>(CUTS + 2);
	auto C_R = std::vector<uint>(CUTS + 2);

	auto primMins = std::vector<double>(N_primitives);
	auto primMaxes = std::vector<double>(N_primitives);

	const Mesh* primitives = this->tree->getPrimitives();
	MeshInfo find_primitive_min_max(primitives);
	for (i = 0; i < N_primitives; i++) 
	{		
		primMins[i] = find_primitive_min_max.getFaceMinCoordinate(primitiveIds[i], axis);
		primMaxes[i] = find_primitive_min_max.getFaceMaxCoordinate(primitiveIds[i], axis);

		for (j = 1; j <= CUTS; j++) {
			C_L[j] += (primMins[i] < bCutPos[j] ? 1 : 0);
			C_R[j] += (primMaxes[i] > bCutPos[j] ? 1 : 0);
		}
	}
	C_L[5] = C_R[0] = N_primitives;
	std::reverse(C_R.begin(), C_R.end()); // store in reverse since C_R(x) is non-increasing

	// ===== Stage 2: Sample range [0, N_primitives] uniformly & count the number of samples within each segment ======

	auto S_L = std::vector<double>(CUTS + 1);
	auto S_R = std::vector<double>(CUTS + 1);

	double ran_s;

	for (i = 0; i < CUTS; i++) {
		ran_s = (double)(i + 1) / (double)(CUTS + 1) * N_primitives;

		for (j = 0; j <= CUTS; j++) {
			S_L[j] += (ran_s > C_L[j] && ran_s < C_L[j + 1] ? 1 : 0);
			S_R[j] += (ran_s > C_R[j] && ran_s < C_R[j + 1] ? 1 : 0);
		}
	}
	std::reverse(S_R.begin(), S_R.end());

	// ==== Stage 3: add more sampling positions to subdivided segments ===========================================

	auto all_splt_L = std::vector<double>(2 * CUTS);
	auto all_splt_R = std::vector<double>(2 * CUTS);
	double segLen = (double)(b - a) / (double)(CUTS + 1);
	uint nSeg_L = 0, nSeg_R = 0;

	for (i = 0; i <= CUTS; i++) {
		if (i > 0) {
			all_splt_L[nSeg_L++] = bCutPos[i];
			all_splt_R[nSeg_R++] = bCutPos[i];
		}
		for (j = 0; j < S_L[i]; j++) {
			all_splt_L[nSeg_L++] = bCutPos[i] + (j + 1) / (S_L[i] + 1) * segLen;
		}
		for (j = 0; j < S_R[i]; j++) {
			all_splt_R[nSeg_R++] = bCutPos[i] + (j + 1) / (S_R[i] + 1) * segLen;
		}
	}

	// Compute surface area heuristic SAH:
	// remaining two dimensions of the child box candidates
	const double boxDim0 = bbox.min[(this->axis + 1) % 3] - bbox.min[(this->axis + 1) % 3];
	const double boxDim1 = bbox.min[(this->axis + 2) % 3] - bbox.min[(this->axis + 2) % 3];

	auto SA_L = std::vector<double>(2 * CUTS);
	auto SA_R = std::vector<double>(2 * CUTS);

	// SA_L(x) = (boxDim_L(x) + boxDim0 + boxDim1) * 2.0 / tot_BoxArea
	// SA_R(x) = (boxDim_R(x) + boxDim0 + boxDim1) * 2.0 / tot_BoxArea

	for (i = 0; i < 2 * CUTS; i++) {
		SA_L[i] = ((all_splt_L[i] - a) + boxDim0 + boxDim1) * 2.0 / tot_BoxArea;
		SA_R[i] = ((b - all_splt_R[i]) + boxDim0 + boxDim1) * 2.0 / tot_BoxArea;
	}

	// ==== Stage 4: RESAMPLE C_L and C_R on all sample points & construct an approximation of cost(x) to minimize
	double min, max;
	auto cost = std::vector<double>(2 * CUTS);

	// cost(x) = C_L(x) * SA_L(x) + C_R(x) * SA_R(x):
	for (i = 0; i < N_primitives; i++) {
		min = primMins[i];
		max = primMaxes[i];

		for (j = 0; j < 2 * CUTS; j++) {
			cost[j] += (min < all_splt_L[j] ? 1 : 0) * SA_L[j] + (max > all_splt_R[j] ? 1 : 0) * SA_R[j];
		}
	}

	// ==== Stage 5: Minimize cost(x) & classify primitives  =====================================================

	float bestSplit = 0.5 * (a + b); // if this loop fails to initialize bestSplit, set it to middle
	double minCost = DBL_MAX;

	for (i = 1; i < 2 * CUTS; i++) {
		if (cost[i] < minCost) {
			minCost = cost[i];
			bestSplit = all_splt_L[i];
		}
	}

	// fill left and right arrays now that best split position is known:
	out_left.reserve(N_primitives);
	out_right.reserve(N_primitives);
	for (i = 0; i < N_primitives; i++) {
		min = primMins[i];
		max = primMaxes[i];

		if (min <= bestSplit) {
			out_left.push_back(primitiveIds[i]);
		}
		if (max >= bestSplit) {
			out_right.push_back(primitiveIds[i]);
		}
	}
	out_left.shrink_to_fit();
	out_right.shrink_to_fit();

	return bestSplit;
}

bool AabbTree::AabbNode::hasEnoughBranching(size_t nLeftPrims, size_t nRightPrims, size_t nPrims)
{
	return nLeftPrims + nRightPrims < 1.5f * nPrims;
}

void AabbTree::AabbNode::filterPrimitives()
{
	Mesh::Point bboxCenter = bbox.getCenter();
	Mesh::Point bboxHalfSize = 0.5f * bbox.getSize();
	std::vector<uint> newPrimitiveIds = {};
	for (uint i = 0; i < this->primitiveIds.size(); i++) 
	{
		const Mesh* mesh = this->tree->getPrimitives();
		Mesh::FaceHandle triangle = mesh->face_handle(this->primitiveIds[i]);
		std::vector<Mesh::Point> vertices;
		vertices.reserve(3);
		for (auto iVertex = mesh->cfv_ccwiter(triangle); iVertex.is_valid(); ++iVertex)
		{
			vertices.push_back(mesh->point(*iVertex));
		}

		if (getTiangleBoxIntersection(&(vertices[0]), &bboxCenter, &bboxHalfSize)) {
			newPrimitiveIds.push_back(this->primitiveIds[i]);
		}
	}
	//cout << "AabbTree::AabbNode::filterPrimitives(): input primitive number = " << this->primitiveIds.size()
	//	<< ", output filtered number = " << newPrimitiveIds.size() << endl;

	this->primitiveIds = newPrimitiveIds;
}

AabbTree::Box3::Box3(Mesh::Point _min, Mesh::Point _max)
{
	min = _min;
	max = _max;
}

bool AabbTree::Box3::isEmpty() const
{
	return (
		(this->min[0] > FLT_MAX && this->min[1] > FLT_MAX && this->min[2] > FLT_MAX) &&
		(this->max[0] < -FLT_MAX && this->max[1] < -FLT_MAX && this->max[2] < -FLT_MAX)
		); // notice min, max initialized to +-1e40;
}

void AabbTree::Box3::expandByFactor(double factor)
{
	Mesh::Point scale = this->getSize();
	Mesh::Point scaled = factor * scale;
	Mesh::Point offset = 0.5 * (scaled - scale);
	this->min = this->min - offset;
	this->max = this->max + offset;
}

Mesh::Point AabbTree::Box3::getCenter() const
{
	if (isEmpty()) {
		return Mesh::Point();
	}

	Mesh::Point center = 0.5 * (min + max);
	return center;
}

Mesh::Point AabbTree::Box3::getSize() const
{
	if (isEmpty()) {
		return Mesh::Point();
	}

	Mesh::Point dims = (max - min);
	return dims;
}

