#include "InterfaceDesign.h"
#include "../MaskInterface/CushionSurface.h"
#include "../MaskInterface/CurveDistanceToFace.h"
#include "../MaskInterface/ConnectFrontNBack.h"
#include "../MeshProcessing/MeshInfo.h"
#include "../MeshProcessing/TriangleMeshDistance.h"
#include "../MeshProcessing/MeshOffset.h"
#include "../Optimization/RailOptimizeInfo.h"
#include "../Optimization/NelderMeadRailOptimize.h"
#include "../Optimization/CrossSectionOptimizeInfo.h"
#include "../Optimization/CushionOptimize.h"
#include "../Utility/ReadRhino.h"
#include <windows.h>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <Discregrid/All>
#include <iostream>
using std::cout;
using std::endl;
InterfaceDesign::InterfaceDesign(const Mesh* _human_face, std::string _base_path): human_face_(_human_face)
{
	base_path_ = _base_path;
	data_path_ = base_path_ + "/data";
}

InterfaceDesign::~InterfaceDesign()
{
}

void InterfaceDesign::Design(Mesh& _cushion_surface, Mesh& _mask_interface, int _max_cushion_optimization_iterations,
	/*rail curve adjustment parameters*/double _objective, double _width, double _curvature, double _symmetry, double _align, double _angle,
	/*FEM simulation parameters*/ 	double _E_module, double _poisson_ratio, int _num_layer, double _thickness, 
	/*cushion optimization paremeters*/double _weight_average_pressure , double _force_distribution, double _area_distribution,
	double _convexity, double _weight_tetrahedral_inversion, double _height_field_weight,
	/*cushion optimization regularization*/ double _weight_width, double _weight_length, double _weight_similarity,double _weight_size) const
{
	// we refer trajectory curve also as rail curve; profile curve also as cross section curve; human face also as human head

	//=======================================================================================//
	//=================================== preprocessing =====================================//
	//=======================================================================================//

	// initialize cushion as a generic cushion, using saved data
	CushionSurface cushion;
	cushion.InitializeAsGenericCushion(data_path_);
	//Mesh generic_cushion = cushion.getSurfaceMesh();

	// tree data structures for face mesh
	AabbTree head_tree;
	head_tree.constructTree(human_face_);
	cushion.setHeadTree(&head_tree);
	MeshInfo head_mesh_info(human_face_);
	std::vector<std::array<double, 3>>* head_vertices = new std::vector<std::array<double, 3>>;
	std::vector<std::array<int, 3>>* head_triangles = new std::vector<std::array<int, 3>>;
	head_mesh_info.getVerticesAndTriangles(*head_vertices, *head_triangles);
	tmd::TriangleMeshDistance* head_tmd = new tmd::TriangleMeshDistance(*head_vertices, *head_triangles);// head distance tree
	// here we use only one human face expression
	int expressions = 1;
	std::vector<const tmd::TriangleMeshDistance*> head_tmds(expressions);
	head_tmds[0] = head_tmd;
	delete head_vertices;
	delete head_triangles;

	//=======================================================================================//
	//====================== cushion optimization initialization ============================//
	//=======================================================================================//

	cout << "======================== start cushion optimization initialization =======================" << endl;

	// 1. 
	// rail curve initialization

	int rail_degree = 3;
	int rail_smoothness = 2;
	// selected polyline vertex index for fitting a spline curve (20 polynomials)
	std::vector<int> segment_parts_idx = { 0, 96, 224, 374, 489, 628, 772, 897, 1028, 1101, 1176, 1245, 1335, 1464, 1597, 1728, 1870, 1990, 2132, 2271 };
	cushion.changeRailByFittingProjectedCushionPolylineOnFace(segment_parts_idx, rail_degree, rail_smoothness);

	// rail curve adjustment	
	Mesh connector;// connector distance tree
	OpenMesh::IO::read_mesh(connector, data_path_ + "/connector/connector.obj");
	MeshInfo connector_mesh_info(&connector);
	std::vector<std::array<double, 3>>* connector_vertices = new std::vector<std::array<double, 3>>;
	std::vector<std::array<int, 3>>* connector_triangles = new std::vector<std::array<int, 3>>;
	connector_mesh_info.getVerticesAndTriangles(*connector_vertices, *connector_triangles);
	tmd::TriangleMeshDistance* connector_tmd = new tmd::TriangleMeshDistance(*connector_vertices, *connector_triangles);

	// provide optimization common information independent of optimization algorithm
	RailOptimizeInfo* rail_info = new RailOptimizeInfo;
	rail_info->SetHeadSdf(head_tmds);
	rail_info->SetInitialRail(&(cushion.getRail()));
	rail_info->SetConnectorDistanceTree(connector_tmd);
	rail_info->FindBmcRailToConnectorDistance(cushion.getOriginalRail());
	rail_info->SetWeights(_objective, _width, _curvature, _symmetry, _align, _angle);

	NelderMeadRailOptimize* optimize_rail = new NelderMeadRailOptimize;
	optimize_rail->SetRailInfoAssistant(rail_info);
	Rail* optimal_rail = new Rail;
	optimal_rail->setRail(cushion.getRail().getRail());
	double optimal_f;
	int count_f_evaluation = 0;
	int num_restart = 0;
	double min_f_variance = 0.1;
	double step = 4;
	int convergence_check = 1000;
	int max_f_evaluation = 5000;
	optimize_rail->Optimize(*optimal_rail, &optimal_f, min_f_variance, step, convergence_check, max_f_evaluation, &count_f_evaluation, &num_restart);

	cushion.setRail(optimal_rail->getRail());

	delete optimize_rail;
	delete rail_info;
	delete optimal_rail;
	delete connector_tmd;
	delete connector_vertices;
	delete connector_triangles;

	// 2. 
	// cross section initialization

	double frame_relaxation_radius = 20.0;
	int sample_u = 50;
	int sample_v = 10;
	cushion.computeSampleFrames(frame_relaxation_radius, sample_u);// the u parameter refers to the rail curve parameter; v parameter refers to the 
	// 15 key parameters, uniformly along rail curve
	std::vector<double> sample_u_F2 = { 1.0 / 300.0, 18.0 / 300.0, 36.0 / 300.0, 60.0 / 300.0, 84.0 / 300.0, 108.0 / 300.0, 131.0 / 300.0, 146.0 / 300.0,
		160.0 / 300.0, 175.0 / 300.0, 197.0 / 300.0, 226.0 / 300.0, 253.0 / 300.0,  276.0 / 300.0, 294.0 / 300.0 };
	int cross_section_half1_degree = 2;
	int cross_section_half2_degree = 2;
	cushion.setInitialCustomKeyCrossSections_ChangableDegree(sample_u_F2, cross_section_half1_degree, cross_section_half2_degree);
	cushion.ConstructCrossSectionSurface(3);
	cushion.setFixedBoundary(data_path_ + "/connector/outside_boundary.txt", 101);



	//==========================================================================================//
	//=============================== cushion optimization =====================================//
	//==========================================================================================//

	cout << "=============================== start cushion optimization ===============================" << endl;

	// metric based on FEM simulation

	// 1.
	// construct signed distance field for face mesh

	// create interpenetration to mimick a human wearing the mask
	CurveDistanceToFace find_distance;// for calculating human face to connector distance
	find_distance.setHead(head_tmds[0]);
	std::vector<Mesh::Point> outside_boundary;
	ReadRhino read_boundary;
	read_boundary.Read3DPoints_Brackets(data_path_ + "/connector/outside_boundary.txt", outside_boundary, true);
	double connector_to_head_distance = find_distance.MinDistance(outside_boundary);
	double push_cushion_into_face = connector_to_head_distance / 4.0;
	Mesh human_face_pushed = *human_face_;
	for (auto& iVertex : human_face_pushed.vertices())
	{
		Mesh::Point translated_point = human_face_pushed.point(iVertex);
		translated_point[2] += push_cushion_into_face;
		human_face_pushed.set_point(iVertex, translated_point);
	}

	// construct translated head sdf
	MeshInfo pushed_head_mesh_info(&human_face_pushed);
	std::vector<std::array<double, 3>>* pushed_head_vertices = new std::vector<std::array<double, 3>>;
	std::vector<std::array<int, 3>>* pushed_head_triangles = new std::vector<std::array<int, 3>>;
	pushed_head_mesh_info.getVerticesAndTriangles(*pushed_head_vertices, *pushed_head_triangles);
	tmd::TriangleMeshDistance* pushed_head_tmd = new tmd::TriangleMeshDistance(*pushed_head_vertices, *pushed_head_triangles);

	// SDF in domain x[-100,100] y[-70,98] z[-38,43] 
	// this domian selection is based on the connector position(which is fixed and human faces should aligned to it)
	Eigen::AlignedBox3d domain;
	domain.max() = Eigen::Vector3d(100, 98, 43);
	domain.min() = Eigen::Vector3d(-100, -70, -38);
	std::array<unsigned int, 3> resolution = { {50, 50, 50} };
	Discregrid::CubicLagrangeDiscreteGrid* discrete_grid = new Discregrid::CubicLagrangeDiscreteGrid(domain, resolution);
	auto func = Discregrid::DiscreteGrid::ContinuousFunction{};
	func = [pushed_head_tmd](Eigen::Vector3d const& xi) {return pushed_head_tmd->signed_distance(xi).distance; };
	std::cout << "Generate human face SDF...";
	int sdf_index = discrete_grid->addFunction(func, true);

	std::vector<int> sdf_indexs(expressions);
	std::vector<const Discregrid::CubicLagrangeDiscreteGrid*> discrete_grids(expressions);
	std::vector<const tmd::TriangleMeshDistance*> pushed_head_tmds(expressions);
	sdf_indexs[0] = sdf_index;
	discrete_grids[0] = discrete_grid;
	pushed_head_tmds[0] = pushed_head_tmd;
	delete pushed_head_vertices;
	delete pushed_head_triangles;

	// 2.
	//  evaluation parameters

	// objective function parameters
	double weight_average_pressure = _weight_average_pressure;
	double force_distribution = _force_distribution;
	double area_distribution = _area_distribution;
	double convexity = _convexity;
	double weight_tetrahedral_inversion = _weight_tetrahedral_inversion;
	double weight_curvature = 0.0;// set as 0 since we do not use curvature to test self-intersection now
	double height_field_weight = _height_field_weight;

	double weight_width = _weight_width;
	double global_length_variance_weight = _weight_length;
	double shape_regularization = _weight_similarity;
	double size_regularization = _weight_size;

	double scale = 100;
	double force_epsilon = 0.0;// a threshold, not used anymore in the current implementation

	// simulation parameters
	double E = _E_module;
	double poisson_ratio = _poisson_ratio;
	double collision_wieght = 1;
	int max_iterations = 0;
	double epsilon_gradient = 1e-20;
	double epsilon_function = 0;
	double epsilon_step = 0;

	// tetrahedral cushion
	int num_layer = _num_layer;
	double thickness = _thickness;

	// 3.
	// optimization

	// prepare cushion_ optimization
	CrossSectionOptimizeInfo* cross_section_info = nullptr;
	cross_section_info = new CrossSectionOptimizeInfo(
		E, poisson_ratio,
		sdf_indexs, discrete_grids, head_tmds, collision_wieght, max_iterations, epsilon_gradient, epsilon_function, epsilon_step,
		weight_average_pressure, force_distribution, area_distribution,
		weight_tetrahedral_inversion,
		scale, force_epsilon,
		num_layer, thickness,
		&cushion, sample_v, cross_section_half1_degree, cross_section_half2_degree);

	// gradient based optimization
	cross_section_info->SetCurvatureBarrierWeight(weight_curvature);
	cross_section_info->SetHeightFieldWeight(height_field_weight);
	cross_section_info->SetWidthRegularizationWeight(weight_width);
	cross_section_info->SetGlobalLengthVarianceWeight(global_length_variance_weight);
	cross_section_info->SetShapeRegularizationWeight(shape_regularization);
	cross_section_info->SetSizeRegularizationWeight(size_regularization);
	cross_section_info->SetConvexityWeight(convexity);
	CushionOptimize* optimize_it = new CushionOptimize(cross_section_info);

	cout << "Cushion Optimization Begins ..." << endl;
	LARGE_INTEGER t1, t2, tc;
	QueryPerformanceFrequency(&tc);
	QueryPerformanceCounter(&t1);
	optimize_it->BfgsOptimize(cushion, _max_cushion_optimization_iterations, 0.01);
	QueryPerformanceCounter(&t2);
	cout << "CrossSection Design time = " << (double)(t2.QuadPart - t1.QuadPart) / (double)tc.QuadPart / 60.0 << "min" << endl;
	delete optimize_it;

	if (cross_section_info != nullptr)
	{
		delete cross_section_info;
		cross_section_info = nullptr;
	}
	for (int iExpression = 0; iExpression < expressions; ++iExpression)
	{
		delete discrete_grids[iExpression];
		delete pushed_head_tmds[iExpression];
		delete head_tmds[iExpression];
	}

	//==========================================================================================//
	//=====================================  output ============================================//
	//==========================================================================================//

	cout << "============================= assemble cushion with connector ============================" << endl;

	// get a denser sample of u and v
	cushion.computeSampleFrames(frame_relaxation_radius, 300);
	cushion.computeDiscreteCrossSections_NewFrames(50);
	cushion.constructSurfaceMesh();
	_cushion_surface = cushion.getSurfaceMesh();


	// offset cushion surface, then CSG union with connector to get mask interface

	double offset_length = 3.0;

	// 1. 
	// offset cushion to get a watertight cushion
	Mesh watertight_cushion;
	Mesh* watertight_cushion0 = new Mesh();
	Mesh* watertight_cushion_libigl = new Mesh();
	MeshOffset get_watertight_cushion;
	get_watertight_cushion.ThickenSurface_Libigl(_cushion_surface, *watertight_cushion_libigl, offset_length);
	get_watertight_cushion.ThickenSurface(_cushion_surface, *watertight_cushion0, offset_length + 1);
	get_watertight_cushion.FlipNormal(*watertight_cushion0);
	ConnectFrontNBack union_it;
	union_it.CsgIntersect(*watertight_cushion0, *watertight_cushion_libigl, watertight_cushion);
	delete watertight_cushion0;
	delete watertight_cushion_libigl;


	// 2. 
	// create a bridge between cushion and connector

	// get the connecting boundary on the connector
	std::vector<Mesh::Point> inside_points;
	read_boundary.Read3DPoints_Brackets(data_path_ + "/connector/inside_boundary.txt", inside_points, true);

	Mesh intermediate;
	ConnectFrontNBack connect_it;
	connect_it.CreateIntermediatePart(inside_points, outside_boundary, _cushion_surface, offset_length, intermediate);

	Mesh connector_n_intermediate;
	connect_it.CsgConnect(connector, intermediate, connector_n_intermediate);

	// 4. 
	// final adapter
	connect_it.CsgConnect(connector_n_intermediate, watertight_cushion, _mask_interface);
	

	cout << "success" << endl;
}
