#pragma once
//#include "../MeshViewer/MeshDefinition.h"
#include "../MeshProcessing/TetrahedralMesh.h"
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K2;
typedef CGAL::Surface_mesh<K2::Point_3> MeshCgal;

// check CGAL doc in 3D Mesh Generation
// tutorial file: Mesh_3/mesh_polyhedral_domain.cpp
// Domain
typedef CGAL::Polyhedron_3<K2> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K2> Mesh_domain;
#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default, Concurrency_tag>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

class Tetrahedralize
{
public:
	Tetrahedralize();
	~Tetrahedralize();

	void DoIt(const MeshCgal& _surface_mesh, TetrahedralMesh& _tetrahedral_mesh, double _cell_radius_edge_ratio, double _cell_size) const;

	void CgalTriangulation2TetrahedralMesh(const Tr& _cgal, TetrahedralMesh& _tetrahedral_mesh) const;
	void CgalC3t3TetrahedralMesh(const C3t3& _cgal, TetrahedralMesh& _tetrahedral_mesh) const;
};

