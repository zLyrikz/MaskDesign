#include "Tetrahedralize.h"
#include <map>
#include <CGAL/Mesh_criteria_3.h>
//#include <CGAL/boost/graph/helpers.h>
#include <CGAL/make_mesh_3.h>
//#include <CGAL/refine_mesh_3.h>
#include <Eigen/Core>
#include <iostream>
using std::cout;
using std::endl;


// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;


Tetrahedralize::Tetrahedralize()
{
}

Tetrahedralize::~Tetrahedralize()
{
}

void Tetrahedralize::DoIt(const MeshCgal& _surface_mesh, TetrahedralMesh& _tetrahedral_mesh, double _cell_radius_edge_ratio, double _cell_size) const
{
    // Create input polyhedron
    Polyhedron polyhedron;
    CGAL::copy_face_graph(_surface_mesh, polyhedron);
    if (!CGAL::is_triangle_mesh(polyhedron)) {
        std::cerr << "Input geometry is not triangulated." << std::endl;
        return;
    }
    // Create domain
    Mesh_domain domain(polyhedron);

    // Set tetrahedron size (keep cell_radius_edge_ratio), ignore facets
    Mesh_criteria criteria(CGAL::parameters::cell_radius_edge_ratio = _cell_radius_edge_ratio, CGAL::parameters::cell_size = _cell_size);
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, CGAL::parameters::no_perturb(), CGAL::parameters::no_exude());
    // Mesh refinement (and make the output manifold)
    CGAL::refine_mesh_3(c3t3, domain, criteria, CGAL::parameters::manifold());

    ////test output
    //std::ofstream medit_file("out.mesh");
    //c3t3.output_to_medit(medit_file);
    //medit_file.close();

    CgalC3t3TetrahedralMesh(c3t3, _tetrahedral_mesh);

}

void Tetrahedralize::CgalTriangulation2TetrahedralMesh(const Tr& _cgal, TetrahedralMesh& _tetrahedral_mesh) const
{
    Eigen::MatrixX3d vertex(_cgal.number_of_vertices(), 3);// vertex coordinates
    Eigen::MatrixX4i topology_tetrahedra(_cgal.number_of_finite_cells(),4);// vertex index of each tetrahedron
    Eigen::MatrixX3i topology_boundary(_cgal.number_of_cells() - _cgal.number_of_finite_cells(),3);

    // will create a vertex index map, since cgal triangulation doesn't offer global vertex index
    std::map<const Tr::Vertex_handle, int> vertex_index;
    
    // get vertex
    int count = 0;
    for (auto& iV : _cgal.finite_vertex_handles())
    {
        vertex.row(count) = Eigen::Vector3d((*iV).point().x(), (*iV).point().y(), (*iV).point().z());
        vertex_index.insert(std::make_pair(iV, count));
        ++count;
    }
    assert(vertex_index.size() == _cgal.number_of_vertices());

    // get cell
    int count_cell = 0;
    for (auto& iC : _cgal.finite_cell_handles())
    {
        Eigen::Vector4i tetrahedron;
        
        // get the vertex index
        for (int iV = 0; iV < 4; ++iV)
        {
            const Tr::Vertex_handle& v = iC->vertex(iV);
            std::map<const Tr::Vertex_handle, int>::iterator map_iter = vertex_index.find(v);
            if (map_iter != vertex_index.end())
            {
                tetrahedron(iV) = map_iter->second;
            }
            else
            {
                cout << "ERROR 1 FROM Tetrahedralize::CgalTriangulation2TetrahedralMesh vertex not found" << endl;
            }
        }
        topology_tetrahedra.row(count_cell) = tetrahedron;
        ++count_cell;
    }

    // get boundary faces
    int count_face = 0;
    for (auto& iF : _cgal.finite_facets())
    {
        if (_cgal.is_infinite(iF))
        {
            // this finite facet is in a infinite cell, so it's a boundary facet
            Eigen::Vector3i triangle;
            // get the vertex index
            int count_vertex_added = 0;
            for (int iV = 0; iV < 4; ++iV)
            {
                if (iV != iF.second)
                {
                    // this vertex is not the infinite vertex in the cell
                    const Tr::Vertex_handle& v = iF.first->vertex(iV);
                    std::map<const Tr::Vertex_handle, int>::iterator map_iter = vertex_index.find(v);
                    if (map_iter != vertex_index.end())
                    {
                        //TODO: how to maintain consistent orientation
                        triangle(count_vertex_added) = map_iter->second;
                        ++count_vertex_added;
                    }
                    else
                    {
                        cout << "ERROR 2 FROM Tetrahedralize::CgalTriangulation2TetrahedralMesh vertex not found" << endl;
                    }
                }

            }
            assert(count_face < _cgal.number_of_cells() - _cgal.number_of_finite_cells());

            topology_boundary.row(count_face) = triangle;
            ++count_face;
        }

    }

    _tetrahedral_mesh.SetMesh(vertex, topology_tetrahedra, topology_boundary);
}

void Tetrahedralize::CgalC3t3TetrahedralMesh(const C3t3& _cgal, TetrahedralMesh& _tetrahedral_mesh) const
{
    // output then read...
    std::ofstream medit_file("Tetrahedralize.mesh");
    _cgal.output_to_medit(medit_file);
    medit_file.close();
    _tetrahedral_mesh.SetMesh("Tetrahedralize.mesh");

}
