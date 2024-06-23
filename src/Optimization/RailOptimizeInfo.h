#pragma once
#include "../MaskInterface/Rail.h"
#include "../MeshProcessing/TriangleMeshDistance.h"

// transfer between rail data and general raw data for an optimization algorithm
// provide rail objective function and feasible check


// TODO prvide deep copy functions
// we have dynamic allocate pointer member

class RailOptimizeInfo
{
public:
    RailOptimizeInfo();
	~RailOptimizeInfo();


    //void SetHead(const MeshCgal* _head);
    // multiple head models represent different expressions
    void SetHeadSdf(const std::vector<const tmd::TriangleMeshDistance*>& _head_sdfs);
    void SetInitialRail(const Rail* _initial_rail);
    void SetBoxConstraintRadius(double _box_edge_radius);

    void SetWeights(double _objective = 1e3, double _cushion_width = 2e-4, double _strain_energy = 10.0 * 1e-6,
        double _symmetry_energy = 0.01, double _size_regularity = 1e8, double _angle_energy = 80.0);

    void SetConnectorDistanceTree(const tmd::TriangleMeshDistance* _connector);
    void FindBmcRailToConnectorDistance(const Rail& _bmc_rail);
public:

    // input x size = spline space dim * point dim; (spline space dim = K)
    // control_points[iPoint][iDim] = _x[iPoint * point_dim + iDim]
    void GetRailFromX(Rail& _rail, const double* _x) const;
    void GetXFromRail(const Rail& _rail, double* _x) const;

    int GetDimX() const;
    const double* GetInitialX() const;
public:
    double ObjectiveFunction(const double* _x);
    bool IsFeasible(const double* _x) const;

private:
    // approximate average curvature
    // we assume rail is degree 3, smoothness 2
    double StrainEnergy(Rail& _rail) const;
    double StrainEnergyIntegrand(double _x, Rail& _rail) const;// (D^2 r)^2 

    // top point has max y value, bottom has min
    void FindRailTopBottomPoint(const std::vector<Mesh::Point>& _sample_rail_points, int& _top_id, int& _bottom_id) const;
    double SymmetryEnergy(const std::vector<Mesh::Point>& _sample_rail_points, int _top_id, int _bottom_id, int _sample_y) const;
    double SizeRegularity(double _max_y, double _min_y) const;
    double WideTopAngleEnergy(const double* _x) const;
    double WideTopEnergy(const double* _x) const;
    double WideUpperCurve(const double* _x) const;
    double DistanceToConnector(const std::vector<Mesh::Point>& _sample_rail_points) const;
    double CushionWidth(const std::vector<Mesh::Point>& _sample_rail_points) const;
    // a simple square box constriant
    bool IsFeasible_BoxConstraint(const double* _initial_x, const double* _current_x, const double& _box_edge_radius) const;

    void SetRailControlPointInfo(int _point_num, int _point_dim);


private:
    // special head data structure for evaluation
    //const MeshCgal* head_;
    std::vector<const tmd::TriangleMeshDistance*> head_sdfs_;
    Rail rail_;// store initial rail, and will be changed during objective function evaluation
    double* initial_x_;

    // dim x of optimization = num_control_point_ * control_point_dim_
    int num_control_point_;
    int control_point_dim_;
    int dim_x_;
    double initial_max_y_;//top and bottom y value of the rail curve
    double initial_min_y_;

    double box_edge_radius_;

    int sample_num_;
    const tmd::TriangleMeshDistance* connector_;
    double distance_to_connector_;// the original BMC rail distance to connector

    // weights
    double objective_;
    double cushion_width_; 
    double strain_energy_;
    double symmetry_energy_;
    double size_regularity_;
    double angle_energy_;
};

