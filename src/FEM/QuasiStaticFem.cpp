#include "QuasiStaticFem.h"
#include "WrapperAlglibFunction.h"
#include "../Utility/TypeConvert.h"
#include <windows.h>
#include <finitediff.hpp>
#include <iostream>
#include <fstream>
using std::cout;
using std::endl;
QuasiStaticFem::QuasiStaticFem()
{
    collision_penalty_function_ = std::bind(&QuasiStaticFem::CollisionPenalty, this, std::placeholders::_1);
    min_sdf_0_function_ = std::bind(&QuasiStaticFem::MinSdf0, this, std::placeholders::_1);
    sdf_function_ = std::bind(&QuasiStaticFem::Sdf, this, std::placeholders::_1);
    sdf_gradient_ = std::bind(&QuasiStaticFem::SdfGradient, this, std::placeholders::_1);
}

QuasiStaticFem::~QuasiStaticFem()
{
}

void QuasiStaticFem::SetInitialConfiguration(const Eigen::VectorXd& _initial_nodes)
{
    initial_nodes_ = _initial_nodes;
    dof_ = initial_nodes_.rows();
    num_node_ = dof_ / 3;
    if (dof_ % 3 != 0)
    {
        cout << "[ERROR QuasiStaticFem::SetInitialConfiguration] initial_nodes_ size % 3 != 0" << endl;
    }
}

void QuasiStaticFem::SetInitialDisplacement(const Eigen::VectorXd& _initial_displacement)
{
    initial_displacement_ = _initial_displacement;
}

void QuasiStaticFem::SetCollisionSdf(const tmd::TriangleMeshDistance* _other_objects, int _sdf_function_index, const Discregrid::CubicLagrangeDiscreteGrid* _real_sdf)
{
	sdf_other_objects_ = _other_objects;
    sdf_function_index_ = _sdf_function_index;
    real_sdf_ = _real_sdf;
}

void QuasiStaticFem::SetCollisionWeight(double _collision_weight)
{
	collision_weight_ = _collision_weight;
}

void QuasiStaticFem::CgOptimize(alglib::ae_int_t _max_iteration, double _epsilon_gradient, double _epsilon_fucntion, double _epsilon_step)
{
    // First, we create optimizer object and tune its properties
    //
    TypeConvert type_converter;
    alglib::real_1d_array initial_displacement;
    initial_displacement.setlength(initial_displacement_.rows());
    type_converter.EigenVector2AlglibArray(initial_displacement_, initial_displacement);


    alglib::mincgstate state;// the state of the algorithm
    alglib::mincgcreate(initial_displacement, state);

    alglib::mincgsetcond(state, _epsilon_gradient, _epsilon_fucntion, _epsilon_step, _max_iteration);

    // variable scale for stopping criteria
    alglib::real_1d_array scale;
    scale.setlength(dof_);
    for (int iDof = 0; iDof < dof_; ++iDof)
    {
        scale[iDof] = 5;
    }
    alglib::mincgsetscale(state, scale);

    //
    // Optimize and evaluate results
    // we have to change our class member function to a normal function
    // check https://www.zhihu.com/question/411373289
    WrapperAlglibFunction::broker_ObjcetiveFunction_alglib = 
        std::bind(&QuasiStaticFem::ObjcetiveFunction_alglib, this, 
            std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);

    alglib::mincgoptimize(state, WrapperAlglibFunction::wrapper_ObjcetiveFunction_alglib);

    alglib::mincgreport report;
    alglib::real_1d_array final_displacement;
    mincgresults(state, final_displacement, report);
    type_converter.AlglibArray2EigenVector(final_displacement, final_displacement_);

    // result report
    // TerminationType field contains completion code, which can be :
    //     -8    internal integrity control detected  infinite or NAN  values  in
    //    function / gradient.Abnormal termination signalled.
    //    1    relative function improvement is no more than EpsF.
    //    2    relative step is no more than EpsX.
    //    4    gradient norm is no more than EpsG
    //    5    MaxIts steps was taken
    //    7    stopping conditions are too stringent,
    //    further improvement is impossible,
    //    X contains best point found so far.
    //    8    terminated by user who called mincgrequesttermination().X contains
    //    point which was "current accepted" when  termination  request  was
    //    submitted.
    cout << "total number of inner iterations " << report.iterationscount << endl;
    cout << "number of gradient evaluations " << report.nfev << endl;
    cout << "termination type " << report.terminationtype << endl;
    // get the final objective function value
    double final_objective_function_value = 0.0;
    alglib::real_1d_array final_objective_gradient_value;
    final_objective_gradient_value.setlength(dof_);
    ObjcetiveFunction_alglib(final_displacement, final_objective_function_value, final_objective_gradient_value, nullptr);
}

void QuasiStaticFem::LbfgsOptimize(alglib::ae_int_t _max_iteration, double _epsilon_gradient, double _epsilon_fucntion, double _epsilon_step, bool _print_result)
{
    alglib::minlbfgsstate state;// the state of the algorithm
    int M = 5;
    double scale = 5;
    LbfgsOptimizeSetUp(state, M, scale, _max_iteration, _epsilon_gradient, _epsilon_fucntion, _epsilon_step);

    //
    // Optimize and evaluate results
    // 
    // we have to change our class member function to a normal function
    // check https://www.zhihu.com/question/411373289
    WrapperAlglibFunction::broker_ObjcetiveFunction_alglib =
        std::bind(&QuasiStaticFem::ObjcetiveFunction_alglib, this,
            std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);

    // initial objective function value
    // get the final objective function value
    TypeConvert type_converter0;
    alglib::real_1d_array initial_displacement;
    initial_displacement.setlength(initial_displacement_.rows());
    type_converter0.EigenVector2AlglibArray(initial_displacement_, initial_displacement);
    double initial_objective_function_value = 0.0;
    alglib::real_1d_array initial_objective_gradient_value;
    initial_objective_gradient_value.setlength(dof_);
    ObjcetiveFunction_alglib(initial_displacement, initial_objective_function_value, initial_objective_gradient_value, nullptr);
    //cout << "initial_objective_function_value=" << initial_objective_function_value << endl;
    //cout << "initial displacement norm=" << initial_displacement_.norm() << endl;

    alglib::minlbfgsoptimize(state, WrapperAlglibFunction::wrapper_ObjcetiveFunction_alglib);

    alglib::minlbfgsreport report;
    alglib::real_1d_array final_displacement;
    minlbfgsresults(state, final_displacement, report);
    TypeConvert type_converter;
    type_converter.AlglibArray2EigenVector(final_displacement, final_displacement_);

    if (_print_result)
    {
        // result report
        // TerminationType field contains completion code, which can be :
        //     -8    internal integrity control detected  infinite or NAN  values  in
        //    function / gradient.Abnormal termination signalled.
        //    1    relative function improvement is no more than EpsF.
        //    2    relative step is no more than EpsX.
        //    4    gradient norm is no more than EpsG
        //    5    MaxIts steps was taken
        //    7    stopping conditions are too stringent,
        //    further improvement is impossible,
        //    X contains best point found so far.
        //    8    terminated by user who called mincgrequesttermination().X contains
        //    point which was "current accepted" when  termination  request  was
        //    submitted.
        //
        cout << "total number of inner iterations " << report.iterationscount << endl;
        cout << "number of gradient evaluations " << report.nfev << endl;
        cout << "termination type " << report.terminationtype << endl;
    }

    // get the final objective function value
    double final_objective_function_value = 0.0;
    alglib::real_1d_array final_objective_gradient_value;
    final_objective_gradient_value.setlength(dof_);
    ObjcetiveFunction_alglib(final_displacement, final_objective_function_value, final_objective_gradient_value, nullptr);
}

void QuasiStaticFem::LbfgsPenaltyOptimize(alglib::ae_int_t _max_outer_iteration, double _max_collision, double _penalty_increase, 
    alglib::ae_int_t _max_inner_iteration, double _epsilon_gradient, double _epsilon_fucntion, double _epsilon_step, int _M,
    std::vector<Eigen::VectorXd>& _final_displacement_each_step)
{
    alglib::minlbfgsstate state;// the state of the algorithm
    double scale = 5;
    LbfgsOptimizeSetUp(state, _M, scale, _max_inner_iteration, _epsilon_gradient, _epsilon_fucntion, _epsilon_step);

    //
    // Optimize and evaluate results
    // 
    // we have to change our class member function to a normal function
    // check https://www.zhihu.com/question/411373289
    WrapperAlglibFunction::broker_ObjcetiveFunction_alglib =
        std::bind(&QuasiStaticFem::ObjcetiveFunction_alglib, this,
            std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);

    // first iteration
    alglib::real_1d_array final_displacement;
    alglib::minlbfgsoptimize(state, WrapperAlglibFunction::wrapper_ObjcetiveFunction_alglib);

    alglib::minlbfgsreport report;
    alglib::minlbfgsresults(state, final_displacement, report);

    int total_inner_iteration = report.iterationscount;
    int total_gradient_iteration = report.nfev;
    std::vector<int> termination_type{int(report.terminationtype)};
    int outer_iteration_num = 1;

    // more iteration
    for (int iOuter = 1; iOuter < _max_outer_iteration; ++iOuter)
    {
        // If no collision already, then we end outer iteration
        // 
        _final_displacement_each_step.emplace_back(Eigen::VectorXd());
        TypeConvert type_converter;
        type_converter.AlglibArray2EigenVector(final_displacement, _final_displacement_each_step.back());

        if (CollisionPenalty(_final_displacement_each_step.back()) < _max_collision)
        {
            break;
        }


        // increase the collision penalty weight
        collision_weight_ *= _penalty_increase;

        // algorithm is restarted
        alglib::minlbfgsrestartfrom(state, final_displacement);

        alglib::minlbfgsoptimize(state, WrapperAlglibFunction::wrapper_ObjcetiveFunction_alglib);
        alglib::minlbfgsresults(state, final_displacement, report);

        total_inner_iteration += report.iterationscount;
        total_gradient_iteration += report.nfev;
        termination_type.push_back(report.terminationtype);
        outer_iteration_num += 1;


    }

    TypeConvert type_converter;
    type_converter.AlglibArray2EigenVector(final_displacement, final_displacement_);
    _final_displacement_each_step.emplace_back(final_displacement_);

    // result report
    // TerminationType field contains completion code, which can be :
    //     -8    internal integrity control detected  infinite or NAN  values  in
    //    function / gradient.Abnormal termination signalled.
    //    1    relative function improvement is no more than EpsF.
    //    2    relative step is no more than EpsX.
    //    4    gradient norm is no more than EpsG
    //    5    MaxIts steps was taken
    //    7    stopping conditions are too stringent,
    //    further improvement is impossible,
    //    X contains best point found so far.
    //    8    terminated by user who called mincgrequesttermination().X contains
    //    point which was "current accepted" when  termination  request  was
    //    submitted.
    //
    cout << "total number of outer iterations " << outer_iteration_num << endl;
    cout << "total number of inner iterations " << total_inner_iteration << endl;
    cout << "number of gradient evaluations " << total_gradient_iteration << endl;
    cout << "termination type ";
    for (int iInner = 0; iInner < termination_type.size(); ++iInner)
    {
        cout << termination_type[iInner] << " ";

    }
    cout << endl;

    // get the final objective function value
    double final_objective_function_value = 0.0;
    alglib::real_1d_array final_objective_gradient_value;
    final_objective_gradient_value.setlength(dof_);
    ObjcetiveFunction_alglib(final_displacement, final_objective_function_value, final_objective_gradient_value, nullptr);
}

std::function<double(const Eigen::VectorXd&, Eigen::VectorXd&)>* QuasiStaticFem::GetElasticEnergyValueAndGradient()
{
    return &elastic_energy_value_and_gradient_;
}

std::function<void(const Eigen::VectorXd&)>* QuasiStaticFem::GetInfoFunction()
{
    return &info_function_;
}

void QuasiStaticFem::ObjcetiveFunction_alglib(const alglib::real_1d_array& _displacement, double& _function_value, alglib::real_1d_array& _gradient_value, void* _ptr) const
{
    LARGE_INTEGER t11, t2, tc1, t3, t4, t5, t6;
    QueryPerformanceFrequency(&tc1);
    QueryPerformanceCounter(&t11);


    TypeConvert type_converter;
    Eigen::VectorXd displacement;
    type_converter.AlglibArray2EigenVector(_displacement, displacement);


    _function_value = 0.0;
    Eigen::VectorXd gradient;
    gradient.setZero(displacement.rows());

    // collision
    double collison_penalty = 0.0;
    Eigen::VectorXd collision_gradient;
    CollisionValueAndGradient(displacement, collison_penalty, collision_gradient);
    _function_value += collision_weight_ * collison_penalty;
    gradient += collision_weight_ * collision_gradient;

    QueryPerformanceCounter(&t2);
    //cout <<"collison time=" << (t2.QuadPart - t11.QuadPart) / (double)tc1.QuadPart * 1000 << "ms" << endl;

    // energy
    Eigen::VectorXd energy_gradient;
    double elastic_energy = elastic_energy_value_and_gradient_(displacement, energy_gradient);
    _function_value += elastic_energy;
    gradient += energy_gradient;

    QueryPerformanceCounter(&t3);

    // save simulation process
    if (0)
    {
        static int total_iteration = 0;
        static int iterations = 0;
        int sample = 250;// include two end points
        double min_displce = 0;
        double max_displce = 583.42;// 584.875
        std::vector<double> displce_middle(sample);
        for (int iSample = 0; iSample < sample; ++iSample)
        {
            displce_middle[iSample] = iSample * (max_displce - min_displce) / double(sample - 1) + min_displce;
        }
        displce_middle.push_back(displce_middle.back() + 1);
        static std::vector<double> f_middle{ 6.95985e+06 };// this value larger than initial f
        std::ofstream output_displace("C:\\Users\\15539\\Desktop\\MaskDesign\\Video\\data\\simulation\\displace_free.txt", std::ios::app);
        if (iterations < sample && _function_value < f_middle[iterations])
        {
            double displace_norm = displacement.norm();
            if (displace_norm >= displce_middle[iterations] /*&& displace_norm < displce_middle[iterations + 1]*/)
            {
                cout << "sample " << iterations << ", total iterations=" << total_iteration << ", f=" << _function_value << endl;
                ++iterations;
                f_middle.push_back(_function_value);
                output_displace << displacement.transpose() << endl;
            }
        }
        ++total_iteration;
    }

    //cout << "gradient=\n" << gradient << endl;

    // gradient output
    type_converter.EigenVector2AlglibArray(gradient, _gradient_value);
}

void QuasiStaticFem::ConstraintOptimizationFunctions_alglib(const alglib::real_1d_array& _displacement, alglib::real_1d_array& _functions, alglib::real_2d_array& _jacobian, void* ptr)
{
    // TODO update sdf_other_objects_ to real_sdf
    TypeConvert type_converter;
    Eigen::VectorXd displacement;
    type_converter.AlglibArray2EigenVector(_displacement, displacement);

    // i-th collision constraints = (-1) * sdf(xi) <= 0 (note the -1)
    Eigen::VectorXd nodes;
    GetNodeFromDisplacement(displacement, nodes);
    fd::AccuracyOrder accuracy = fd::SECOND;// gradient finite difference approximation order
    double h = 1e-8;                        // gradient finite difference approximation step
    for (int iNode = 0; iNode < num_node_; ++iNode)
    {
        tmd::Result result = sdf_other_objects_->signed_distance(
            { nodes[space_dim_ * iNode],
            nodes[space_dim_ * iNode + 1],
            nodes[space_dim_ * iNode + 2] });
        
        _functions[iNode + 1] = -result.distance;


        // gradient

        Eigen::Vector3d node(nodes[space_dim_ * iNode], nodes[space_dim_ * iNode + 1], nodes[space_dim_ * iNode + 2]);
        Eigen::VectorXd gradient_node;// size = 3
        fd::finite_gradient(node, sdf_function_, gradient_node, accuracy, h);       

        // put the gradient at current to the correct position in the whole gradient
        for (int iX = 0; iX < dof_; ++iX)
        {
            _jacobian[iNode + 1][iX] = 0;
        }
        for (int iCoordinate = 0; iCoordinate < 3; ++iCoordinate)
        {
            _jacobian[iNode + 1][space_dim_ * iNode + iCoordinate] = -gradient_node[iCoordinate];
        }
    }

    Eigen::VectorXd energy_gradient;
    double elastic_energy = elastic_energy_value_and_gradient_(displacement, energy_gradient);
    _functions[0] = elastic_energy;
    for (int iX = 0; iX < dof_; ++iX)
    {
        _jacobian[0][iX] = energy_gradient(iX);
    }

    cout << "_function_value = energy =" << elastic_energy << endl;
}

void QuasiStaticFem::GetNodeFromDisplacement(const Eigen::VectorXd& _displacement, Eigen::VectorXd& _nodes) const
{
    _nodes = initial_nodes_;
    _nodes += _displacement;
}

void QuasiStaticFem::GetFinalDisplacement(Eigen::VectorXd& _final_displacement) const
{
    _final_displacement = final_displacement_;
}

void QuasiStaticFem::FinalCollisionHessian(Eigen::SparseMatrix<double>& _hessian) const
{
    CollisionHessian(final_displacement_, _hessian);
    _hessian *= collision_weight_;
}

void QuasiStaticFem::LbfgsOptimizeSetUp(alglib::minlbfgsstate& _state, int _M, double _scale, alglib::ae_int_t _max_iteration, double _epsilon_gradient, double _epsilon_fucntion, double _epsilon_step)
{
    // the procedure follows alglib doc on L-BFGS https://www.alglib.net/optimization/lbfgsandcg.php

    // First, we create optimizer object and tune its properties
    //
    TypeConvert type_converter;
    alglib::real_1d_array initial_displacement;
    initial_displacement.setlength(initial_displacement_.rows());
    type_converter.EigenVector2AlglibArray(initial_displacement_, initial_displacement);

    // generates the approximation of an inverse Hessian matrix by
    // using information about the last M steps of the algorithm
    alglib::minlbfgscreate(dof_, _M, initial_displacement, _state);

    alglib::minlbfgssetcond(_state, _epsilon_gradient, _epsilon_fucntion, _epsilon_step, _max_iteration);

    // variable scale for stopping criteria
    alglib::real_1d_array scale;
    scale.setlength(dof_);
    for (int iDof = 0; iDof < dof_; ++iDof)
    {
        scale[iDof] = _scale;
    }
    alglib::minlbfgssetscale(_state, scale);

}

double QuasiStaticFem::CollisionPenalty(const Eigen::VectorXd& _displacement) const
{
    // TODO update sdf_other_objects_ to real_sdf

    Eigen::VectorXd nodes;
    GetNodeFromDisplacement(_displacement, nodes);


    double penalty = 0.0;
    for (int iNode = 0; iNode < num_node_; ++iNode)
    {
        tmd::Result result = sdf_other_objects_->signed_distance(
            { nodes[space_dim_ * iNode],
            nodes[space_dim_ * iNode + 1],
            nodes[space_dim_ * iNode + 2] });

        // ( min{sdf(x), 0} )^2
        if (result.distance < 0)
        {
            penalty += result.distance * result.distance;

        }

    }
    return penalty;
}

void QuasiStaticFem::CollisionValueAndGradient(const Eigen::VectorXd& _displacement, double& _value, Eigen::VectorXd& _gradient) const
{
    if (real_sdf_ == nullptr)
    {
        cout << "ERROR !!!!!!!!!!!!!! QuasiStaticFem::CollisionValueAndGradient real_sdf_=null" << endl;
        return;
    }

    Eigen::VectorXd nodes;
    GetNodeFromDisplacement(_displacement, nodes);

    _value = 0.0;
    _gradient.setZero(_displacement.rows());
    fd::AccuracyOrder accuracy = fd::SECOND;// gradient finite difference approximation order
    double h = 1e-4;                        // gradient finite difference approximation step
    for (int iNode = 0; iNode < num_node_; ++iNode)
    {
        double distance = real_sdf_->interpolate(sdf_function_index_, { nodes[space_dim_ * iNode], nodes[space_dim_ * iNode + 1], nodes[space_dim_ * iNode + 2] });

        // ( min{sdf(x), 0} )^2
        if (distance < 0.0)
        {
            // value
            _value += distance * distance;

            // gradient w.r.t node equal to displacement
            // need to use  min{sdf(x), 0} to compute the gradient when sdf near 0, 
            // and if sdf is far from 0, just use sdf(x)

            Eigen::Vector3d gradient_node;// size = 3
            real_sdf_->interpolate(sdf_function_index_, { nodes[space_dim_ * iNode], nodes[space_dim_ * iNode + 1], nodes[space_dim_ * iNode + 2] }, &gradient_node);

            gradient_node *= (2.0 * distance);

            // put the gradient at current to the correct position in the whole gradient
            for (int iCoordinate = 0; iCoordinate < 3; ++iCoordinate)
            {
                _gradient[space_dim_ * iNode + iCoordinate] = gradient_node[iCoordinate];
            }
        }

    }
}

double QuasiStaticFem::MinSdf0(const Eigen::Vector3d& _node) const
{
    const tmd::Result& result = sdf_other_objects_->signed_distance({ _node[0], _node[1], _node[2] });

    // min{sdf(x), 0}
    if (result.distance < 0)
    {
        return result.distance;
    }
    else
    {
        return 0.0;
    }
    
}

double QuasiStaticFem::Sdf(const Eigen::Vector3d& _node) const
{
    // TODO update sdf_other_objects_ to real_sdf

    const tmd::Result& result = sdf_other_objects_->signed_distance({ _node[0], _node[1], _node[2] });
    return result.distance;
}

Eigen::Vector3d QuasiStaticFem::SdfGradient(const Eigen::Vector3d& _node) const
{
    Eigen::Vector3d gradient_node;
    real_sdf_->interpolate(sdf_function_index_,
        { _node[0],  _node[1],  _node[2] }, &gradient_node);
    return gradient_node;
}

void QuasiStaticFem::CollisionHessian(const Eigen::VectorXd& _displacement, Eigen::SparseMatrix<double>& _hessian) const
{
    if (real_sdf_ == nullptr)
    {
        cout << "ERROR !!!!!!!!!!!!!! QuasiStaticFem::CollisionValueAndGradient real_sdf_=null" << endl;
        return;
    }
    Eigen::VectorXd nodes;
    GetNodeFromDisplacement(_displacement, nodes);

    std::vector<Eigen::Triplet<double>> triplet_list;
    triplet_list.reserve(9 * num_node_);
    for (int iNode = 0; iNode < num_node_; ++iNode)
    {
        //tmd::Result result = sdf_other_objects_->signed_distance(
        //    { nodes[space_dim_ * iNode],
        //    nodes[space_dim_ * iNode + 1],
        //    nodes[space_dim_ * iNode + 2] });
        double distance = real_sdf_->interpolate(sdf_function_index_, { nodes[space_dim_ * iNode], nodes[space_dim_ * iNode + 1], nodes[space_dim_ * iNode + 2] });


        // ( min{sdf(x), 0} )^2
        // denote f(x)=min{sdf(x), 0}
        // hessian ( min{sdf(x), 0} )^2 = 2 * (Df^T * Df + f*D^2f) = 2*(term1 + term2)
        // term1 = Df^T * Df
        // term2 = f*D^2f
        if (distance < 0.0)
        {
            Eigen::Vector3d node(nodes[space_dim_ * iNode], nodes[space_dim_ * iNode + 1], nodes[space_dim_ * iNode + 2]);

            // term1
            //Eigen::VectorXd gradient_node;// gradient of min{sdf(x), 0}, size = 3
            //fd::finite_gradient(node, min_sdf_0_function_, gradient_node, fd::SECOND, 1.0e-5);
            //Eigen::MatrixXd hessian_term1 = gradient_node * gradient_node.transpose();// 3x3 matrix
            Eigen::Vector3d gradient_node;
            real_sdf_->interpolate(sdf_function_index_, 
                { nodes[space_dim_ * iNode], nodes[space_dim_ * iNode + 1], nodes[space_dim_ * iNode + 2] }, &gradient_node);
            Eigen::MatrixXd hessian_term1 = gradient_node * gradient_node.transpose();// 3x3 matrix


            // term2
            Eigen::MatrixXd hessian_term2;// 3x3 matrix
            //fd::finite_hessian(node, min_sdf_0_function_, hessian_term2, fd::SECOND, 1.0e-5);//1e-5
            fd::finite_jacobian(node, sdf_gradient_, hessian_term2, fd::SECOND, 1.0e-5);// hessian should be symmetric(so transpose doesn't matter)
            hessian_term2 *= distance;


            // put the hessian at current to the correct position in the whole hessian
            for (int iRow = 0; iRow < 3; ++iRow)
            {
                for (int iCol = 0; iCol < 3; ++iCol)
                {
                    double value = 2.0 * (hessian_term1(iRow, iCol) + hessian_term2(iRow, iCol));
                    triplet_list.push_back(Eigen::Triplet<double>(space_dim_ * iNode + iRow, space_dim_ * iNode + iCol, value));
                }
            }
        }

    }
    _hessian.resize(dof_, dof_);
    _hessian.setFromTriplets(triplet_list.begin(), triplet_list.end());
}
