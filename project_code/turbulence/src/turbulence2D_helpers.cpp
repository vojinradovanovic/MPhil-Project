#include <cassert>
#include <cmath>
#include <iostream>

#include "turbulence2D.hpp"

// protected

void Turbulence2D::calc_G_field_incompressible(std::vector<double>& G_field)
{
    // clear G field
    G_field.clear();
    G_field.resize(num_cells);

    // calculate grad_velocity_transpose_field
    std::vector<Tensor<Point>> grad_velocity_transpose_field;
    calc_grad_velocity_transpose_field(grad_velocity_transpose_field);

    // calculate symmetric velocity gradient
    std::vector<Tensor<Point>> sym_velocity_grad_field;
    calc_sym_velocity_grad_field(sym_velocity_grad_field, grad_velocity_transpose_field);

    // calculate G field
    for(int cell_id = 0; cell_id < num_cells; cell_id++)
    {
        G_field[cell_id] = 2 * nu_turb_field[cell_id] * sym_velocity_grad_field[cell_id]
                                                        .double_dot(grad_velocity_transpose_field[cell_id]);

                            // same as computing sym_velocity_grad_field : grad_velocity_field

        // k is excluded from computing G, as, in incompressible flow, 2/3 k I : grad(u) = 2/3 k div(u) = 0
    }
}

void Turbulence2D::calc_sym_velocity_grad_field(std::vector<Tensor<Point>>& sym_velocity_grad_field,
                                                const std::vector<Tensor<Point>>& grad_vel_t_field)
{
    // clear sym_velocity_grad_field
    sym_velocity_grad_field.clear();
    sym_velocity_grad_field.resize(num_cells);  // initially clear tensor field

    for(int cell_id = 0; cell_id < num_cells; cell_id++)
    {
        Tensor<Point> grad_vel_t = grad_vel_t_field[cell_id];
        Tensor<Point> grad_vel = grad_vel_t.transpose();

        // add grad_vel and grad_vel_t to sym_velocity_grad_field
        sym_velocity_grad_field[cell_id].add(grad_vel);
        sym_velocity_grad_field[cell_id].add(grad_vel_t);

        // multiply field by 0.5
        sym_velocity_grad_field[cell_id].scale_by(0.5);
    }
}

void Turbulence2D::calc_grad_velocity_transpose_field(std::vector<Tensor<Point>>& grad_vel_t_field)
{
    // initialisations
    std::vector<std::vector<double>> velocity_dim_fields(N_DIMS, std::vector<double>(num_cells));
    std::vector<std::vector<Point>> grad_velocity_dim_fields(N_DIMS, std::vector<Point>(num_cells));

    // populate velocity_dim_fields
    for(Dimensions dim: ALL_DIMS)
        extract_dim_field(velocity_field, velocity_dim_fields[dim], dim);

    // populate grad_velocity_dim_fields
    for(Dimensions dim: ALL_DIMS)
    {
        // isolate boundary conditions of velocity component
        std::vector<double> velocity_dim_field;
        extract_dim_field(velocity_field, velocity_dim_field, dim);

        std::vector<std::vector<double>> velocity_dim_bc_values;
        extract_dim_field(velocity_bc_values, velocity_dim_bc_values, dim);

        calc_grad_field_of_scalar(velocity_dim_fields[dim], grad_velocity_dim_fields[dim],
                                    velocity_bcs, velocity_dim_bc_values, velocity_intpol_method);
    }

    // populate output vector, grad_vel_t_field
    grad_vel_t_field.clear();
    grad_vel_t_field.resize(num_cells); // initially clear tensor field

    for(int cell_id = 0; cell_id < num_cells; cell_id++)
    {
        for(Dimensions dim: ALL_DIMS)
            grad_vel_t_field[cell_id].push_component(grad_velocity_dim_fields[dim][cell_id]);
    }
}



void Turbulence2D::calc_nu_turb_field()
{
    // resize field appropriately
    nu_turb_field.resize(num_cells);

    for(int cell_id = 0; cell_id < num_cells; cell_id++)
        nu_turb_field[cell_id] = C_nu * k_field[cell_id] * k_field[cell_id] / epsilon_field[cell_id];
}

void Turbulence2D::calc_scaled_nu_eff_fields(std::vector<double>& scaled_nu_eff_field,
                                                std::vector<BoundaryConditions>& scaled_nu_eff_bcs,
                                                std::vector<std::vector<double>>& scaled_nu_eff_bc_values, double scaling_coeff)
{
    // resize fields appropriately
    scaled_nu_eff_field.resize(num_cells);
    scaled_nu_eff_bcs.resize(num_patches);
    scaled_nu_eff_bc_values.resize(num_patches);

    for(int patch_id = 0; patch_id < num_patches; patch_id++)
    {
        std::size_t num_patch_faces = patches[patch_id].get_face_ids().size();
        scaled_nu_eff_bc_values[patch_id].resize(num_patch_faces);
    }

    // populate scaled_nu_eff_field
    for(int cell_id = 0; cell_id < num_cells; cell_id++)
    {
        scaled_nu_eff_field[cell_id] = nu_field[cell_id] + nu_turb_field[cell_id] / scaling_coeff;
    }

    if(nu_bcs == nu_turb_bcs)
    {
        // in case of equal boundary conditions, copy scaled_nu_eff_bcs from nu_bcs
        scaled_nu_eff_bcs = nu_bcs;

        // populate scaled_nu_eff_bc_values
        for(int patch_id = 0; patch_id < num_patches; patch_id++)
        {
            std::size_t num_patch_faces = patches[patch_id].get_face_ids().size();

            for(std::size_t pos_in_patch = 0; pos_in_patch < num_patch_faces; pos_in_patch++)
                scaled_nu_eff_bc_values[patch_id][pos_in_patch] = nu_bc_values[patch_id][pos_in_patch]
                                                                    + nu_turb_bc_values[patch_id][pos_in_patch] / scaling_coeff;
        }
    }
    else
    {
        // populate scaled_nu_eff_bc_values with dirichlet values
        for(int patch_id = 0; patch_id < num_patches; patch_id++)
        {
            scaled_nu_eff_bcs[patch_id] = BoundaryConditions::DIRICHLET;

            std::vector<int> patch_face_ids = patches[patch_id].get_face_ids();
            std::size_t num_patch_faces = patch_face_ids.size();
            
            for(std::size_t pos_in_patch = 0; pos_in_patch < num_patch_faces; pos_in_patch++)
            {
                int face_id = patch_face_ids[pos_in_patch];

                // interpolate viscosities
                double nu_interpolated = interpolate_to_face(nu_field, nu_bcs, nu_bc_values,
                                                                face_id, diffusivity_intpol_method);
                double nu_turb_interpolated = interpolate_to_face(nu_turb_field, nu_turb_bcs, nu_turb_bc_values,
                                                                    face_id, diffusivity_intpol_method);

                // set dirichlet boundary condition for nu_eff
                scaled_nu_eff_bc_values[patch_id][pos_in_patch] = nu_interpolated + nu_turb_interpolated / scaling_coeff;
            }
        }
    }
}



void Turbulence2D::update_k_field_vars(std::function<void ()>& k_update_func)
{
    // update k
    k_update_func();

    // enforce minimum value
    enforce_min_field_value(k_field, EPSILON);

    // no dependent fields to update
}

void Turbulence2D::update_epsilon_field_vars(std::function<void ()>& epsilon_update_func)
{
    // update epsilon
    epsilon_update_func();

    // enforce minimum value
    enforce_min_field_value(epsilon_field, EPSILON);

    // no dependent fields to update
}



void Turbulence2D::calc_k_bcs()
{
    // loop through velocity bcs and assign k accordingly
    for(int patch_id = 0; patch_id < num_patches; patch_id++)
    {
        switch (velocity_bcs[patch_id])
        {
            case BoundaryConditions::DIRICHLET:
            {
                if(has_zero_dirichlet_velocity(patch_id))
                    k_bcs[patch_id] = BoundaryConditions::NEUMANN;  // avoid having k=0 Dirichlet bc
                else
                    k_bcs[patch_id] = BoundaryConditions::DIRICHLET;

                break;
            }

            case BoundaryConditions::NEUMANN:
            {
                k_bcs[patch_id] = BoundaryConditions::NEUMANN;
                break;
            }
            
            default:
            {
                std::cerr<<"NO PROCEDURE FOR DETERMINING k FROM VELOCITY BOUNDARY CONDITIONS ON PATCH "
                            <<patches[patch_id].get_name()<<std::endl;
                std::abort();
                break;
            }
        }
    }

    for(int patch_id = 0; patch_id < num_patches; patch_id++)
    {
        std::size_t num_patch_faces = patches[patch_id].get_face_ids().size();

        switch(k_bcs[patch_id])
        {
            case BoundaryConditions::DIRICHLET:
            {
                for(std::size_t pos_in_patch = 0; pos_in_patch < num_patch_faces; pos_in_patch++)
                {
                    /* These Dirichlet boundary conditions use I_1D, unlike
                        I_3D (simply I in the report) to set k_bc_values
                        
                        The difference is that I_1D requires 1.5 * temp * temp
                            instead of 0.5 * temp * temp
                    */
                    double I_1D = 0.05;    // turbulence intensity estimate
                    double temp = velocity_bc_values[patch_id][pos_in_patch].calc_modulus() * I_1D;

                    // k is determined via formula
                    k_bc_values[patch_id][pos_in_patch] = 1.5 * temp * temp;
                }

                break;
            }

            case BoundaryConditions::NEUMANN:
            {
                for(std::size_t pos_in_patch = 0; pos_in_patch < num_patch_faces; pos_in_patch++)
                {
                    k_bc_values[patch_id][pos_in_patch] = 0;
                }

                break;
            }
            
            default:
            {
                std::cerr<<"BOUNDARY CONDITION NOT IMPLEMENTED FOR K"<<std::endl;
                std::abort();
                break;
            }
        }
    }
}

void Turbulence2D::calc_epsilon_bcs()
{
    // inherit bcs from k
    epsilon_bcs = k_bcs;

    for(int patch_id = 0; patch_id < num_patches; patch_id++)
    {
        std::size_t num_patch_faces = patches[patch_id].get_face_ids().size();

        switch(epsilon_bcs[patch_id])
        {
            case BoundaryConditions::DIRICHLET:
            {
                double l = 0.1 * char_length;  // taken as 10% of pipe diameter

                for(std::size_t pos_in_patch = 0; pos_in_patch < num_patch_faces; pos_in_patch++)
                {
                    // epsilon is determined via formula
                    epsilon_bc_values[patch_id][pos_in_patch] = C_nu * pow(k_bc_values[patch_id][pos_in_patch], 1.5) / l;
                }

                break;
            }

            case BoundaryConditions::NEUMANN:
            {
                for(std::size_t pos_in_patch = 0; pos_in_patch < num_patch_faces; pos_in_patch++)
                {
                    epsilon_bc_values[patch_id][pos_in_patch] = 0;
                }

                break;
            }
            
            default:
            {
                std::cerr<<"BOUNDARY CONDITION NOT IMPLEMENTED FOR EPSILON"<<std::endl;
                std::abort();
                break;
            }
        }
    }
}

bool Turbulence2D::has_zero_dirichlet_velocity(int patch_id)
{
    std::size_t num_patch_faces = patches[patch_id].get_face_ids().size();

    for(std::size_t pos_in_patch = 0; pos_in_patch < num_patch_faces; pos_in_patch++)
    {
        if(velocity_bc_values[patch_id][pos_in_patch].calc_modulus() == 0)
            return true;
    }
    return false;
}

double Turbulence2D::calc_field_average(std::vector<double>& field)
{
    double average = 0;

    for(double value: field)
        average += value;

    average /= field.size();

    return average;
}



void Turbulence2D::enforce_min_field_value(std::vector<double>& field, const double min_value, const bool raise_warning)
{
    // limit field to the zero tolerance value if less
    for(int cell_id = 0; cell_id < num_cells; cell_id++)
    {
        if(field[cell_id] < min_value)
        {
            if(raise_warning)
                std::cerr<<"WARNING: FIELD WAS CALCULATED AS LESS THAN "<<min_value<<" AT CELL "<<cell_id<<"!"<<std::endl;

            field[cell_id] = min_value;
        }
    }
}



void Turbulence2D::apply_wall_function(const std::function<void (int, std::size_t, int, double, double)>& wall_func,
                                        const std::function<void ()>& preset_func)
{
    // preparation
    if(preset_func)
        preset_func();

    // loop through wall faces and apply wall_func to all
    // cells with centres sufficiently close to the wall
    for(int patch_id = 0; patch_id < num_patches; patch_id++)
    {
        if(patches[patch_id].get_name().find("Wall") != std::string::npos)   // if this is a wall
        {
            std::vector<int> patch_face_ids = patches[patch_id].get_face_ids();
            std::size_t num_patch_faces = patch_face_ids.size();

            for(std::size_t pos_in_patch = 0; pos_in_patch < num_patch_faces; pos_in_patch++)
            {
                // get owner and face data
                int face_id = patch_face_ids[pos_in_patch];
                int face_owner_id = faces[face_id].get_owner_id();
                Point owner_centroid = cells[face_owner_id].get_centroid();
                Point face_centroid = faces[face_id].get_centroid();
                double dist_owner_face = Point::distance(owner_centroid, face_centroid).calc_modulus();

                // calculate normalised distance;
                // NOTE: using old value of k_field, as intended
                double norm_dist_owner_face = pow(C_nu, 0.25) * sqrt(k_field[face_owner_id]) / nu_field[face_owner_id]
                                                * dist_owner_face;

                if(norm_dist_owner_face <= inertial_sublayer_min)
                {
                    std::cerr<<"ERROR: CELL "<<face_owner_id<<" IS IN THE VISCOUS SUBLAYER OF "
                                <<patches[patch_id].get_name()<<" (FACE "<<face_id<<")!"<<std::endl;
                    std::abort();
                }
                // else if(norm_dist_owner_face >= inertial_sublayer_max)
                // {
                //     std::cerr<<"ERROR: CELL "<<face_owner_id<<" IS TOO FAR AWAY FROM THE INERTIAL SUBLAYER OF "
                //                 <<patches[patch_id].get_name()<<" (FACE "<<face_id<<")!"<<std::endl;
                //     std::cout<<"norm_dist_owner_face: "<<norm_dist_owner_face<<std::endl;
                //     std::cout<<"dist_owner_face: "<<dist_owner_face<<std::endl;
                //     std::cout<<"inertial_sublayer_max: "<<inertial_sublayer_max<<std::endl;
                //     std::abort();
                // }
                else
                {
                    wall_func(patch_id, pos_in_patch, face_owner_id, dist_owner_face, norm_dist_owner_face);
                }
            }
        }
    }
}



void Turbulence2D::update_residuals(double& init_residual, std::array<double, RESIDUAL_NO>& residuals, LinearSystem& system, std::vector<double>& field)
{
    // update previous residuals
    for(size_t i=0; i<RESIDUAL_NO-1; i++)
        residuals[i] = residuals[i+1];

    // save current residual
    residuals[RESIDUAL_NO-1] = system.calc_residual(field);

    // update initial residual
    if(init_residual == 0)
        init_residual = residuals[RESIDUAL_NO-1];
}

bool Turbulence2D::check_outer_convergence(double& init_residual, std::array<double, RESIDUAL_NO>& residuals)
{
    // custom convergence check

    bool is_convergent = true;

    for(size_t i=0; i<RESIDUAL_NO-1; i++)
    {
        is_convergent = is_convergent && ( fabs(residuals[i] - residuals[i+1]) / init_residual <= OUTER_REL_TOL);
    }

    return is_convergent;
}