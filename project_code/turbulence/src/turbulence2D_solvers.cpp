#include <cmath>
#include <functional>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "turbulence2D.hpp"

// public

void Turbulence2D::solve_k_epsilon()
{
    // declarations of fields from previous time step
    std::vector<Point> velocity_field_old;
    std::vector<double> k_field_old;
    std::vector<double> epsilon_field_old;

    // declaration of the pressure field prior to solving the pressure system
    std::vector<double> pressure_field_pre;

    // declaration of auxilliary fields
    std::vector<double> G_field;
    std::unordered_map<int, double> epsilon_wall_field;

    // declarations of update functions
    std::function<void ()> velocity_update_func;
    std::function<void ()> pressure_update_func;

    // initialise epsilon matrix const row functions
    std::function<void (double, int)> apply_overload_row = [this](double epsilon_wall, int face_owner_id)
        {
            epsilon_system.overload_row(epsilon_wall, face_owner_id);
        };

    std::function<void (double, int)> apply_set_const_row = [this](double epsilon_wall, int face_owner_id)
        {
            std::vector<int> neighbour_ids = get_neighbouring_cell_ids(face_owner_id);
            epsilon_system.set_const_row(epsilon_wall, face_owner_id, neighbour_ids);
        };

    // initialisations of wall functions and related preset functions
    std::function<void ()> preset_function_G = [this, &G_field]
        {
            // preset all G values to 0 in wall patches
            for(int patch_id = 0; patch_id < num_patches; patch_id++)
            {
                if(patches[patch_id].get_name().find("Wall") != std::string::npos)   // if this is a wall
                {
                    for(int face_id: patches[patch_id].get_face_ids())  // loop through the faces
                    {
                        // get face owner centroid
                        int face_owner_id = faces[face_id].get_owner_id();

                        // set G field to 0
                        G_field[face_owner_id] = 0;
                    }
                }
            }
        };

    std::function<void (int, std::size_t, int, double, double)> wall_function_G =
        [this, &G_field] (int, std::size_t, int face_owner_id, double dist_owner_face, double )
        {
            // calculate the gradient at the wall - assume linear gradient
            // grad_(normal_vector) (velocity) = (velocity[cell_centre] - 0) / |normal_vector|
            Point grad_velocity = velocity_field[face_owner_id];
            grad_velocity.scale_by(1./dist_owner_face);

            // G field
            // the result is additive
            G_field[face_owner_id] += wall_cell_weights[face_owner_id]
                                        * (nu_field[face_owner_id] + nu_turb_field[face_owner_id])
                                        * grad_velocity.calc_modulus() * pow(C_nu, 0.25) * sqrt(k_field[face_owner_id])
                                        / (kappa * dist_owner_face);
        };

    std::function<void ()> preset_function_epsilon = [this, &epsilon_wall_field]
        {
            // preset all epsilon values to 0 in wall-adjacent cells
            for(auto member_map: wall_cell_weights)
                epsilon_wall_field[member_map.first] = 0;
        };

    std::function<void (int, std::size_t, int, double, double)> wall_function_epsilon =
        [this, &epsilon_wall_field] (int, std::size_t, int face_owner_id, double dist_owner_face, double )
        {
            // calculate epsilon_wall_field additively
            epsilon_wall_field[face_owner_id] += wall_cell_weights[face_owner_id]
                                                    * pow(C_nu, 0.75) * pow(k_field[face_owner_id], 1.5)
                                                    / (kappa * dist_owner_face);
        };

    std::function<void (int, std::size_t, int, double, double)> wall_function_nu_turb = 
        [this] (int patch_id, std::size_t pos_in_patch, int face_owner_id, double , double norm_dist_owner_face)
        {
            // streamline in future implementations - needs to be done only once as nu_field is constant throughout the flow

            // nu wall
            nu_turb_bcs[patch_id] = BoundaryConditions::DIRICHLET;
            nu_turb_bc_values[patch_id][pos_in_patch] = nu_field[face_owner_id] * (
                                                norm_dist_owner_face * kappa / (log(E * norm_dist_owner_face)) - 1
                                            );
        };

    // defining helper functions
    std::function<void ()> pressure_post_processing = [this]()
        {
            // resolve Neumann boundary condition if needed
            std::vector<BoundaryConditions> neumann_bcs(num_patches, BoundaryConditions::NEUMANN);
            if(pressure_bcs == neumann_bcs)
            {
                int fixed_cell_id = 0;
                std::vector<int> neighbour_ids = get_neighbouring_cell_ids(fixed_cell_id);

                pressure_system.set_const_row(P_DIRICHLET(), fixed_cell_id, neighbour_ids);  // reslove by using default dirichlet value
            }
        };
    
    // initial residual declarations
    std::array<double, N_DIMS> init_velocity_residuals;
    double init_pressure_residual;
    double init_epsilon_residual;
    double init_k_residual;

    // residual declarations
    std::array<std::array<double, RESIDUAL_NO>, N_DIMS> velocity_residuals;
    std::array<double, RESIDUAL_NO> pressure_residuals;
    std::array<double, RESIDUAL_NO> epsilon_residuals;
    std::array<double, RESIDUAL_NO> k_residuals;

    // other initialisations
    std::ofstream output_file("../data/turbulence_end.dat");
    bool corr_flux_check_cells = true;

    // relative tolerance for pressure
    double pressure_tol = 1e-5;

    // numbers of operations in each assembly function
    int num_operations_velocity = 3;    // 3 off-diagonal sweeps:   1 convection and 2 diffusions
    int num_operations_pressure = 1;    // 1 off-diagonal sweep:    1 diffusion
    int num_operations_k = 2;           // 2 off-diagonal sweeps:   1 convection and 1 diffusion
    int num_operations_epsilon = 2;     // 2 off-diagonal sweeps:   1 convection and 1 diffusion



    // write initial conditions
    write_output(output_file);

    for(int i=1; i<=N_TIME_STEPS; i++)
    {
        std::cout<<"################### t="<<i*D_T<<" ###################"<<std::endl;

        // initialise old fields from the previous iteration (or initial guesses)
        velocity_field_old = velocity_field;
        k_field_old = k_field;
        epsilon_field_old = epsilon_field;

        // initialise convergence variable
        bool convergence_reached = false;

        // set initial residuals to 0 for consistency
        for(Dimensions dim: ALL_DIMS)
            init_velocity_residuals[dim] = 0;

        init_pressure_residual = 0;
        init_epsilon_residual = 0;
        init_k_residual = 0;

        // set residuals to 0
        for(size_t i=0; i<RESIDUAL_NO; i++)
        {
            for(Dimensions dim: ALL_DIMS)
                velocity_residuals[dim][i] = 0;

            pressure_residuals[i] = 0;
            epsilon_residuals[i] = 0;
            k_residuals[i] = 0;
        }

        // iterate until convergence for t_i is reached or until N_ITERATIONS
        for(int n=0; !convergence_reached && n<N_ITERATIONS; n++)
        {
            std::cout<<"------------------- n="<<n<<" -------------------"<<std::endl;

            //  ------------------- PRE-SET VARIABLES -------------------

            // clear matrix coefficients
            for(Dimensions dim: ALL_DIMS)
                velocity_systems[dim].clear();
            
            pressure_system.clear();
            k_system.clear();
            epsilon_system.clear();

            //  ------------------- ASSEMBLE VELOCITY SYSTEMS WITH GRAD P SOURCE -------------------

            // assemble velocity systems
            std::function<void ()> velocity_assembly =
                [this, &velocity_field_old] () 
                {
                    assemble_velocity_systems(velocity_field_old);

                    // add -grad(p) to rhs
                    for(Dimensions dim: ALL_DIMS)
                    {
                        // extract grad(p) [dim]
                        std::vector<double> grad_p_dim_field;
                        extract_dim_field(grad_p_field, grad_p_dim_field, dim);

                        calc_collocated_source_coeffs(velocity_systems[dim], -1, grad_p_dim_field);
                    }
                };
            assemble_matrix(velocity_systems, num_operations_velocity, velocity_assembly);

            // apply implicit underrelaxation
            for(Dimensions dim: ALL_DIMS)
            {
                std::vector<double> velocity_dim_field;
                extract_dim_field(velocity_field, velocity_dim_field, dim);

                underrelax_system_implicit(velocity_systems[dim], velocity_dim_field, velocity_underrel_coeff);
            }

            // ------------------- SOLVE FOR VELOCITY -------------------

            // update velocity residuals
            for(Dimensions dim: ALL_DIMS)
            {
                // extract velocity field
                std::vector<double> velocity_dim_field;
                extract_dim_field(velocity_field, velocity_dim_field, dim);

                // update init_residual and residuals
                update_residuals(init_velocity_residuals[dim], velocity_residuals[dim], velocity_systems[dim], velocity_dim_field);
            }

            // perform the solve function
            velocity_update_func = [this]() { solve_velocity(); };
            update_velocity_field_vars(velocity_update_func);

            // ------------------- REMOVE grad(p) FROM velocity_systems DISCRETISATION -------------------

            // add +grad p to rhs, removing it from the equation
            for(Dimensions dim: ALL_DIMS)
            {
                // extract grad(p) [dim]
                std::vector<double> grad_p_dim_field;
                extract_dim_field(grad_p_field, grad_p_dim_field, dim);

                calc_collocated_source_coeffs(velocity_systems[dim], 1, grad_p_dim_field);  // +1 for eq_sign
            }

            // ------------------- GET DIFFUSIVITY FIELD FOR PRESSURE -------------------
            std::vector<double> diffusivity_field_p = calc_pressure_diffusivity_field();

            // ------------------- ADJUST VELOCITY AND FLUXES TO PSEUDOVELOCITY -------------------
            velocity_update_func = [this, &diffusivity_field_p]() { calc_pseudovelocity_field(diffusivity_field_p); };

            // below is bookkeeping work; currently unphysical
            update_velocity_field_vars(velocity_update_func);
            calc_fluxes(velocity_intpol_method);

            // ------------------- SOLVE FOR PRESSURE -------------------

            // assemble pressure matrix
            std::function<void ()> pressure_assembly =
                [this, &diffusivity_field_p]()
                {
                    PVCoupling2D::assemble_pressure_system(diffusivity_field_p,
                                                            diffusivity_field_p_bcs,
                                                            diffusivity_field_p_bc_values);
                };
            assemble_matrix(pressure_system, num_operations_pressure, pressure_assembly, pressure_post_processing);

            // update pressure residuals
            update_residuals(init_pressure_residual, pressure_residuals, pressure_system, pressure_field);

            // book keep for the next iteration
            pressure_field_pre = pressure_field;

            // solve for pressure
            pressure_update_func = [this, pressure_tol]() { solve_pressure(pressure_tol); };

            update_pressure_field_vars(pressure_update_func);

            // ------------------- UPDATE FLUXES TO SATISFY INCOMPRESSIBILITY -------------------

            // correct the face flux: f_corr = f_pre - (1/ap)_f sf * (grad p)f
            update_fluxes_with_pressure(diffusivity_field_p, diffusivity_field_p_bcs, diffusivity_field_p_bc_values);

            check_fluxes(corr_flux_check_cells);

            // ------------------- UNDERRELAX PRESSURE FOR BETTER PHYSICAL VALUE -------------------

            pressure_update_func = [this, &pressure_field_pre]() {
                underrelax_field_explicit(pressure_field, pressure_field_pre, pressure_underrel_coeff);
            };

            update_pressure_field_vars(pressure_update_func);

            // ------------------- CORRECT VELOCITIES - WILL MATCH FLUXES AT CONVERGENCE -------------------

            // correct cell centre velocities with underrelaxed pressure: u_corr = u_pre - 1/ap grad p
            velocity_update_func = [this, &diffusivity_field_p]() { update_velocity_field_with_grad_p(diffusivity_field_p); };
            update_velocity_field_vars(velocity_update_func);

            // ------------------- CALCULATE G FIELD -------------------

            calc_G_field_incompressible(G_field);

            // apply wall function
            apply_wall_function(wall_function_G, preset_function_G);

            // ------------------- ASSEMBLE EPSILON SYSTEM -------------------

            // assemble epsilon matrix
            std::function<void ()> epsilon_assembly =
                [this, &epsilon_field_old, &G_field]()
                {
                    assemble_epsilon_system(epsilon_field_old, G_field);
                };
            assemble_matrix(epsilon_system, num_operations_epsilon, epsilon_assembly);

            // apply implicit underrelaxation to epsilon
            underrelax_system_implicit(epsilon_system, epsilon_field, epsilon_underrel_coeff);

            // apply wall function to epsilon_wall_field
            apply_wall_function(wall_function_epsilon, preset_function_epsilon);

            // propagate epsilon_wall_field throught the epsilon system
            for(auto member_map: epsilon_wall_field)
                apply_set_const_row(member_map.second, member_map.first);

            // ------------------- SOLVE FOR EPSILON -------------------

            // update epsilon residuals
            update_residuals(init_epsilon_residual, epsilon_residuals, epsilon_system, epsilon_field);

            // solve for epsilon
            std::function<void ()> epsilon_update_func = [this]() { solve_epsilon(); };
            update_epsilon_field_vars(epsilon_update_func);                 // min_value enforced in update method

            // ------------------- ASSEMBLE K SYSTEM -------------------

            // assemble k matrix
            std::function<void ()> k_assembly =
                [this, &k_field_old, &G_field]()
                {
                    assemble_k_system(k_field_old, G_field);
                };
            assemble_matrix(k_system, num_operations_k, k_assembly);

            // apply implicit underrelaxation to k
            underrelax_system_implicit(k_system, k_field, k_underrel_coeff);

            // ------------------- SOLVE FOR K -------------------

            // update k residuals
            update_residuals(init_k_residual, k_residuals, k_system, k_field);

            // solve for k
            std::function<void ()> k_update_func = [this]() { solve_k(); };
            update_k_field_vars(k_update_func);                             // min_value enforced in update method

            // ------------------- UPDATE TURBULENT VISCOSITY -------------------

            calc_nu_turb_field();

            // apply wall function
            apply_wall_function(wall_function_nu_turb);

            // ------------------- DETERMINE IF CONVERGENCE IS REACHED -------------------

            convergence_reached = check_outer_convergence(init_velocity_residuals[X], velocity_residuals[X])
                                    && check_outer_convergence(init_velocity_residuals[Y], velocity_residuals[Y])
                                    && check_outer_convergence(init_pressure_residual, pressure_residuals)
                                    && check_outer_convergence(init_epsilon_residual, epsilon_residuals)
                                    && check_outer_convergence(init_k_residual, k_residuals);
        }

        // ------------------- WRITE OUTPUT -------------------
        write_output(output_file);
    }
}



// protected

void Turbulence2D::solve_k(double rel_tol)
{
    std::cout<<"----- k system -----"<<std::endl;

    k_system.solve_Gauss_Seidel(k_field, rel_tol);
}

void Turbulence2D::solve_epsilon(double rel_tol)
{
    std::cout<<"----- epsilon system -----"<<std::endl;

    epsilon_system.solve_Gauss_Seidel(epsilon_field, rel_tol);
}