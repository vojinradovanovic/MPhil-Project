#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "pvCoupling2D.hpp"

// public

void PVCoupling2D::solve_simple()
{
    // solve with SIMPLE algorithm

    // initialisations
    std::vector<double> pressure_field_old = pressure_field;
    std::ofstream output_file("../../data/simple_end.dat");
    // std::ofstream output_file_velocity("../../data/simple_velocity.dat");
    // std::ofstream output_v_matrices("../../data/v_matrix_data.dat");
    // std::ofstream output_v_vectors("../../data/v_vector_data.dat");
    // std::ofstream output_p_matrix("../../data/p_matrix_data.dat");
    // std::ofstream output_p_vectors("../../data/p_vector_data.dat");

    // defining parameters
    bool pre_flux_check_cells = false;
    bool corr_flux_check_cells = true;

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

    // numbers of operations in each assembly function
    int num_operations_velocity = 2;    // 2 off-diagonal sweeps:   1 convection and 1 diffusion
    int num_operations_pressure = 1;    // 1 off-diagonal sweep:    1 diffusion

    for(int n=0; n<N_ITERATIONS; n++)
    {
        std::cout<<"------------------- n="<<n<<" -------------------"<<std::endl;

        //  ------------------- PRE-SET VARIABLES -------------------
        // clear matrix coefficients
        for(Dimensions dim: ALL_DIMS)
            velocity_systems[dim].clear();
        
        pressure_system.clear();

        //  ------------------- DISCRETISATION -------------------

        // CONVECTION           DIFFUSION
        // div(u_new*u_old) - div(nu * grad(u)) = - grad(p), p = P/rho is the kinematic pressure

        // assemble velocity systems
            std::function<void ()> velocity_assembly =
                [this] () 
                {
                    assemble_velocity_systems();

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

        // result is a sparse matrix A, s.t. A*u=b

        // ------------------- SOLVE FOR VELOCITY -------------------
        
        // std::function<void ()> velocity_update_func =
        //     [this, &output_v_matrices, &output_v_vectors]() { solve_velocity(output_v_matrices, output_v_vectors); };
        std::function<void ()> velocity_update_func =
            [this]() { solve_velocity(); };

        update_velocity_field_vars(velocity_update_func);

        // // write velocity output
        // write_output(output_file_velocity);

        // face flux from above velocity field, s_f*u_f with linear interpolation on u
        std::vector<double> old_fluxes = fluxes;
        calc_fluxes(velocity_intpol_method);    // only calculated for output, need to use old_fluxes to
        check_fluxes(pre_flux_check_cells);     // assemble H operator later
        fluxes = old_fluxes;    // restore old flux values

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
                assemble_pressure_system(diffusivity_field_p, diffusivity_field_p_bcs, diffusivity_field_p_bc_values);
            };
        assemble_matrix(pressure_system, num_operations_pressure, pressure_assembly, pressure_post_processing);

        // book keep for the next iteration
        pressure_field_old = pressure_field;  

        // solve for pressure
        // std::function<void ()> pressure_update_func =
        //     [this, &output_p_matrix, &output_p_vectors]() { solve_pressure(output_p_matrix, output_p_vectors); };
        std::function<void ()> pressure_update_func =
            [this]() { solve_pressure(); };

        update_pressure_field_vars(pressure_update_func);

        // ------------------- UPDATE FLUXES TO SATISFY INCOMPRESSIBILITY -------------------

        // correct the face flux: f_corr = f_pre - (1/ap)_f sf * (grad p)f
        update_fluxes_with_pressure(diffusivity_field_p, diffusivity_field_p_bcs, diffusivity_field_p_bc_values);

        check_fluxes(corr_flux_check_cells);

        // ------------------- UNDERRELAX PRESSURE FOR BETTER PHYSICAL VALUE -------------------
        pressure_update_func = [this, &pressure_field_old]() {
            underrelax_field_explicit(pressure_field, pressure_field_old, pressure_underrel_coeff);
        };

        update_pressure_field_vars(pressure_update_func);

        // ------------------- CORRECT VELOCITIES - WILL MATCH FLUXES AT CONVERGENCE -------------------

        // correct cell centre velocities with underrelaxed pressure: u_corr = u_pre - 1/ap grad p
        velocity_update_func = [this, &diffusivity_field_p]() { update_velocity_field_with_grad_p(diffusivity_field_p); };
        update_velocity_field_vars(velocity_update_func);

        // ------------------- WRITE OUTPUT -------------------
        write_output(output_file);
    }
}

void PVCoupling2D::solve_diffusion()
{
    // solve Poisson equation

    // solving velocity_systems[X]

    //  ------------------- PRE-SET VARIABLES -------------------
    // clear matrix coefficients
    for(Dimensions dim: ALL_DIMS)
        velocity_systems[dim].clear();

    //  ------------------- DISCRETISATION -------------------

    // DIFFUSION
    // div(nu * grad(v_x)) = RHS

    //  ------------------- RHS = -2*pi^2*sin(pi*10*x)*sin(pi*10*y) -------------------

    double PI = 4 * atan(1);

    for(int i=0; i<num_cells; i++)
    {
        // get cell centroid coordinates
        std::array<double, N_DIMS> centroid_coords = cells[i].get_centroid().get_coordinates();

        // rhs_term = source_at_cell * cell_volume
        double rhs_term = -2*PI*PI
                            *sin(10*PI*centroid_coords[X])*sin(10*PI*centroid_coords[Y])
                            * cells[i].get_volume();

        velocity_systems[X].add_to_rhs_vector(rhs_term, i);
    }

    //  ------------------- DIFFUSION ON V_X (V_Y will also be computed but not acted on) -------------------

    // assemble epsilon matrix
    std::function<void ()> velocity_x_assembly =
        [this]()
        {
            // extract velocity boundary conditions
            std::vector<std::vector<double>> velocity_x_bc_values;
            extract_dim_field(velocity_bc_values, velocity_x_bc_values, X);

            // use field of ones
            nu_field = std::vector<double> (num_cells, 1);
            calc_diffusion_coeffs(velocity_systems[X], 1, velocity_bcs, velocity_x_bc_values,
                                    nu_field, nu_bcs, nu_bc_values, diffusivity_intpol_method);
        };
    assemble_matrix(velocity_systems[X], 1, velocity_x_assembly);


    // result is a sparse matrix A, s.t. A*u=b

    // ------------------- SOLVE FOR V_X -------------------

    std::function<void ()> velocity_update_func = [this]() { solve_velocity_x_diffusion(); };
    update_velocity_field_vars(velocity_update_func);

    // write output
    std::ofstream output_file("../../data/diffusion.dat");
    write_output(output_file);
}

void PVCoupling2D::solve_convection()
{
    // solve (material derivative of velocity) = 0

    // defining parameters
    bool pre_flux_check_cells = false;
    std::ofstream output_file_velocity("../../data/convection.dat");
    // std::ofstream output_v_matrices("../../data/v_matrix_data.dat");
    // std::ofstream output_v_vectors("../../data/v_vector_data.dat");

    // write initial conditions
    write_output(output_file_velocity);

    for(int n=0; n<N_TIME_STEPS; n++)
    {
        std::cout<<"------------------- n="<<n<<" -------------------"<<std::endl;

        //  ------------------- PRE-SET VARIABLES -------------------
        // clear matrix coefficients
        for(Dimensions dim: ALL_DIMS)
            velocity_systems[dim].clear();

        //  ------------------- DISCRETISATION -------------------

        // INERTIA  CONVECTION
        // du/dt +  div(u_new*u_old) = 0

        // assemble velocity systems
        std::function<void ()> velocity_assembly =
            [this] () 
            {
                for(Dimensions dim: ALL_DIMS)
                {
                    // set up scalar field before calculating coefficients
                    std::vector<double> velocity_dim_field;
                    extract_dim_field(velocity_field, velocity_dim_field, dim);

                    std::vector<std::vector<double>> velocity_dim_bc_values;
                    extract_dim_field(velocity_bc_values, velocity_dim_bc_values, dim);

                    calc_time_der_coeffs(velocity_systems[dim], 1, velocity_dim_field);

                    // assemble convection matrix coefficients with the upwind differencing method
                    calc_convection_coeffs(velocity_systems[dim], 1,
                                            velocity_bcs, velocity_dim_bc_values,
                                            convection_intpol_method);      // +1 as the equation is +div(u_new*u_old)
                }
            };
        assemble_matrix(velocity_systems, 1, velocity_assembly);

        // result is a sparse matrix A, s.t. A*u=b

        // ------------------- SOLVE FOR VELOCITY -------------------

        // std::function<void ()> velocity_update_func =
        //     [this, &output_v_matrices, &output_v_vectors]() {
        //         solve_velocity(output_v_matrices, output_v_vectors);
        //     };
        std::function<void ()> velocity_update_func =
            [this]() {
                solve_velocity();
            };

        update_velocity_field_vars(velocity_update_func);

        // face flux from above velocity field, s_f*u_f with linear interpolation on u
        calc_fluxes(velocity_intpol_method);
        check_fluxes(pre_flux_check_cells);

        write_output(output_file_velocity);

        // ------------------- ASSESS THAT DIAGONAL COEFFICIENTS ARE THE SAME -------------------
        std::array<std::vector<double>, N_DIMS> diag_coeffs;
        for(Dimensions dim: ALL_DIMS)
            diag_coeffs[dim] = velocity_systems[dim].get_diagonal_coeffs();

        assert(diag_coeffs[X] == diag_coeffs[Y]);
    }
}

void PVCoupling2D::solve_advection_diffusion()
{
    // solve advection diffusion equation

    // defining parameters
    bool pre_flux_check_cells = false;
    std::ofstream output_file_velocity("../../data/advection_diffusion.dat");
    // std::ofstream output_v_matrices("../../data/v_matrix_data.dat");
    // std::ofstream output_v_vectors("../../data/v_vector_data.dat");

    for(int n=0; n<1; n++)
    {
        std::cout<<"------------------- n="<<n<<" -------------------"<<std::endl;

        //  ------------------- PRE-SET VARIABLES -------------------
        // clear matrix coefficients
        for(Dimensions dim: ALL_DIMS)
            velocity_systems[dim].clear();

        //  ------------------- DISCRETISATION -------------------

        // assemble velocity systems
        std::function<void ()> velocity_assembly =
            [this] () 
            {
                for(Dimensions dim: ALL_DIMS)
                {
                    // get field values for velocity[dim]
                    std::vector<double> velocity_dim_field;
                    extract_dim_field(velocity_field, velocity_dim_field, dim);

                    std::vector<std::vector<double>> velocity_dim_bc_values;
                    extract_dim_field(velocity_bc_values, velocity_dim_bc_values, dim);

                    // discretise equation:
                    // CONVECTION           DIFFUSION
                    // div(u_new*u_old) - div(nu * grad(u)) = 0

                    // steady state: d/dt = 0

                    // assemble convection matrix coefficients with the upwind differencing method
                    calc_convection_coeffs(velocity_systems[dim], 1,
                                            velocity_bcs, velocity_dim_bc_values,
                                            convection_intpol_method);    // +1 as the equation is +div(u_new*u_old)

                    // assemble diffusion matrix in standard way
                    calc_diffusion_coeffs(velocity_systems[dim], -1,
                                            velocity_bcs, velocity_dim_bc_values,
                                            nu_field, nu_bcs, nu_bc_values, diffusivity_intpol_method);  // -1 as the equation is -div(nu*grad(u))
                }
            };
        assemble_matrix(velocity_systems, 2, velocity_assembly);

        // ------------------- SOLVE FOR VELOCITY -------------------

        // std::function<void ()> velocity_update_func =
        //     [this, &output_v_matrices, &output_v_vectors]() {
        //         solve_velocity(output_v_matrices, output_v_vectors);
        //     };
        std::function<void ()> velocity_update_func =
            [this]() {
                solve_velocity();
            };

        update_velocity_field_vars(velocity_update_func);

        // face flux from above velocity field, s_f*u_f with linear interpolation on u
        calc_fluxes(velocity_intpol_method);
        check_fluxes(pre_flux_check_cells);

        // write output
        write_output(output_file_velocity);

        // ------------------- ASSESS THAT DIAGONAL COEFFICIENTS ARE THE SAME -------------------
        std::array<std::vector<double>, N_DIMS> diag_coeffs;
        for(Dimensions dim: ALL_DIMS)
            diag_coeffs[dim] = velocity_systems[dim].get_diagonal_coeffs();

        assert(diag_coeffs[X] == diag_coeffs[Y]);
    }
}



// protected

void PVCoupling2D::solve_velocity(double tol)
{
    std::array<std::string, N_DIMS> dim_names {"X", "Y"};

    // solve each velocity_systems[dim]
    for(Dimensions dim: ALL_DIMS)
    {
        std::cout<<"----- velocity "<<dim_names[dim]<<" system -----"<<std::endl;

        // set up initial field for solver
        std::vector<double> velocity_dim_field;
        extract_dim_field(velocity_field, velocity_dim_field, dim);

        // solve for velocity_dim
        velocity_systems[dim].solve_Gauss_Seidel(velocity_dim_field, tol);

        // update velocity
        for(int cell_id = 0; cell_id < num_cells; cell_id++)
        {
            std::array<double, N_DIMS> velocity_coords = velocity_field[cell_id].get_coordinates();
            velocity_coords[dim] = velocity_dim_field[cell_id];

            velocity_field[cell_id].set_coordinates(velocity_coords);
        }
    }
}

void PVCoupling2D::solve_velocity(std::ofstream& output_v_matrices, std::ofstream& output_v_vectors, double tol)
{
    // initialisations
    std::array<std::string, N_DIMS> dim_names {"X", "Y"};

    // set up residual fields
    std::vector<std::vector<double>> velocity_residual_fields_before(N_DIMS, std::vector<double>(num_cells)),
                                    velocity_residual_fields_after(N_DIMS, std::vector<double>(num_cells));

    // solve each velocity_systems[dim]
    for(Dimensions dim: ALL_DIMS)
    {
        // set up initial iteration for solver
        std::vector<double> velocity_dim_field;
        extract_dim_field(velocity_field, velocity_dim_field, dim);

        // set up residual field before
        for(int cell_id = 0; cell_id < num_cells; cell_id++)
        {
            double residual_dim = velocity_systems[dim].calc_residual(velocity_dim_field, cell_id);
            velocity_residual_fields_before[dim][cell_id] = residual_dim;
        }

        std::cout<<"----- velocity "<<dim_names[dim]<<" system -----"<<std::endl;

        // solve for velocity
        velocity_systems[dim].solve_Gauss_Seidel(velocity_dim_field, tol);

        // update velocity
        for(int cell_id = 0; cell_id < num_cells; cell_id++)
        {
            std::array<double, N_DIMS> velocity_coords = velocity_field[cell_id].get_coordinates();
            velocity_coords[dim] = velocity_dim_field[cell_id];

            velocity_field[cell_id].set_coordinates(velocity_coords);
        }

        // set up residual field after
        for(int cell_id = 0; cell_id < num_cells; cell_id++)
        {
            double residual_dim = velocity_systems[dim].calc_residual(velocity_dim_field, cell_id);
            velocity_residual_fields_after[dim][cell_id] = residual_dim;
        }
    }

    // set up vector for writing function
    std::vector<LinearSystem> vector_of_velocity_systems {velocity_systems[X], velocity_systems[Y]};

    // write to output files
    write_linear_systems_data(output_v_matrices, output_v_vectors, vector_of_velocity_systems,
                                velocity_residual_fields_before,
                                velocity_residual_fields_after);
}



void PVCoupling2D::solve_pressure(double tol)
{
    std::cout<<"----- pressure system -----"<<std::endl;

    pressure_system.solve_Gauss_Seidel(pressure_field, tol);
}

void PVCoupling2D::solve_pressure(std::ofstream& output_p_matrix, std::ofstream& output_p_vectors, double tol)
{
    // set up residual fields
    std::vector<std::vector<double>> pressure_residual_fields_before{std::vector<double>(num_cells)},
                                    pressure_residual_fields_after{std::vector<double>(num_cells)};

    // set up residual field before
    for(int cell_id = 0; cell_id < num_cells; cell_id++)
    {
        double residual = pressure_system.calc_residual(pressure_field, cell_id);
        pressure_residual_fields_before[0][cell_id] = residual;
    }

    std::cout<<"----- pressure system -----"<<std::endl;

    pressure_system.solve_Gauss_Seidel(pressure_field, tol);

    // set up residual field after
    for(int cell_id = 0; cell_id < num_cells; cell_id++)
    {
        double residual = pressure_system.calc_residual(pressure_field, cell_id);
        pressure_residual_fields_after[0][cell_id] = residual;
    }

    // set up vector for output function
    std::vector<LinearSystem> vector_of_pressure_systems {pressure_system};

    // write to output files
    write_linear_systems_data(output_p_matrix, output_p_vectors, vector_of_pressure_systems,
                                pressure_residual_fields_before,
                                pressure_residual_fields_after);
}



void PVCoupling2D::query_cell(LinearSystem& lin_system, int cell_id)
{
    std::cout<<"###### CELL QUERY: CELL "<<cell_id<<" #####"<<std::endl;

    // get matrix coefficients
    lin_system.query_cell(cell_id);

    // get cell centre values
    std::cout<<"===== cell centre values ====="<<std::endl;

    // velocity
    std::array<double, N_DIMS> velocity_coords = velocity_field[cell_id].get_coordinates();
    std::cout<<"Velocity: ";
    for(Dimensions dim: ALL_DIMS)
        std::cout<<velocity_coords[dim]<<" ";
    std::cout<<std::endl;

    // pressure
    std::cout<<"Pressure: "<<pressure_field[cell_id]<<std::endl;

    // grad p
    std::array<double, N_DIMS> grad_p_coords = grad_p_field[cell_id].get_coordinates();
    std::cout<<"Grad p: ";
    for(Dimensions dim: ALL_DIMS)
        std::cout<<grad_p_coords[dim]<<" ";
    std::cout<<std::endl;

    // get face values
    std::cout<<"===== cell face values ====="<<std::endl;
    for(int face_id: cells[cell_id].get_face_ids())
    {
        int face_owner_id = faces[face_id].get_owner_id();
        int face_neighbour_id = faces[face_id].get_neighbour_id();

        std::cout<<"----- Face "<<face_id<<", owner "<<face_owner_id<<", neighbour "<<face_neighbour_id<<" -----"<<std::endl;

        // sf * velocity_f
        std::cout<<"Face flux: "<<fluxes[face_id]<<std::endl;

        // face pressure
        std::cout<<"Interpolated face pressure: "
                    <<interpolate_to_face(pressure_field, pressure_bcs, pressure_bc_values, face_id, pressure_intpol_method)
                    <<std::endl;
    }

    std::cout<<"######################"<<std::endl;
}



void PVCoupling2D::solve_velocity_x_diffusion(double tol)
{
    // set up initial iteration for solver
    std::vector<double> velocity_x_field;
    extract_dim_field(velocity_field, velocity_x_field, X);

    // solve for velocity
    velocity_systems[X].solve_Gauss_Seidel(velocity_x_field, tol);

    // update velocity
    for(int cell_id = 0; cell_id < num_cells; cell_id++)
    {
        std::array<double, N_DIMS> velocity_coords = velocity_field[cell_id].get_coordinates();
        velocity_coords[X] = velocity_x_field[cell_id];

        velocity_field[cell_id].set_coordinates(velocity_coords);
    }
}