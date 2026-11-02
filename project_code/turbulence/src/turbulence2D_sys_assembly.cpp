#include <array>
#include <cassert>
#include <vector>

#include "turbulence2D.hpp"

// protected

void Turbulence2D::assemble_k_system(std::vector<double>& k_field_old, std::vector<double>& G_field)
{
    // ------------------- CALCULATE PARAMETERS -------------------
    // initialisations
    std::vector<double> scaled_nu_eff_field;
    std::vector<BoundaryConditions> scaled_nu_eff_bcs;
    std::vector<std::vector<double>> scaled_nu_eff_bc_values;

    // calculate scaled nu_eff fields
    calc_scaled_nu_eff_fields(scaled_nu_eff_field, scaled_nu_eff_bcs, scaled_nu_eff_bc_values, sigma_k);

    // ------------------- ASSEMBLE K SYSTEM -------------------
    // TIME DER     CONVECTION      DIFFUSION
    // dk/dt    + div(u_bar*k) - div(nu_eff * grad(k)) = G - epsilon

    // assemble time derivative
    calc_time_der_coeffs(k_system, 1, k_field_old);

    // assemble convection matrix coefficients with the upwind differencing method
    calc_convection_coeffs(k_system, 1, k_bcs, k_bc_values, convection_intpol_method);

    // assemble diffusion matrix
    calc_diffusion_coeffs(k_system, -1,
                            k_bcs, k_bc_values,
                            scaled_nu_eff_field, scaled_nu_eff_bcs, scaled_nu_eff_bc_values,
                            diffusivity_intpol_method);

    // assemble rhs terms
    calc_k_rhs_coeffs(G_field);     // k_system is implied in the function
}

void Turbulence2D::assemble_epsilon_system(std::vector<double>& epsilon_field_old, std::vector<double>& G_field)
{
    // ------------------- CALCULATE SCALED NU_EFF FIELDS -------------------
    // initialisations
    std::vector<double> scaled_nu_eff_field;
    std::vector<BoundaryConditions> scaled_nu_eff_bcs;
    std::vector<std::vector<double>> scaled_nu_eff_bc_values;

    // calculate fields
    calc_scaled_nu_eff_fields(scaled_nu_eff_field, scaled_nu_eff_bcs, scaled_nu_eff_bc_values, sigma_epsilon);

    // ------------------- ASSEMBLE EPSILON SYSTEM -------------------
    // TIME DER     CONVECTION      DIFFUSION
    // deps/dt  + div(u_bar*eps) - div(nu_eff * grad(eps)) = C1*G*eps/k - C2*eps^2/k

    // assemble time derivative
    calc_time_der_coeffs(epsilon_system, 1, epsilon_field_old);

    // assemble convection matrix coefficients with the upwind differencing method
    calc_convection_coeffs(epsilon_system, 1, epsilon_bcs, epsilon_bc_values, convection_intpol_method);

    // assemble diffusion matrix
    calc_diffusion_coeffs(epsilon_system, -1,
                            epsilon_bcs, epsilon_bc_values,
                            scaled_nu_eff_field, scaled_nu_eff_bcs, scaled_nu_eff_bc_values,
                            diffusivity_intpol_method);

    // as the source field is of the form C1*G*eps/k - C2*eps^2/k, to boost diagonal dominance we move
    // C2*eps^2/k to the LHS, linearising eps^2 as eps_old*eps_new

    // assemble rhs coeffs
    calc_epsilon_rhs_coeffs(G_field);    // epsilon_system is implied in the function
}



void Turbulence2D::assemble_velocity_systems(std::vector<Point>& velocity_field_old)
{
    // ------------------- CALCULATE HELPER FIELDS -------------------
    // declarations
    std::vector<double> nu_eff_field;
    std::vector<BoundaryConditions> nu_eff_bcs;
    std::vector<std::vector<double>> nu_eff_bc_values;
    std::vector<Point> grad_k_field;
    std::vector<Tensor<Point>> grad_velocity_transpose_field;

    // calculate scaled nu_eff fields
    calc_scaled_nu_eff_fields(nu_eff_field, nu_eff_bcs, nu_eff_bc_values, 1);   // no scaling: nu_eff=nu+nu_turb

    // initialise grad(k) field
    grad_k_field.resize(num_cells);
    calc_grad_field_of_scalar(k_field, grad_k_field, k_bcs, k_bc_values, k_intpol_method);

    // scale grad(k) by 2/3
    for(int cell_id = 0; cell_id < num_cells; cell_id++)
        grad_k_field[cell_id].scale_by(2./3);

    // calculate grad_velocity_field
    calc_grad_velocity_transpose_field(grad_velocity_transpose_field);

    // other initialisations
    std::vector<BoundaryConditions> grad_velocity_transpose_bcs (num_patches, BoundaryConditions::NEUMANN);
    Point zero_point;
    std::vector<std::vector<Point>> grad_velocity_transpose_bc_values(num_patches);
    InterpolationMethods grad_velocity_transpose_intpol_method = InterpolationMethods::LINEAR;

    for(int patch_id = 0; patch_id < num_patches; patch_id++)
    {
        std::size_t num_patch_faces = patches[patch_id].get_face_ids().size();
        grad_velocity_transpose_bc_values[patch_id].resize(num_patch_faces, zero_point);
    }

    // ------------------- ASSEMBLE VELOCITY SYSTEMS -------------------
    for(Dimensions dim: ALL_DIMS)
    {
        // get field values for dim
        std::vector<double> velocity_dim_field_old;
        extract_dim_field(velocity_field_old, velocity_dim_field_old, dim);

        std::vector<std::vector<double>> velocity_dim_bc_values;
        extract_dim_field(velocity_bc_values, velocity_dim_bc_values, dim);

        std::vector<double> grad_k_dim_field;
        extract_dim_field(grad_k_field, grad_k_dim_field, dim);

        std::vector<Point> reynolds_source_field;
        reynolds_source_field.resize(num_cells);
        for(int cell_id = 0; cell_id < num_cells; cell_id++)
        {
            reynolds_source_field[cell_id] = grad_velocity_transpose_field[cell_id].get_component(dim);
            reynolds_source_field[cell_id].scale_by(nu_turb_field[cell_id]);
        }

        // discretise equation:
        // TIME DER         CONVECTION           DIFFUSION
        // d(u_new)/dt + div(u_new*u_old) - div(nu_eff * grad(u_new)) = 0

        // assemble time derivative
        calc_time_der_coeffs(velocity_systems[dim], 1, velocity_dim_field_old);

        // assemble convection matrix coefficients with the upwind differencing method
        calc_convection_coeffs(velocity_systems[dim], 1,
                                velocity_bcs, velocity_dim_bc_values,
                                convection_intpol_method);    // +1 as the equation is +div(u_new*u_old)

        // assemble diffusion matrix
        calc_diffusion_coeffs(velocity_systems[dim], -1,
                                velocity_bcs, velocity_dim_bc_values,
                                nu_eff_field, nu_eff_bcs, nu_eff_bc_values,
                                diffusivity_intpol_method);  // -1 as the equation is -div(nu*grad(u))

        // assemble Reynolds stress tensor term: diffusion component and source component
        // R = 2/3 * k * I - nu_turb (grad(u) + grad(u)^T);

        // LHS: note that -grad(nu_turb*grad(u)) was already added to LHS when calculating nu_eff_field

        // RHS: grad(nu_turb * grad(u)^T) - 2/3 grad(k)

        // add the rest to source terms on RHS
        calc_divergence_source_coeffs(velocity_systems[dim], 1, reynolds_source_field,
                                        grad_velocity_transpose_bcs, grad_velocity_transpose_bc_values,
                                        grad_velocity_transpose_intpol_method);
        calc_collocated_source_coeffs(velocity_systems[dim], -1, grad_k_dim_field);

        // pressure source added outside of function
    }
}



void Turbulence2D::calc_k_rhs_coeffs(std::vector<double>& G_field)
{
    // the source is G - epsilon, reframed as G - (epsilon/k)*k, hence
    // epsilon is moved to the LHS to boost diagonal dominance

    // initialisations
    std::vector<double> k_Sp_field(num_cells, 0);

    // populate k_Sp_field
    for(int cell_id = 0; cell_id < num_cells; cell_id++)
    {
        double Sp = - epsilon_field[cell_id] / k_field[cell_id];    // note the minus sign here

        // Sp is guaranteed negative because we enforce epsilon>0 and k>0 after solving
        
        k_Sp_field[cell_id] = Sp;
    }

    // assemble RHS terms
    calc_collocated_source_coeffs(k_system, 1, G_field);
    calc_Sp_term_coeffs(k_system, k_Sp_field);
}

void Turbulence2D::calc_epsilon_rhs_coeffs(std::vector<double>& G_field)
{
    // initialisations
    std::vector<double> epsilon_source_field(num_cells, 0);
    std::vector<double> epsilon_Sp_field(num_cells, 0); // according to the source discretisation: S(field)=S_u+S_p*field
    std::array<double, 2> mult_terms;                   // used for examining source terms

    // populate epsilon_source_field (and diagonal_coeffs if required)
    for(int cell_id = 0; cell_id < num_cells; cell_id++)
    {
        // set mult_terms
        mult_terms[0] = C1 * G_field[cell_id] / k_field[cell_id];
        mult_terms[1] = - C2 * epsilon_field[cell_id] / k_field[cell_id];   // note the - sign in the definition;
                                                                            // this is so that the ifs can be generalised

        // loop through mult_terms and add them
        // either to epsilon_source_field or diagonal coeffs
        for(double mult: mult_terms)
        {
            if(mult >= 0)
            {
                epsilon_source_field[cell_id] += mult * epsilon_field[cell_id];
            }
            else
            {
                epsilon_Sp_field[cell_id] += mult;  // no *epsilon_field, it goes to the diagonal terms
            }
        }
    }

    // assemble rhs terms
    calc_collocated_source_coeffs(epsilon_system, 1, epsilon_source_field);
    calc_Sp_term_coeffs(epsilon_system, epsilon_Sp_field);
}

void Turbulence2D::calc_Sp_term_coeffs(LinearSystem& system, std::vector<double>& Sp_field)
{
    // note that Sp has negative values, so we need to add it
    // to diagonal coefficients with a - sign
    for(int cell_id = 0; cell_id < num_cells; cell_id++)
    {
        // the term is multiplied by the volume as it would have been if it stayed on the RHS
        system.add_to_diagonal_coeffs(- Sp_field[cell_id] * cells[cell_id].get_volume(), cell_id);
    }
}