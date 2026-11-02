#include <array>
#include <cassert>
#include <iostream>
#include <vector>

#include "pvCoupling2D.hpp"

// protected

void PVCoupling2D::underrelax_system_implicit(LinearSystem& system, std::vector<double> scalar_field,
                                                double underrel_coeff)
{
    double diag_underrel_coeff = 1./underrel_coeff;
    double rhs_underrel_coeff = (1 - underrel_coeff) / underrel_coeff;

    // get diagonal coefficients before underrelaxation
    std::vector<double> old_diag_coeffs = system.get_diagonal_coeffs();

    for(int cell_id = 0; cell_id < num_cells; cell_id++)
    {
        // scale diagonal coefficient
        system.scale_diagonal_coeff(diag_underrel_coeff, cell_id);

        // add rhs underrelaxation term
        system.add_to_rhs_vector(old_diag_coeffs[cell_id] * rhs_underrel_coeff * scalar_field[cell_id],
                                    cell_id);
    }
}

void PVCoupling2D::underrelax_field_explicit(std::vector<double>& scalar_field, std::vector<double>& scalar_field_old,
                                                double underrel_coeff)
{
    // assertion precaution
    assert(scalar_field.size() == scalar_field_old.size());

    for(int cell_id = 0; cell_id < num_cells; cell_id++)
        {
            scalar_field[cell_id] = scalar_field_old[cell_id]
                                        + underrel_coeff
                                            * (scalar_field[cell_id] - scalar_field_old[cell_id]);
        }
}



void PVCoupling2D::calc_time_der_coeffs(LinearSystem& system, int eq_sign, const std::vector<double>& old_field)
{
    for(int cell_id = 0; cell_id < num_cells; cell_id++)
    {
        // diagonal
        system.add_to_diagonal_coeffs(eq_sign * cells[cell_id].get_volume() / D_T,
                                                cell_id);

        // rhs coefficient
        system.add_to_rhs_vector(eq_sign * cells[cell_id].get_volume() / D_T * old_field[cell_id],
                                                cell_id);
    }
}



void PVCoupling2D::calc_convection_coeffs(LinearSystem& system, int eq_sign,
                                            const std::vector<BoundaryConditions>& field_bcs,
                                            const std::vector<std::vector<double>>& field_bc_values,
                                            InterpolationMethods intpol_method)
{
    switch (intpol_method)
    {
        case InterpolationMethods::UPWIND:
        {
            calc_convection_coeffs_upwind(system, eq_sign, field_bcs, field_bc_values);
            break;
        }
        
        default:
        {
            std::cerr<<"calc_convection_coeffs(): INTERPOLATION METHOD NOT IMPLEMENTED!"<<std::endl;
            std::abort();
            break;
        }
    }
}

void PVCoupling2D::calc_convection_coeffs_upwind(LinearSystem& system, int eq_sign,
                                                    const std::vector<BoundaryConditions>& field_bcs,
                                                    const std::vector<std::vector<double>>& field_bc_values)
{
    // loop through patches to calculate coefficients at boundary faces
    for(int patch_id = 0; patch_id < num_patches; patch_id++)
    {
        std::vector<int> patch_face_ids = patches[patch_id].get_face_ids();
        std::size_t num_patch_faces = patch_face_ids.size();

        switch(field_bcs[patch_id])
        {
            case BoundaryConditions::DIRICHLET:
            {
                std::vector<double> boundary_field = field_bc_values[patch_id];

                for(std::size_t pos_in_patch = 0; pos_in_patch < num_patch_faces; pos_in_patch++)
                {
                    int face_id = patch_face_ids[pos_in_patch];
                    int face_owner_id = faces[face_id].get_owner_id();

                    system.add_to_rhs_vector(- eq_sign * fluxes[face_id]*boundary_field[pos_in_patch],
                                                                face_owner_id);
                }

                break;
            }

            case BoundaryConditions::NEUMANN:
            {
                std::vector<double> flux_unit_area = field_bc_values[patch_id];

                for(std::size_t pos_in_patch = 0; pos_in_patch < num_patch_faces; pos_in_patch++)
                {
                    int face_id = patch_face_ids[pos_in_patch];
                    int face_owner_id = faces[face_id].get_owner_id();
                    Point owner_centroid = cells[face_owner_id].get_centroid();
                    Point face_centroid = faces[face_id].get_centroid();
                    double dist_owner_face = Point::distance(owner_centroid, face_centroid).calc_modulus();

                    double diag_term = eq_sign * fluxes[face_id];

                    // boundary_field = owner_field + dist_owner_face * flux_unit_area;
                    // the modulus of the distance is used as the distance vector is parallel to the face area vector

                    // add diagonal and rhs term
                    // diagonal
                    system.add_to_diagonal_coeffs(diag_term, face_owner_id);

                    // rhs coefficient
                    system.add_to_rhs_vector(- diag_term * dist_owner_face * flux_unit_area[pos_in_patch],
                                                            face_owner_id);
                }

                break;
            }

            default:
            {
                std::cerr<<"calc_convection_coeffs_upwind(): BOUNDARY CONDITION NOT IMPLEMENTED!"<<std::endl;
                std::abort();
                break;
            }
        }
    }

    // loop through faces to calculate the rest of the coefficients
    for(int face_id = 0; face_id < num_faces; face_id++)
    {
        int face_owner_id = faces[face_id].get_owner_id();
        int face_neighbour_id = faces[face_id].get_neighbour_id();

        if(face_neighbour_id == -1) // if at boundary, skip (already added in previous for)
            continue;
        else                        // if not, add diagonal and off-diagonal coefficients
        {
            int flow_sign;

            if(fluxes[face_id] > 0)         // flowing from owner to neighbour
            {
                flow_sign = 1;

                double diag_term = eq_sign * fluxes[face_id] * flow_sign;

                std::array<int, MATRIX_DIMS> matrix_indices {face_neighbour_id,
                                                                face_owner_id};

                // diagonal
                system.add_to_diagonal_coeffs(diag_term, matrix_indices[1]);

                // off-diagonal
                system.add_to_off_diagonal_coeffs(- diag_term, matrix_indices);
            }
            else if(fluxes[face_id] < 0)    // flowing from neighbour to owner
            {
                flow_sign = -1;

                double diag_term = eq_sign * fluxes[face_id] * flow_sign;

                std::array<int, MATRIX_DIMS> matrix_indices {face_owner_id,
                                                                face_neighbour_id};

                // diagonal
                system.add_to_diagonal_coeffs(diag_term, matrix_indices[1]);

                // off-diagonal
                system.add_to_off_diagonal_coeffs(- diag_term, matrix_indices);
            }

            // if flux is 0 do not add any coefficients
        }
    }
}



void PVCoupling2D::calc_diffusion_coeffs(LinearSystem& system, int eq_sign,
                                            const std::vector<BoundaryConditions>& field_bcs,
                                            const std::vector<std::vector<double>>& field_bc_values,
                                            const std::vector<double>& diffusivity_field,
                                            const std::vector<BoundaryConditions>& diffusivity_field_bcs,
                                            const std::vector<std::vector<double>>& diffusivity_field_bc_values,
                                            InterpolationMethods diffusivity_field_intpol_method)
{
    // loop through patches to calculate coefficients at boundary faces
    for(int patch_id = 0; patch_id < num_patches; patch_id++)
    {
        std::vector<int> patch_face_ids = patches[patch_id].get_face_ids();
        std::size_t num_patch_faces = patch_face_ids.size();

        switch(field_bcs[patch_id])
        {
            case BoundaryConditions::DIRICHLET:
            {
                // calculate boundary field
                std::vector<double> boundary_field = field_bc_values[patch_id];

                for(std::size_t pos_in_patch = 0; pos_in_patch < num_patch_faces; pos_in_patch++)
                {
                    // initialisations
                    int face_id = patch_face_ids[pos_in_patch];
                    int face_owner_id = faces[face_id].get_owner_id();
                    Point owner_centroid = cells[face_owner_id].get_centroid();
                    Point face_centroid = faces[face_id].get_centroid();
                    double face_area = faces[face_id].get_area();

                    // calculate distance between owner and face centroid
                    double dist_owner_face = Point::distance(owner_centroid, face_centroid).calc_modulus();

                    double diffusivity_face = interpolate_to_face(diffusivity_field, diffusivity_field_bcs,
                                                                    diffusivity_field_bc_values, face_id,
                                                                    diffusivity_field_intpol_method);

                    // matrix data
                    double diag_term = - diffusivity_face
                                            * face_area / dist_owner_face
                                            * eq_sign;

                    // diagonal
                    system.add_to_diagonal_coeffs(diag_term, face_owner_id);

                    // rhs coefficient
                    system.add_to_rhs_vector(diag_term * boundary_field[pos_in_patch],
                                                            face_owner_id);
                }

                break;
            }

            case BoundaryConditions::NEUMANN:
            {
                std::vector<double> flux_unit_area = field_bc_values[patch_id];

                for(std::size_t pos_in_patch = 0; pos_in_patch < num_patch_faces; pos_in_patch++)
                {
                    // initialisations
                    int face_id = patch_face_ids[pos_in_patch];
                    int face_owner_id = faces[face_id].get_owner_id();
                    double face_area = faces[face_id].get_area();

                    double diffusivity_face = interpolate_to_face(diffusivity_field, diffusivity_field_bcs,
                                                                    diffusivity_field_bc_values, face_id,
                                                                    diffusivity_field_intpol_method);

                    double rhs_term = - diffusivity_face
                                        * face_area * flux_unit_area[pos_in_patch]
                                        * eq_sign;

                    // rhs coefficient
                    system.add_to_rhs_vector(rhs_term, face_owner_id);
                }

                break;
            }

            default:
            {
                std::cerr<<"calc_diffusion_coeffs(): BOUNDARY CONDITION NOT IMPLEMENTED!"<<std::endl;
                std::abort();
                break;
            }
        }
    }

    // loop through faces to calculate the rest of the coefficients
    for(int face_id = 0; face_id < num_faces; face_id++)
    {
        // initialisations
        int face_owner_id = faces[face_id].get_owner_id();
        int face_neighbour_id = faces[face_id].get_neighbour_id();

        if(face_neighbour_id == -1) // if at boundary, skip (already added in previous for)
            continue;
        else
        {
            // more initialisations
            double face_area = faces[face_id].get_area();
            Point owner_centroid = cells[face_owner_id].get_centroid();
            Point neighbour_centroid = cells[face_neighbour_id].get_centroid();

            // calculate distance between owner and neighbour
            double dist_owner_neighbour = Point::distance(owner_centroid, neighbour_centroid).calc_modulus();

            double diffusivity_face = interpolate_to_face(diffusivity_field, diffusivity_field_bcs,
                                                        diffusivity_field_bc_values, face_id, diffusivity_field_intpol_method);

            // matrix data
            double diag_term = - diffusivity_face
                                    * face_area / dist_owner_neighbour
                                    * eq_sign;

            std::array<int, MATRIX_DIMS> matrix_indices {face_owner_id, face_neighbour_id};
            // reverse matrix_indices for symmetric input
            std::array<int, MATRIX_DIMS> reverse_matrix_indices {matrix_indices[1], matrix_indices[0]};

            // diagonal contributions
            system.add_to_diagonal_coeffs(diag_term, matrix_indices[0]);

            system.add_to_diagonal_coeffs(diag_term, matrix_indices[1]);

            // off-diagonal contributions - negative of diagonal
            system.add_to_off_diagonal_coeffs(- diag_term, matrix_indices);

            system.add_to_off_diagonal_coeffs(- diag_term, reverse_matrix_indices);
        }
    }
}



void PVCoupling2D::calc_collocated_source_coeffs(LinearSystem& system, int eq_sign, const std::vector<double>& source_field)
{
    for(int cell_id = 0; cell_id < num_cells; cell_id++)
    {
        // initialisations
        double cell_volume = cells[cell_id].get_volume();

        // add to rhs coeffs
        system.add_to_rhs_vector(eq_sign * cell_volume * source_field[cell_id], cell_id);
    }
}



void PVCoupling2D::assemble_velocity_systems()
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
}

void PVCoupling2D::assemble_pressure_system(const std::vector<double>& diffusivity_field,
                                            const std::vector<BoundaryConditions>& diffusivity_field_bcs,
                                            const std::vector<std::vector<double>>& diffusivity_field_bc_values)
{
    // div(1/apu * grad p) = div(1/apu * H(u))
    calc_diffusion_coeffs(pressure_system, 1,
                            pressure_bcs, pressure_bc_values,
                            diffusivity_field, diffusivity_field_bcs, diffusivity_field_bc_values,
                            diffusivity_intpol_method);   // grad_p is a calculus operation here, not discretised
    calc_flux_source_coeffs(pressure_system, 1, fluxes);
}



void PVCoupling2D::calc_flux_source_coeffs(LinearSystem& system, int eq_sign, const std::vector<double>& flux_field)
{
    // assertion precaution
    assert(flux_field.size() == (size_t) num_faces);

    for(int face_id = 0; face_id < num_faces; face_id++)
    {
        int face_owner_id = faces[face_id].get_owner_id();
        int face_neighbour_id = faces[face_id].get_neighbour_id();

        double rhs_term = eq_sign * flux_field[face_id];

        system.add_to_rhs_vector(rhs_term, face_owner_id);

        if(face_neighbour_id != -1) // if at boundary, add RHS coefficient
        {
            system.add_to_rhs_vector(- rhs_term, face_neighbour_id);
        }
    }
}

void PVCoupling2D::calc_divergence_source_coeffs(LinearSystem& system, int eq_sign, const std::vector<Point>& source_field,
                                                    const std::vector<BoundaryConditions>& source_field_bcs,
                                                    const std::vector<std::vector<Point>>& source_field_bc_values,
                                                    InterpolationMethods& source_field_intpol_method)
{
    // assertion precaution
    assert(source_field.size() == (size_t) num_cells);

    for(int face_id = 0; face_id < num_faces; face_id++)
    {
        int face_owner_id = faces[face_id].get_owner_id();
        int face_neighbour_id = faces[face_id].get_neighbour_id();

        // interpolate to face
        Point face_source = interpolate_to_face(source_field, source_field_bcs, source_field_bc_values,
                                                face_id, source_field_intpol_method);
        double rhs_term = eq_sign * faces[face_id].get_area_vector().dot_product_with(face_source);

        system.add_to_rhs_vector(rhs_term, face_owner_id);

        if(face_neighbour_id != -1) // if at boundary, add RHS coefficient
        {
            system.add_to_rhs_vector(- rhs_term, face_neighbour_id);
        }
    }
}



void PVCoupling2D::assemble_matrix(LinearSystem& system, int num_of_operations, std::function<void ()>& raw_assembly_function,
                                    const std::function<void ()>& off_diagonal_postprocessing)
{
    // reserve appropriate space for the raw off-diagonal vectors
    system.reserve_coeffs_raw(num_of_operations * num_internal_faces);

    // assemble the matrix
    raw_assembly_function();
    
    // populate CompressedRowMatrix-es
    system.populate_off_diagonal_coefficients();

    // do postprocessing on assembled off-diagonal coefficients
    if(off_diagonal_postprocessing)
        off_diagonal_postprocessing();
}

void PVCoupling2D::assemble_matrix(std::array<LinearSystem, N_DIMS>& systems,
                                    int num_of_operations, std::function<void ()>& raw_assembly_function,
                                    const std::function<void ()>& off_diagonal_postprocessing)
{
    // reserve appropriate space for the raw off-diagonal vectors
    for(Dimensions dim: ALL_DIMS)
        systems[dim].reserve_coeffs_raw(num_of_operations * num_internal_faces);

    // assemble the matrix
    raw_assembly_function();
    
    // populate CompressedRowMatrix-es
    for(Dimensions dim: ALL_DIMS)
        systems[dim].populate_off_diagonal_coefficients();

    // do postprocessing on assembled off-diagonal coefficients
    if(off_diagonal_postprocessing)
        off_diagonal_postprocessing();
}