#include <cassert>
#include <cmath>
#include <iostream>

#include "pvCoupling2D.hpp"

// protected

void PVCoupling2D::calc_char_length()
{
    for(int patch_id = 0; patch_id < num_patches; patch_id++)
    {
        if(patches[patch_id].get_name().find("Inlet") != std::string::npos)
        {
            char_length = 0;

            // calculate char_length in 2D as a sum of all the 2D face areas ("face lengths")
            for(int face_id: patches[patch_id].get_face_ids())
                char_length += faces[face_id].get_area();

            break;
        }
    }
}

void PVCoupling2D::calc_nu_field()
{
    // get velocity modulus at the inlet, otherwise use 1
    double avg_velocity_modulus_inlet = 1;

    for(int patch_id = 0; patch_id < num_patches; patch_id++)
    {
        if(patches[patch_id].get_name().find("Inlet") != std::string::npos)
        {
            // (re) initialisations
            std::vector<int> patch_face_ids = patches[patch_id].get_face_ids();
            std::size_t num_patch_faces = patch_face_ids.size();

            avg_velocity_modulus_inlet = 0;

            // sum all boundary velocities
            for(std::size_t pos_in_patch = 0; pos_in_patch < num_patch_faces; pos_in_patch++)
                avg_velocity_modulus_inlet += velocity_bc_values[patch_id][pos_in_patch].calc_modulus();

            avg_velocity_modulus_inlet /= num_patch_faces;

            break;
        }
    }

    // calculate nu and populate nu_field
    nu = char_length * avg_velocity_modulus_inlet / Re;

    nu_field = std::vector<double>(num_cells, nu);
}

void PVCoupling2D::calc_fluxes(InterpolationMethods intpol_method)
{
    for(int face_id = 0; face_id < num_faces; face_id++)
    {
        // interpolate velocity field to face_id
        Point face_velocity = interpolate_to_face(velocity_field, velocity_bcs, velocity_bc_values, face_id, intpol_method);

        // compute the face flux;
        // sign is positive if pointing from owner to neighbour
        fluxes[face_id] = faces[face_id].get_area_vector().dot_product_with(face_velocity);
    }
}

void PVCoupling2D::check_fluxes(bool check_cells)
{
    // initialisations
    double total_boundary_flux = 0;

    // check fluxes on all boundaries sum to zero
    for(int patch_id = 0; patch_id < num_patches; patch_id++)
    {
        for(int face_id: patches[patch_id].get_face_ids())
            total_boundary_flux += fluxes[face_id];
    }

    if(fabs(total_boundary_flux) > FLUX_TOL)
    {
        std::cerr<<"WARNING! BOUNDARY HAS NON-ZERO TOTAL FLUX: "<<total_boundary_flux<<std::endl;
        std::cerr<<"TOLERANCE: "<<FLUX_TOL<<std::endl;
    }

    // check all fluxes sum to zero, as well as in each cell
    if(check_cells)
    {
        // initialisations
        double total_mesh_flux = 0;

        for(int cell_id = 0; cell_id < num_cells; cell_id++)
        {
            double total_cell_flux = 0;

            for(int face_id: cells[cell_id].get_face_ids())
                total_cell_flux += fluxes[face_id]
                                    * ((cell_id == faces[face_id].get_owner_id())? 1 : -1);

            total_mesh_flux += total_cell_flux;
        }

        if(fabs(total_mesh_flux) > FLUX_TOL)
        {
            std::cerr<<"WARNING! MESH HAS NON-ZERO TOTAL FLUX: "<<total_mesh_flux<<std::endl;
            std::cerr<<"TOLERANCE: "<<FLUX_TOL<<std::endl;
        }
    }
}

void PVCoupling2D::update_fluxes_with_pressure(const std::vector<double>& diffusivity_field,
                                                const std::vector<BoundaryConditions>& diffusivity_field_bcs,
                                                const std::vector<std::vector<double>>& diffusivity_field_bc_values)
{
    for(int face_id = 0; face_id < num_faces; face_id++)
        {
            Point face_grad_p = calc_face_grad_p(face_id);
            double update = faces[face_id].get_area_vector().dot_product_with(face_grad_p);  // sf * (grad p)f

            if(faces[face_id].get_neighbour_id() == -1)
                update *= diffusivity_field[faces[face_id].get_owner_id()]; // diffusivity_field_p = 1/ap by definition
            else
                update *= interpolate_to_face(diffusivity_field, diffusivity_field_bcs, diffusivity_field_bc_values,
                                                face_id, InterpolationMethods::LINEAR);
            
            fluxes[face_id] -= update;
        }
}

void PVCoupling2D::update_velocity_field_with_grad_p(std::vector<double>& diffusivity_field)
{
    // correct cell centre velocities: u_corr = u_pre - 1/ap grad p
    for(int cell_id = 0; cell_id < num_cells; cell_id++)
    {
        // move velocity field by delta
        Point delta = grad_p_field[cell_id];
        delta.scale_by(-diffusivity_field[cell_id]);    // diffusivity_field_p = 1/ap by definition
        velocity_field[cell_id].move_by(delta);
    }
}



void PVCoupling2D::extract_dim_field(std::vector<Point>& field, std::vector<double>& dim_field, Dimensions dim)
{
    // initialisations
    size_t field_size = field.size();

    // resize dim_field before assignment
    dim_field.resize(field_size);

    // extract dim_field
    for(size_t i = 0; i < field_size; i++)
        dim_field[i] = field[i].get_coordinates()[dim];
}

void PVCoupling2D::extract_dim_field(std::vector<std::vector<Point>>& field, std::vector<std::vector<double>>& dim_field, Dimensions dim)
{
    // initialisations
    size_t field_size = field.size();

    // resize dim_field before assignment
    dim_field.resize(field_size);

    // extract dim_field
    for(size_t i = 0; i < field_size; i++)
        extract_dim_field(field[i], dim_field[i], dim);
}

void PVCoupling2D::calc_pseudovelocity_field(const std::vector<double>& diffusivity_field)
{
    // declare H(u) operator, where H(u)=r(u)-LU(u)
    std::vector<std::array<double, N_DIMS>> H_velocity_coords;
    H_velocity_coords.resize(num_cells);

    // calculate H(u) operator
    for(Dimensions dim: ALL_DIMS)
    {
        // get velocity coordinates
        std::vector<double> velocity_dim_field;
        extract_dim_field(velocity_field, velocity_dim_field, dim);

        std::vector<double> H_velocity_dim = velocity_systems[dim].calc_H_operator(velocity_dim_field);

        for(int cell_id = 0; cell_id < num_cells; cell_id++)
            H_velocity_coords[cell_id][dim] += H_velocity_dim[cell_id];
    }

    // initialise velocity_field to H(u)
    for(int cell_id = 0; cell_id < num_cells; cell_id++)
        velocity_field[cell_id] = Point(H_velocity_coords[cell_id]);

    // scale velocity_field by diffusivity field to obtain
    // u_pseudo = 1/a_p * H(u)
    for(int cell_id = 0; cell_id < num_cells; cell_id++)
        velocity_field[cell_id].scale_by(diffusivity_field[cell_id]);
}



void PVCoupling2D::calc_kinematic_pressure_field()
{
    // assumption, the pressure field currently contains dynamic pressure

    // convert dynamic to kinematic pressure
    for(int cell_id = 0; cell_id < num_cells; cell_id++)
    {
        pressure_field[cell_id] /= rho;     // kinematic_p = dynamic_p / rho
    }
}

std::vector<double> PVCoupling2D::calc_pressure_diffusivity_field()
{
    // assess that diagonal coefficients are the same
    std::array<std::vector<double>, N_DIMS> diag_coeffs;
    for(Dimensions dim: ALL_DIMS)
        diag_coeffs[dim] = velocity_systems[dim].get_diagonal_coeffs();

    assert(diag_coeffs[X] == diag_coeffs[Y]);

    // assemble diffusivity field for p
    std::vector<double> diffusivity_field_p = velocity_systems[X].get_diagonal_coeffs();    // above check confirms
                                                                                            // they are the same
    int diffusivity_field_p_size = diffusivity_field_p.size();
    for(int i=0; i<diffusivity_field_p_size; i++)
        diffusivity_field_p[i] = 1./diffusivity_field_p[i];

    return diffusivity_field_p;
}

Point PVCoupling2D::calc_face_grad_p(int face_id)
{
    // initialisations
    int face_owner_id = faces[face_id].get_owner_id();
    Point owner_centroid = cells[face_owner_id].get_centroid();

    int face_neighbour_id = faces[face_id].get_neighbour_id();

    // calculate unit area vector
    Point unit_area_vector = faces[face_id].get_area_vector();
    unit_area_vector.scale_by(1./unit_area_vector.calc_modulus());

    // initialise face_grad_p as the unit area vector
    Point face_grad_p = unit_area_vector;

    if(face_neighbour_id == -1)     // if at boundary
    {
        int patch_id, pos_in_patch;
        calc_patch_and_position(face_id, patch_id, pos_in_patch);

        // apply boundary conditions based on patch
        switch (pressure_bcs[patch_id])
        {
            case BoundaryConditions::DIRICHLET:
            {
                Point face_centroid = faces[face_id].get_centroid();

                double boundary_pressure = pressure_bc_values[patch_id][pos_in_patch];

                // scale face_grad_p by the gradient modulus, linear assumption
                face_grad_p.scale_by( (boundary_pressure - pressure_field[face_owner_id])
                                    / Point::distance(owner_centroid, face_centroid)
                                        .calc_modulus());

                break;
            }

            case BoundaryConditions::NEUMANN:
            {
                double flux_unit_area = pressure_bc_values[patch_id][pos_in_patch];

                // scale face_grad_p according to its boundary condition, which is n \cdot face_grad_p
                face_grad_p.scale_by(flux_unit_area);

                break;
            }
            
            default:
            {
                std::cerr<<"calc_pressure_gradient(): BOUNDARY CONDITION NOT IMPLEMENTED!"<<std::endl;
                std::abort();
                break;
            }
        }
    }
    else
    {
        // more initialisations
        Point neighbour_centroid = cells[face_neighbour_id].get_centroid();

        // scale face_grad_p by the gradient modulus, linear assumption
        face_grad_p.scale_by( (pressure_field[face_neighbour_id] - pressure_field[face_owner_id])
                            / Point::distance(owner_centroid, neighbour_centroid)
                                .calc_modulus());
    }

    return face_grad_p;
}

void PVCoupling2D::calc_grad_field_of_scalar(const std::vector<double>& scalar_field, std::vector<Point>& grad_field,
                                                const std::vector<BoundaryConditions>& field_bcs,
                                                const std::vector<std::vector<double>>& field_bc_vals, InterpolationMethods field_intpol_method)
{
    // assertion precaution
    assert(scalar_field.size() == (size_t) num_cells);
    assert(grad_field.size() == (size_t) num_cells);
    assert(field_bcs.size() == (size_t) num_patches);
    assert(field_bc_vals.size() == (size_t) num_patches);

    // empty grad_field
    grad_field = std::vector<Point> (num_cells);

    for(int face_id = 0; face_id < num_faces; face_id++)
    {
        // initialisations
        int face_owner_id = faces[face_id].get_owner_id();
        int face_neighbour_id = faces[face_id].get_neighbour_id();
        Point face_area_vector = faces[face_id].get_area_vector();

        // update owner cell
        // initialise sf * scalar f_term as face area vector
        Point owner_update = face_area_vector;
        double scalar_at_face = interpolate_to_face(scalar_field, field_bcs, field_bc_vals, face_id, field_intpol_method);
        owner_update.scale_by(scalar_at_face);

        grad_field[face_owner_id].move_by(owner_update);

        // update neighbour cell if possible
        if(face_neighbour_id != -1)
        {
            Point neighbour_update = owner_update;
            neighbour_update.negate(); // negation to orient it away from neighbour centroid

            grad_field[face_neighbour_id].move_by(neighbour_update);
        }
    }

    // scale down all grad field by cell volume
    for(int cell_id = 0; cell_id < num_cells; cell_id++)
    {
        grad_field[cell_id].scale_by(1./cells[cell_id].get_volume());
    }
}



void PVCoupling2D::update_velocity_field_vars(std::function<void ()>& velocity_update_func)
{
    // update velocity
    velocity_update_func();

    // WARNING - calc_face_velocity_field not implemented!
}

void PVCoupling2D::update_pressure_field_vars(std::function<void ()>& pressure_update_func)
{
    // update pressure
    pressure_update_func();

    // calculate grad p at cells
    calc_grad_field_of_scalar(pressure_field, grad_p_field, pressure_bcs, pressure_bc_values, pressure_intpol_method);
}



double PVCoupling2D::interpolate_to_face(const std::vector<double>& field, const std::vector<BoundaryConditions>& field_bcs,
                                            const std::vector<std::vector<double>>& field_bc_values,
                                            int face_id, InterpolationMethods intpol_method)
{
    // initialisations
    int face_owner_id = faces[face_id].get_owner_id();
    int face_neighbour_id = faces[face_id].get_neighbour_id();
    double result;      // initialised to 0

    if(face_neighbour_id == -1)     // if at boundary
    {
        int patch_id, pos_in_patch;
        calc_patch_and_position(face_id, patch_id, pos_in_patch);

        // apply boundary conditions based on patch
        switch(field_bcs[patch_id])
        {
            case BoundaryConditions::DIRICHLET:
            {
                result = field_bc_values[patch_id][pos_in_patch];     // already at face

                break;
            }

            case BoundaryConditions::NEUMANN:
            {
                // get owner and face data
                Point owner_centroid = cells[face_owner_id].get_centroid();
                Point face_centroid = faces[face_id].get_centroid();
                double dist_owner_face = Point::distance(owner_centroid, face_centroid).calc_modulus();

                switch(intpol_method)
                {
                    case InterpolationMethods::LINEAR:
                    {
                        result = field[face_owner_id] + dist_owner_face * field_bc_values[patch_id][pos_in_patch];

                        break;
                    }

                    default:
                    {
                        std::cerr<<"interpolate_to_face(double): INTERPOLATION METHOD NOT IMPLEMENTED!"<<std::endl;
                        std::abort();
                        break;
                    }
                }

                break;
            }

            default:
            {
                std::cerr<<"interpolate_to_face(double): BOUNDARY CONDITION NOT IMPLEMENTED!"<<std::endl;
                std::abort();
                break;
            }
        }
    }
    else
    {
        // more initialisations
        Point owner_centroid = cells[face_owner_id].get_centroid();
        Point neighbour_centroid = cells[face_neighbour_id].get_centroid();

        double owner_value = field[face_owner_id];
        double neighbour_value = field[face_neighbour_id];

        switch (intpol_method)
        {
            case InterpolationMethods::LINEAR:
            {
                // more initialisations
                double dist_owner_neighbour = Point::distance(owner_centroid, neighbour_centroid).calc_modulus();
                Point face_centroid = faces[face_id].get_centroid();

                // owner_ratio = |fN| / |ON|
                double dist_neighbour_face = Point::distance(neighbour_centroid, face_centroid).calc_modulus();
                double owner_ratio = dist_neighbour_face / dist_owner_neighbour;

                result = owner_ratio * owner_value + (1-owner_ratio) * neighbour_value;

                break;
            }
                    
            default:
            {
                std::cerr<<"interpolate_to_face(double): INTERPOLATION METHOD NOT IMPLEMENTED!"<<std::endl;
                std::abort();
                break;
            }
        }
    }

    return result;
}

Point PVCoupling2D::interpolate_to_face(const std::vector<Point>& field, const std::vector<BoundaryConditions>& field_bcs,
                                            const std::vector<std::vector<Point>>& field_bc_values,
                                            int face_id, InterpolationMethods intpol_method)
{
    // initialisations
    int face_owner_id = faces[face_id].get_owner_id();
    int face_neighbour_id = faces[face_id].get_neighbour_id();
    Point result;       // initialised to origin

    if(face_neighbour_id == -1)     // if at boundary
    {
        int patch_id, pos_in_patch;
        calc_patch_and_position(face_id, patch_id, pos_in_patch);

        // apply boundary conditions based on patch
        switch(field_bcs[patch_id])
        {
            case BoundaryConditions::DIRICHLET:
            {
                result = field_bc_values[patch_id][pos_in_patch];     // already at face

                break;
            }

            case BoundaryConditions::NEUMANN:
            {
                // get owner and face data
                Point owner_centroid = cells[face_owner_id].get_centroid();
                Point face_centroid = faces[face_id].get_centroid();
                double dist_owner_face = Point::distance(owner_centroid, face_centroid).calc_modulus();

                Point owner_field = field[face_owner_id];

                switch(intpol_method)
                {
                    case InterpolationMethods::LINEAR:
                    {
                        // result = owner_field + dist_owner_face * (n \cdot grad(field))_face_id
                        result = field_bc_values[patch_id][pos_in_patch];
                        result.scale_by(dist_owner_face);
                        result.move_by(owner_field);

                        break;
                    }

                    default:
                    {
                        std::cerr<<"interpolate_to_face(Point): INTERPOLATION METHOD NOT IMPLEMENTED!"<<std::endl;
                        std::abort();
                        break;
                    }
                }

                break;
            }

            default:
            {
                std::cerr<<"interpolate_to_face(Point): BOUNDARY CONDITION NOT IMPLEMENTED!"<<std::endl;
                std::abort();
                break;
            }
        }
    }
    else
    {
        // more initialisations
        Point owner_centroid = cells[face_owner_id].get_centroid();
        Point neighbour_centroid = cells[face_neighbour_id].get_centroid();
        double dist_owner_neighbour = Point::distance(owner_centroid, neighbour_centroid).calc_modulus();
        Point face_centroid = faces[face_id].get_centroid();

        Point owner_value = field[face_owner_id];
        Point neighbour_value = field[face_neighbour_id];

        switch (intpol_method)
        {
            case InterpolationMethods::LINEAR:
            {
                // owner_ratio = |fN| / |ON|
                double dist_neighbour_face = Point::distance(neighbour_centroid, face_centroid).calc_modulus();
                double owner_ratio = dist_neighbour_face / dist_owner_neighbour;

                owner_value.scale_by(owner_ratio);
                neighbour_value.scale_by(1-owner_ratio);

                result.move_by(owner_value);
                result.move_by(neighbour_value);

                break;
            }
                    
            default:
            {
                std::cerr<<"interpolate_to_face(Point): INTERPOLATION METHOD NOT IMPLEMENTED!"<<std::endl;
                std::abort();
                break;
            }
        }
    }

    return result;
}



void PVCoupling2D::calc_patch_and_position(int face_id, int& patch_id, int& position_in_patch)
{
    // loop through patches to find which one face_id belongs to
    for(int i_patch = 0; i_patch < num_patches; i_patch++)
    {
        std::vector<int> patch_face_ids = patches[i_patch].get_face_ids();
        int num_patch_faces = patch_face_ids.size();

        for(int i_pos = 0; i_pos < num_patch_faces; i_pos++)
        {
            if(face_id == patch_face_ids[i_pos])
            {
                patch_id = i_patch;
                position_in_patch = i_pos;

                return;
            }
        }
    }

    // set numbers to -1 if patch not found
    patch_id = -1;
    position_in_patch = -1;
}