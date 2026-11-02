#include <array>
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "pvCoupling2D.hpp"

// protected

void PVCoupling2D::read_scalar_field(std::ifstream& input_file, std::vector<double>& scalar_field)
{
    // declarations
    std::string line;

    // char vectors for trimming
    std::vector<char> lchars{' ', '\t', '('},
                        rchars{'\r', ')'},
                        EOL_chars{'\r'};    // to accommodate for "\r\n" for new line in Windows

    // read number of points
    std::getline(input_file, line);
    
    assert(num_cells == std::stoi(line));   // the number of field points should be the same as the number of cells

    // check for open bracket to begin read
    std::getline(input_file, line);
    assert(rtrim(line, EOL_chars) == "(");

    for(int i=0; i<num_cells; i++)
    {
        std::getline(input_file, line);

        // trim line from left and right
        line = trim(line, lchars, rchars);

        std::istringstream isstream(line);

        // read field point
        isstream >> scalar_field[i];

        // assert reaching end of line - it was fully read
        assert(isstream.peek() == -1);
    }

    // check for close bracket to end read
    std::getline(input_file, line);
    assert(rtrim(line, EOL_chars) == ")");

    // close input_file
    input_file.close();
}

void PVCoupling2D::read_vector_field(std::ifstream& input_file, std::vector<Point>& vector_field)
{
    // declarations
    std::string line;
    std::array<double, N_DIMS> coords;

    // char vectors for trimming
    std::vector<char> lchars{' ', '\t', '('},
                        rchars{'\r', ')'},
                        EOL_chars{'\r'};    // to accommodate for "\r\n" for new line in Windows

    // read number of points
    std::getline(input_file, line);
    
    assert(num_cells == std::stoi(line));   // the number of field points should be the same as the number of cells

    // check for open bracket to begin read
    std::getline(input_file, line);
    assert(rtrim(line, EOL_chars) == "(");

    for(int i=0; i<num_cells; i++)
    {
        std::getline(input_file, line);

        // trim line from left and right
        line = trim(line, lchars, rchars);

        std::istringstream isstream(line);

        // read vector coords
        for(Dimensions dim: ALL_DIMS)
            isstream >> coords[dim];

        // assert reaching end of line - it was fully read
        assert(isstream.peek() == -1);

        // add vector to field
        Point vector(coords);
        vector_field[i] = vector;
    }

    // check for close bracket to end read
    std::getline(input_file, line);
    assert(rtrim(line, EOL_chars) == ")");

    // close input_file
    input_file.close();
}

void PVCoupling2D::read_boundary_conds(std::ifstream& boundary_conds_file)
{
    // declarations
    std::string line, boundary_cond_type, boundary_patch, variable;
    bool abort_flag = false, patch_flag;

    // read boundary_conds_file
    while(std::getline(boundary_conds_file, line))
    {
        patch_flag = false; // reset patch_flag to "unfound"

        std::istringstream isstream(line);

        // read boundary patch
        isstream >> boundary_patch;

        // loop through patches to check that boundary_patch is valid,
        // perform operations if found
        for(int patch_id = 0; patch_id < num_patches; patch_id++)
        {
            if(patches[patch_id].get_name() == boundary_patch)
            {
                // mark patch_flag as found
                patch_flag = true;

                // get patch data
                int num_patch_faces = patches[patch_id].get_face_ids().size();

                // read boundary condition
                isstream >> boundary_cond_type;

                if(boundary_cond_type == "Dirichlet")       // CONSTANT DIRICHLET
                {
                    // set boundary condition
                    isstream >> variable;

                    if(variable == "P")
                    {
                        // ----------------- PRESSURE -----------------
                        pressure_bcs[patch_id] = BoundaryConditions::DIRICHLET;

                        // extract Dirichlet pressure
                        double dirichlet_pressure;
                        isstream >> dirichlet_pressure;

                        // apply boundary conditions
                        pressure_bc_values[patch_id] = std::vector(num_patch_faces, dirichlet_pressure);
                        
                        // ----------------- VELOCITY -----------------
                        velocity_bcs[patch_id] = BoundaryConditions::NEUMANN;

                        velocity_bc_values[patch_id] = std::vector(num_patch_faces, DVEL_DN_NEUMANN());
                    }
                    else if(variable == "V")
                    {
                        // ----------------- VELOCITY -----------------
                        velocity_bcs[patch_id] = BoundaryConditions::DIRICHLET;

                        // extract Dirichlet velocity
                        std::array<double, N_DIMS> dirichlet_velocity_coords;

                        for(Dimensions dim: ALL_DIMS)
                            isstream >> dirichlet_velocity_coords[dim];

                        Point dirichlet_velocity(dirichlet_velocity_coords);

                        // apply boundary condition
                        velocity_bc_values[patch_id] = std::vector(num_patch_faces, dirichlet_velocity);

                        // ----------------- PRESSURE -----------------
                        pressure_bcs[patch_id] = BoundaryConditions::NEUMANN;

                        pressure_bc_values[patch_id] = std::vector(num_patch_faces, DP_DN_NEUMANN());
                    }
                    else
                    {
                        std::cerr<<"UNRECOGNISED VARIABLE: "<<variable<<"!"<<std::endl;
                        abort_flag = true;
                        break;
                    }

                    // exit patch loop
                    break;
                }

                else if(boundary_cond_type == "Neumann")
                {
                    // set boundary condition
                    isstream >> variable;

                    if(variable == "P")
                    {
                        // ----------------- PRESSURE -----------------
                        pressure_bcs[patch_id] = BoundaryConditions::NEUMANN;

                        // extract Neumann flux
                        double neumann_flux;
                        isstream >> neumann_flux;

                        // apply boundary conditions
                        pressure_bc_values[patch_id] = std::vector(num_patch_faces, neumann_flux);
                        
                        // ----------------- VELOCITY -----------------
                        velocity_bcs[patch_id] = BoundaryConditions::DIRICHLET;

                        velocity_bc_values[patch_id] = std::vector(num_patch_faces, VEL_DIRICHLET());
                    }
                    else if(variable == "V")
                    {
                        // ----------------- VELOCITY -----------------
                        velocity_bcs[patch_id] = BoundaryConditions::NEUMANN;

                        // extract Neumann flux
                        std::array<double, N_DIMS> neumann_flux_coords;

                        for(Dimensions dim: ALL_DIMS)
                            isstream >> neumann_flux_coords[dim];

                        Point neumann_flux(neumann_flux_coords);

                        // apply boundary condition
                        velocity_bc_values[patch_id] = std::vector(num_patch_faces, neumann_flux);

                        // ----------------- PRESSURE -----------------
                        pressure_bcs[patch_id] = BoundaryConditions::DIRICHLET;

                        pressure_bc_values[patch_id] = std::vector(num_patch_faces, P_DIRICHLET());
                    }
                    else
                    {
                        std::cerr<<"UNRECOGNISED VARIABLE: "<<variable<<"!"<<std::endl;
                        abort_flag = true;
                        break;
                    }

                    // exit patch loop
                    break;
                }

                else
                {
                    std::cerr<<"UNRECOGNISED BOUNDARY CONDITION: "<<boundary_cond_type<<"!"<<std::endl;
                    abort_flag = true;
                }
            }
        }

        if(!patch_flag) // if not found
        {
            std::cerr<<"UNRECOGNISED BOUNDARY PATCH: "<<boundary_patch<<"!"<<std::endl;
            abort_flag = true;
        }

        // assert that end of line is reached
        assert(isstream.peek() == -1);

        // WARNING: MISSING OTHER BOUNDARY CONDITIONS, E.G. MIXED
    }

    // check that all boundary conditions have been assigned
    for(int patch_id = 0; patch_id < num_patches; patch_id++)
    {
        if(pressure_bcs[patch_id] == BoundaryConditions::NONE)
        {
            std::cerr<<"ERROR: PRESSURE BOUNDARY CONDITION NOT ASSIGNED FOR "
                    <<patches[patch_id].get_name()<<std::endl;
            abort_flag = true;
        }
        if(velocity_bcs[patch_id] == BoundaryConditions::NONE)
        {
            std::cerr<<"ERROR: VELOCITY BOUNDARY CONDITION NOT ASSIGNED FOR "
                    <<patches[patch_id].get_name()<<std::endl;
            abort_flag = true;
        }
    }

    // abort if flagged
    if(abort_flag)
        std::abort();
}



void PVCoupling2D::write_output(std::ofstream& output_file)
{
    for(int cell_id = 0; cell_id < num_cells; cell_id++)
    {
        // initialisations
        std::array<double, N_DIMS> centroid_coords = cells[cell_id].get_centroid().get_coordinates();
        std::array<double, N_DIMS> velocity_coords = velocity_field[cell_id].get_coordinates();
        double velocity_modulus = velocity_field[cell_id].calc_modulus();

        // write cell centroid
        for(Dimensions dim: ALL_DIMS)
            output_file<<centroid_coords[dim]<<"\t";

        // write pressure
        output_file<<pressure_field[cell_id]<<"\t";

        // write velocity vector as (centroid-velocity/2, centroid+velocity/2);
        // scale reciprocally wrt modulus
        double grid_scaling = 0.1/15;
        double velocity_scaling = (1 - 1./(2*velocity_modulus + 1)) * grid_scaling;

        // beginning coordinates
        for(Dimensions dim: ALL_DIMS)
            output_file<<centroid_coords[dim] - velocity_coords[dim]/2 / (velocity_modulus+EPSILON) * velocity_scaling<<"\t";

        // ending coordinates
        for(Dimensions dim: ALL_DIMS)
            output_file<<velocity_coords[dim] / (velocity_modulus+EPSILON) * velocity_scaling<<"\t";

        // velocity magnitude
            output_file<<velocity_field[cell_id].calc_modulus()<<"\t";

        output_file<<std::endl;
    }

    output_file<<std::endl<<std::endl;
}

void PVCoupling2D::write_linear_systems_data(std::ofstream& output_file_matrices, std::ofstream& output_file_vectors, std::vector<LinearSystem>& linear_systems,
                                                std::vector<std::vector<double>> residuals_before, std::vector<std::vector<double>> residuals_after)
{
    size_t num_systems = linear_systems.size();

    // assertion precaution
    assert(residuals_before.size() == num_systems);
    assert(residuals_after.size() == num_systems);

    // declarations
    std::vector<std::vector<double>> diagonal_coefficients(num_systems);
    std::vector<std::vector<double>> rhs_vectors(num_systems);
    std::vector<MatrixElmt*> upper_elmt_iters(num_systems);
    std::vector<MatrixElmt*> lower_elmt_iters(num_systems);

    // read linear_systems
    for(size_t i = 0; i < num_systems; i++)
    {
        diagonal_coefficients[i] = linear_systems[i].get_diagonal_coeffs();
        rhs_vectors[i] = linear_systems[i].get_rhs_vector();

        upper_elmt_iters[i] = linear_systems[i].get_upper_coeffs().get_first_elmt();

        lower_elmt_iters[i] = linear_systems[i].get_lower_coeffs().get_first_elmt();
    }

    // ------------------- WRITE MATRIX FILE -------------------
    for(int cell_i = 0; cell_i < num_cells; cell_i++)
    {
        for(int cell_j = 0; cell_j < num_cells; cell_j++)
        {
            // output matrix index
            output_file_matrices<<cell_i<<"\t"<<cell_j<<"\t";

            // output [cell_i, cell_j] matrix coefficient for i-th matrix
            for(size_t i = 0; i < num_systems; i++)
            {
                if(cell_i == cell_j)
                {
                    output_file_matrices<<diagonal_coefficients[i][cell_i]<<"\t";
                }

                else if(cell_i < cell_j)    // upper coefficients
                {
                    bool alloc_flag = false;
                    std::array<int, MATRIX_DIMS> curr_indices {cell_i, cell_j};

                    while(upper_elmt_iters[i])  // loop through iters until curr_indices are reached
                    {
                        std::array<int, MATRIX_DIMS> iter_indices = upper_elmt_iters[i]->indices;

                        for(int k = 0; k < MATRIX_DIMS; k++)
                        {
                            if(curr_indices[k] < iter_indices[k])   // if current indices are before the iter indices
                            {
                                output_file_matrices<<"0"<<"\t";    // curr_indices have zero value and are not in upper_coeffs[i]
                                alloc_flag = true;
                                break;
                            }
                            else if(curr_indices[k] > iter_indices[k])  // if current indices are after the iter indices
                            {
                                break;                                  // break and go to next upper_elmt_iters
                            }
                            else if(k == MATRIX_DIMS-1) // if we reached the final iteration
                            {
                                output_file_matrices<<upper_elmt_iters[i]->value<<"\t"; // then curr_indices == iter_indices and we print the value
                                alloc_flag = true;
                                break;
                            }
                        }

                        if(alloc_flag)
                            break;
                        else
                            upper_elmt_iters[i] = upper_elmt_iters[i]->next_elmt;
                    }

                    if(!alloc_flag)
                        output_file_matrices<<"0"<<"\t";
                }

                else    // lower coefficients
                {
                    bool alloc_flag = false;
                    std::array<int, MATRIX_DIMS> curr_indices {cell_i, cell_j};

                    while(lower_elmt_iters[i])  // loop through iters until curr_indices are reached
                    {
                        std::array<int, MATRIX_DIMS> iter_indices = lower_elmt_iters[i]->indices;

                        for(int k = 0; k < MATRIX_DIMS; k++)
                        {
                            if(curr_indices[k] < iter_indices[k])   // if current indices are before the iter indices
                            {
                                output_file_matrices<<"0"<<"\t";    // curr_indices have zero value and are not in lower_coeffs[i]
                                alloc_flag = true;
                                break;
                            }
                            else if(curr_indices[k] > iter_indices[k])  // if current indices are after the iter indices
                            {
                                break;                                  // break and go to next lower_elmt_iters
                            }
                            else if(k == MATRIX_DIMS-1) // if we reached the final iteration
                            {
                                output_file_matrices<<lower_elmt_iters[i]->value<<"\t"; // then curr_indices == iter_indices and we print the value
                                alloc_flag = true;
                                break;
                            }
                        }

                        if(alloc_flag)
                            break;
                        else
                            lower_elmt_iters[i] = lower_elmt_iters[i]->next_elmt;
                    }

                    if(!alloc_flag)
                        output_file_matrices<<"0"<<"\t";
                }
            }

            output_file_matrices<<std::endl;
        }

        output_file_matrices<<std::endl;
    }

    // new line for separation
    output_file_matrices<<std::endl;

    // ------------------- WRITE VECTOR FILE -------------------
    for(int cell_id = 0; cell_id < num_cells; cell_id++)
    {
        // initialisations
        std::array<double, N_DIMS> centroid_coords = cells[cell_id].get_centroid().get_coordinates();

        // write cell centroid
        for(Dimensions dim: ALL_DIMS)
            output_file_vectors<<centroid_coords[dim]<<"\t";

        // for each linear system, write its rhs and residual vector
        for(size_t i = 0; i < num_systems; i++)
        {
            output_file_vectors<<diagonal_coefficients[i][cell_id]<<"\t";
            output_file_vectors<<rhs_vectors[i][cell_id]<<"\t";
            output_file_vectors<<residuals_before[i][cell_id]<<"\t";
            output_file_vectors<<residuals_after[i][cell_id]<<"\t";
        }

        output_file_vectors<<std::endl;
    }

    // new line for separation
    output_file_vectors<<std::endl<<std::endl;
}