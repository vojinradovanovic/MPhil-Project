#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "pvCoupling2D.hpp"

// public

PVCoupling2D::PVCoupling2D(): FaceMesh2D()
{}

void PVCoupling2D::initialise_params()
{
    // initialise mesh
    initialise();

    // initialise non-scalar variables
    pressure_system = LinearSystem(num_cells);
    pressure_field.resize(num_cells);
    grad_p_field.resize(num_cells);
    pressure_bcs = std::vector<BoundaryConditions> (num_patches, BoundaryConditions::NONE);
    pressure_bc_values.resize(num_patches);
    diffusivity_field_p_bcs = std::vector<BoundaryConditions> (num_patches, BoundaryConditions::NEUMANN);
    diffusivity_field_p_bc_values.resize(num_patches);

    for(Dimensions dim: ALL_DIMS)
        velocity_systems[dim] = LinearSystem(num_cells);

    velocity_field.resize(num_cells);
    velocity_bcs = std::vector<BoundaryConditions> (num_patches, BoundaryConditions::NONE);
    velocity_bc_values.resize(num_patches);

    fluxes.resize(num_faces);

    // initialise constants
    Re = 100;
    rho = 1;
    char_length = 0.1;
    calc_char_length();     // calculate char length if inlet is present
    nu = 0;

    velocity_underrel_coeff = 0.8;
    pressure_underrel_coeff = 0.2;

    // initialise nu_field
    nu_field = std::vector<double>(num_cells, nu);
    nu_bcs = std::vector<BoundaryConditions> (num_patches, BoundaryConditions::NEUMANN);
    nu_bc_values.resize(num_patches);

    // initialise default interpolation methods
    convection_intpol_method = InterpolationMethods::UPWIND;
    pressure_intpol_method = InterpolationMethods::LINEAR;
    velocity_intpol_method = InterpolationMethods::LINEAR;
    diffusivity_intpol_method = InterpolationMethods::LINEAR;

    // other
    num_internal_faces = get_num_faces();
    for(BoundaryPatch p: patches)
    {
        num_internal_faces -= p.get_face_ids().size();
    }

    // resize all bc_values vectors to respective number of faces before writing into them
    for(int patch_id = 0; patch_id < num_patches; patch_id++)
    {
        std::size_t num_patch_faces = patches[patch_id].get_face_ids().size();

        pressure_bc_values[patch_id].resize(num_patch_faces);
        diffusivity_field_p_bc_values[patch_id].resize(num_patch_faces);
        velocity_bc_values[patch_id].resize(num_patch_faces);
        nu_bc_values[patch_id].resize(num_patch_faces);
    }

    std::cerr<<"---------------------------------------------"<<std::endl;
    std::cerr<<"WARNINGS:"<<std::endl;
    std::cerr<<std::endl;
    std::cerr<<"DIFFUSIVITY FIELD nu_field IS INITIALISED AS CONSTANT ("<<nu<<") THROUGHOUT THE ENTIRE MESH!"<<std::endl;
    std::cerr<<"VELOCITY FOR INITIALISING nu_field IS: (1, 0)!"<<std::endl;
    std::cerr<<std::endl;
    std::cerr<<"DIFFUSION COEFFICIENTS COMPUTED ASSUMING ORTHOGONAL MESH!"<<std::endl;
    std::cerr<<"DISTANCE FROM CELL TO FACE COMPUTED ASSUMING ORTHOGONAL MESH!"<<std::endl;
    std::cerr<<std::endl;
    std::cerr<<"VELOCITY DIRICHLET CONDITION ALWAYS IMPLIES PRESSURE 0-NEUMANN CONDITION!"<<std::endl;
    std::cerr<<"VELOCITY NEUMANN CONDITION ALWAYS IMPLIES PRESSURE 0-DIRICHLET CONDITION!"<<std::endl;
    std::cerr<<"PRESSURE DIRICHLET CONDITION ALWAYS IMPLIES VELOCITY 0-NEUMANN CONDITION!"<<std::endl;
    std::cerr<<"PRESSURE NEUMANN CONDITION ALWAYS IMPLIES VELOCITY 1-DIRICHLET CONDITION!"<<std::endl;
    std::cerr<<std::endl;
    std::cerr<<"calc_face_grad_p IS COMPUTED FROM PRESSURE VALUES, NOT grad_p_field VALUES!"<<std::endl;
    std::cerr<<"---------------------------------------------"<<std::endl;
}



void PVCoupling2D::set_conds()
{
    //  ------------------- READ BOUNDARY CONDITIONS -------------------
    // declarations
    std::string boundary_conds_path;
    std::ifstream boundary_conds_file;

    std::cout<<"Write file path for boundary conditions:"<<std::endl;
    std::cin>>boundary_conds_path;
    boundary_conds_file.open(boundary_conds_path);

    read_boundary_conds(boundary_conds_file);

    //  ------------------- READ VARIABLES AND COMPUTE DEPENDENT FIELDS -------------------
    // declarations
    std::string dynamic_pressure_path, velocity_path;
    std::ifstream dynamic_pressure_file, velocity_file;

    // open pressure data file
    std::cout<<"Write file path for dynamic pressure:"<<std::endl;
    std::cin>>dynamic_pressure_path;
    dynamic_pressure_file.open(dynamic_pressure_path);

    // update pressure field and dependents with pressure_file
    std::function<void ()> pressure_update_func =
        [this, &dynamic_pressure_file]() {
                                    // read dynamic pressure
                                    read_scalar_field(dynamic_pressure_file, pressure_field);

                                    // update dynamic to kinematic pressure
                                    calc_kinematic_pressure_field();
                                };

    update_pressure_field_vars(pressure_update_func);

    // read velocity data
    std::cout<<"Write file path for velocity:"<<std::endl;
    std::cin>>velocity_path;
    velocity_file.open(velocity_path);

    // update velocity field and dependents with velocity_file
    std::function<void ()> velocity_update_func =
        [this, &velocity_file]() { read_vector_field(velocity_file, velocity_field); };

    update_velocity_field_vars(velocity_update_func);

    //  ------------------- COMPUTE OTHER FIELDS -------------------
    // calculate molecular viscosity based on inlet velocity
    calc_nu_field();

    // assemble fluxes: compute s_f cdot u_old_f
    calc_fluxes(velocity_intpol_method);
}