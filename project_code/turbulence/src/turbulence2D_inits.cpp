#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "tensor.hpp"
#include "turbulence2D.hpp"

// public

Turbulence2D::Turbulence2D(): PVCoupling2D()
{}

void Turbulence2D::initialise_params()
{
    // initialise params from parent class
    PVCoupling2D::initialise_params();

    // initialise non-scalar variables
    k_system = LinearSystem(num_cells);
    k_field.resize(num_cells, EPSILON);          // initialise field to lowest positive value
    k_bcs = std::vector<BoundaryConditions> (num_patches, BoundaryConditions::NONE);
    k_bc_values.resize(num_patches);

    epsilon_system = LinearSystem(num_cells);
    epsilon_field.resize(num_cells, EPSILON);    // initialise field to lowest positive value
    epsilon_bcs = std::vector<BoundaryConditions> (num_patches, BoundaryConditions::NONE);
    epsilon_bc_values.resize(num_patches);

    // calculate nu_turb_field from the k and epsilon fields
    nu_turb_field.resize(num_cells);
    nu_turb_bcs = std::vector<BoundaryConditions> (num_patches, BoundaryConditions::NEUMANN);
    nu_turb_bc_values.resize(num_patches);

    // initialise flow constants to standard values
    C1 = 1.44;
    C2 = 1.92;
    C_nu = 0.09;
    sigma_k = 1;
    sigma_epsilon = 1.3;

    // constants used elsewhere
    // I - calc_k_bcs - turbulence intensity estimate

    kappa = 0.41;
    E = 9.78;

    inertial_sublayer_min = 11.53;
    inertial_sublayer_max = 250;

    k_underrel_coeff = 0.6;
    epsilon_underrel_coeff = 0.6;

    k_intpol_method = InterpolationMethods::LINEAR;

    // reset Reynolds number
    Re = 22000;

    // resize all bc_values vectors to respective number of faces before writing into them
    for(int patch_id = 0; patch_id < num_patches; patch_id++)
    {
        std::size_t num_patch_faces = patches[patch_id].get_face_ids().size();

        k_bc_values[patch_id].resize(num_patch_faces);
        epsilon_bc_values[patch_id].resize(num_patch_faces);
        nu_turb_bc_values[patch_id].resize(num_patch_faces);
    }

    std::cerr<<"---------------------------------------------"<<std::endl;
    std::cerr<<"K-EPSILON WARNINGS:"<<std::endl;
    std::cerr<<std::endl;
    std::cerr<<"K-EPSILON UNDERRELAXATION COEFFICIENTS SET ARBITRARILY TO 0.6!"<<std::endl;
    std::cerr<<std::endl;
    std::cerr<<"calc_k_bcs() AND calc_epsilon_bcs() CONTAIN MAGIC NUMBERS!"<<std::endl;
    std::cerr<<std::endl;
    std::cerr<<"nu_field INITIALISATIONS INDEPENDENT OF MESH GEOMETRY!"<<std::endl;
    std::cerr<<"k AND epsilon INITIALISATIONS INDEPENDENT OF MESH GEOMETRY!"<<std::endl;
    std::cerr<<std::endl;
    std::cerr<<"WALL FUNCTION ONLY APPLIED TO CELLS WITH BOUNDARY FACES!"<<std::endl;
    std::cerr<<std::endl;
    std::cerr<<"BOUNDARY FACE WEIGHTS COMPUTED AT THE BEGINNING OF THE SIMULATION!"<<std::endl;
    std::cerr<<"ALL WALL FUNCTION VALUES ASSUMED TO HAVE EQUAL WEIGHTS!"<<std::endl;
    std::cerr<<"---------------------------------------------"<<std::endl;
}



void Turbulence2D::set_conds()
{
    // set conditions from parent class
    PVCoupling2D::set_conds();

    // ------------------- SET K AND EPSILON BOUNDARY CONDITIONS -------------------

    calc_k_bcs();       // dependent on velocity
    calc_epsilon_bcs(); // dependent on k

    // ------------------- INITIALISE K FIELD -------------------

    // get k dirichlet boundary value
    double k_inlet_average;
    for(int patch_id = 0; patch_id < num_patches; patch_id++)
    {
        if(patches[patch_id].get_name().find("Inlet") != std::string::npos)
        {
            k_inlet_average = calc_field_average(k_bc_values[patch_id]);
            break;
        }
    }

    // set initial values for k field
    k_field = std::vector<double> (num_cells, k_inlet_average);

    // ------------------- INITIALISE EPSILON FIELD -------------------

    // get epsilon dirichlet boundary value
    double epsilon_inlet_average;
    for(int patch_id = 0; patch_id < num_patches; patch_id++)
    {
        if(patches[patch_id].get_name().find("Inlet") != std::string::npos)
        {
            epsilon_inlet_average = calc_field_average(epsilon_bc_values[patch_id]);
            break;
        }
    }

    // set initial values for epsilon field
    epsilon_field = std::vector<double> (num_cells, epsilon_inlet_average);

    // ------------------- COMPUTE NU_TURB_FIELD -------------------
    calc_nu_turb_field();

    // ------------------- POPULATE WALL_CELL_WEIGHTS -------------------
    calc_wall_cell_weights();
}



// protected

void Turbulence2D::calc_wall_cell_weights()
{
    // initialise
    for(int patch_id = 0; patch_id < num_patches; patch_id++)
    {
        if(patches[patch_id].get_name().find("Wall") != std::string::npos)   // if this is a wall
        {
            for(int face_id: patches[patch_id].get_face_ids())  // loop through the faces
            {
                // get owner and face data
                int face_owner_id = faces[face_id].get_owner_id();

                wall_cell_weights[face_owner_id] = 0;
            }
        }
    }

    // calculate occurences
    for(int patch_id = 0; patch_id < num_patches; patch_id++)
    {
        if(patches[patch_id].get_name().find("Wall") != std::string::npos)   // if this is a wall
        {
            for(int face_id: patches[patch_id].get_face_ids())  // loop through the faces
            {
                // get owner and face data
                int face_owner_id = faces[face_id].get_owner_id();

                wall_cell_weights[face_owner_id]++;
            }
        }
    }

    // calculate weights by inverting occurences (equal weighting)
    for(auto &member: wall_cell_weights)
        member.second = 1./member.second;
}