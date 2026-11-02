#include <array>
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "turbulence2D.hpp"

// protected

void Turbulence2D::write_output(std::ofstream& output_file)
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
        double grid_scaling = 5;
        double velocity_scaling = (1 - 1./(2*velocity_modulus + 1)) * grid_scaling;

        // beginning coordinates
        for(Dimensions dim: ALL_DIMS)
            output_file<<centroid_coords[dim] - velocity_coords[dim]/2 / (velocity_modulus+EPSILON) * velocity_scaling<<"\t";

        // ending coordinates
        for(Dimensions dim: ALL_DIMS)
            output_file<<velocity_coords[dim] / (velocity_modulus+EPSILON) * velocity_scaling<<"\t";

        // velocity components
        for(Dimensions dim: ALL_DIMS)
            output_file << velocity_coords[dim] << "\t";

        // velocity magnitude
        output_file<<velocity_field[cell_id].calc_modulus()<<"\t";

        // write k
        output_file<<k_field[cell_id]<<"\t";

        // write epsilon
        output_file<<epsilon_field[cell_id]<<"\t";

        // write turbulent viscosity
        output_file<<nu_turb_field[cell_id]<<"\t";

        output_file<<std::endl;
    }

    output_file<<std::endl<<std::endl;
}
