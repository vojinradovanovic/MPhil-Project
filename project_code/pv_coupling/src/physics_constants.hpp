#include <vector>

#include "constants.hpp"
#include "point.hpp"

#ifndef PHYSICS_CONSTANTS
#define PHYSICS_CONSTANTS

const int MATRIX_DIMS = 2;  // solver matrices are flattened to 2-dimensional tensors
const double INNER_REL_TOL = 1e-3;
const double OUTER_REL_TOL = 1e-5;
const double FLUX_TOL = 1e-8;
const int N_ITERATIONS = 100;
const int N_TIME_STEPS = 20;

enum class InterpolationMethods {LINEAR, UPWIND};

enum class BoundaryConditions {NONE, DIRICHLET, NEUMANN, DIRICHLET_VARIABLE};

const double D_T = 0.01; // time step

// default boundary conditions;
// constants are defined as functions for uniformity
inline double P_DIRICHLET()
{
    return 1;
}
inline double DP_DN_NEUMANN()
{
    return 0;
}

// these are constants because Point constructors do not take constant arrays
// but if they are not constants then the program returns multiple definition errors,
// even with inline
inline Point VEL_DIRICHLET()
{
    std::array<double, N_DIMS> vel_dirichlet_coords{1, 0};

    return Point(vel_dirichlet_coords);
}
inline Point DVEL_DN_NEUMANN()
{
    std::array<double, N_DIMS> dvel_dn_neumann_coords{0, 0};

    return Point(dvel_dn_neumann_coords);
}

#endif