#include <array>
#include <limits>

#ifndef CONSTANTS
#define CONSTANTS

enum Dimensions {X, Y, Z};

const int N_DIMS = 2;
const std::array<Dimensions, N_DIMS> ALL_DIMS {X, Y};
const double EPSILON = 1e-8;

// 2D-specific constants
const int NUM_FACE_POINTS_2D = 2;
const int POINT_COUNT_VAL_2D = 2;

// 3D-specific constants
const int NUM_EDGE_POINTS_3D = NUM_FACE_POINTS_2D;
const int EDGE_COUNT_VAL_3D = POINT_COUNT_VAL_2D;

#endif