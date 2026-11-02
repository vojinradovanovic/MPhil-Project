#include <array>

#include "constants.hpp"

#ifndef POINT
#define POINT

class Point
{
    std::array<double, N_DIMS> coords;

    public:
        Point();
        Point(std::array<double, N_DIMS>& coords);
        
        std::array<double, N_DIMS> get_coordinates();
        void set_coordinates(std::array<double, N_DIMS>& coords);
        void move_by(Point& p);
        void scale_by(double c);
        void negate();

        double calc_modulus();
        double dot_product_with(Point& p);
        static Point distance(Point& a, Point& b);
};

#endif