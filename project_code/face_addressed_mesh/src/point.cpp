#include <array>
#include <cmath>

#include "point.hpp"

Point::Point()
{
    for(Dimensions dim: ALL_DIMS)
        coords[dim] = 0;
}

Point::Point(std::array<double, N_DIMS>& coords)
{
    this->coords = coords;
}



std::array<double, N_DIMS> Point::get_coordinates()
{
    return coords;
}

void Point::set_coordinates(std::array<double, N_DIMS>& coords)
{
    this->coords = coords;
}

void Point::move_by(Point& p)
{
    std::array<double, N_DIMS> p_coords = p.get_coordinates();

    for(Dimensions dim: ALL_DIMS)
        coords[dim] += p_coords[dim];
}

void Point::scale_by(double c)
{
    for(Dimensions dim: ALL_DIMS)
        coords[dim] *= c;
}

void Point::negate()
{
    scale_by(-1);
}



double Point::calc_modulus()
{
    double modulus_sq = 0;

    for(Dimensions dim: ALL_DIMS)
        modulus_sq += coords[dim]*coords[dim];

    return sqrt(modulus_sq);
}

double Point::dot_product_with(Point& p)
{
    // initialisations
    std::array<double, N_DIMS> p_coords = p.get_coordinates();
    double dot_product = 0;

    for(Dimensions dim: ALL_DIMS)
        dot_product += coords[dim] * p_coords[dim];

    return dot_product;
}

Point Point::distance(Point& a, Point& b)
{
    // distance from a to b is equivalent to b-a
    Point dist = a;
    dist.negate();
    dist.move_by(b);

    return dist;
}