#include <array>
#include <cassert>
#include <cmath>
#include <vector>

#include "face.hpp"

Face::Face(std::vector<int>& point_ids)
{
    this->point_ids = point_ids;

    // initialise owner_id and neighbour_id to -1;
    owner_id = -1;
    neighbour_id = -1;
}



std::vector<int> Face::get_point_ids()
{
    return point_ids;
}

int Face::get_owner_id()
{
    return owner_id;
}

int Face::get_neighbour_id()
{
    return neighbour_id;
}

double Face::get_area()
{
    return area;
}

Point Face::get_area_vector()
{
    return area_vector;
}

Point Face::get_centroid()
{
    return centroid;
}



void Face::set_area_vector(Point& area_vector)
{
    this->area_vector = area_vector;
    area = area_vector.calc_modulus();
}

void Face::set_centroid(Point& centroid)
{
    this->centroid = centroid;
}

void Face::set_owner_id(int owner_id)
{
    this->owner_id = owner_id;
}

void Face::set_neighbour_id(int neighbour_id)
{
    this->neighbour_id = neighbour_id;
}