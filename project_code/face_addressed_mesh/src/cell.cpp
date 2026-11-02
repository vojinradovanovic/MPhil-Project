#include <array>
#include <cassert>
#include <cmath>
#include <vector>

#include "cell.hpp"

Cell::Cell(std::vector<int>& face_ids)
{
    this->face_ids = face_ids;
}



std::vector<int> Cell::get_face_ids()
{
    return face_ids;
}

double Cell::get_volume()
{
    return volume;
}

Point Cell::get_centroid()
{
    return centroid;
}



void Cell::set_volume(double volume)
{
    this->volume = volume;
}

void Cell::set_centroid(Point& centroid)
{
    this->centroid = centroid;
}