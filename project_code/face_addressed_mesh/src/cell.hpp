#include <vector>

#include "point.hpp"

#ifndef CELL
#define CELL

class Cell
{
    std::vector<int> face_ids;
    double volume;
    Point centroid;

    public:
        Cell(std::vector<int>& face_ids);

        std::vector<int> get_face_ids();
        double get_volume();
        Point get_centroid();

        void set_volume(double volume);
        void set_centroid(Point& centroid);
};

#endif