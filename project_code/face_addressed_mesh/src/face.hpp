#include <array>
#include <vector>

#include "point.hpp"

#ifndef FACE
#define FACE

class Face
{
    std::vector<int> point_ids;
    int owner_id;
    int neighbour_id;
    double area;
    Point area_vector;
    Point centroid;

    public:
        Face(std::vector<int>& point_ids);

        std::vector<int> get_point_ids();
        int get_owner_id();
        int get_neighbour_id();
        double get_area();
        Point get_area_vector();
        Point get_centroid();

        void set_area_vector(Point& area_vector);
        void set_centroid(Point& centroid);
        void set_owner_id(int owner_id);
        void set_neighbour_id(int neighbour_id);
};

#endif