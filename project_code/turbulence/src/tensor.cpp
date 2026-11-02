#include <cassert>

#include "point.hpp"
#include "tensor.hpp"

template<>
Tensor<Point> Tensor<Point>::transpose()
{
    // only possible if the tensor components have the same dimension as Point
    assert(components.size() == (size_t) N_DIMS);

    // more initialisations
    std::array<std::array<double, N_DIMS>, N_DIMS> coordinate_matrix;
    std::vector<Point> transposed_points;
    transposed_points.resize(N_DIMS);

    // extract coordinates
    for(int i=0; i<N_DIMS; i++)
        coordinate_matrix[i] = components[i].get_coordinates();

    // swap components
    for(int i=0; i<N_DIMS; i++)
    {
        std::array<double, N_DIMS> transposed_row;
        for(int j=0; j<N_DIMS; j++)
            transposed_row[j] = coordinate_matrix[j][i];

        transposed_points[i] = Point(transposed_row);
    }

    return Tensor(transposed_points);
}

template<>
void Tensor<Point>::add(Tensor<Point>& tensor)
{
    std::vector<Point> tensor_components = tensor.get_components();

    size_t num_components = components.size();

    for(size_t i=0; i<num_components; i++)
        components[i].move_by(tensor_components[i]);
}

template<>
void Tensor<Point>::scale_by(double scaling_factor)
{
    size_t num_components = components.size();

    for(size_t i=0; i<num_components; i++)
        components[i].scale_by(scaling_factor);
}

template<>
double Tensor<Point>::double_dot(Tensor<Point>& tensor)
{
    // initialisations
    std::vector<Point> tensor_components = tensor.get_components();
    size_t num_components = components.size();
    double result = 0;

    for(size_t i=0; i<num_components; i++)
        result += components[i].dot_product_with(tensor_components[i]);

    return result;
}