#include <array>
#include <iostream>

#include "linear_system.hpp"

int main()
{
    LinearSystem test(3);

    /* 3x3 matrix:
     * 6 2 -4   1   -2
     * 5 6 -2 * 2 = 11
     * 5 2 10   3   0
     */

    std::array<std::array<double, 3>, 3> matrix;
    matrix[0] = std::array<double, 3>{6, 2, -4};
    matrix[1] = std::array<double, 3>{5, 6, -2};
    matrix[2] = std::array<double, 3>{5, 2, 10};

    std::array<double, 3> rhs {-2, 11, 39};

    for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    {
        if(i==j)
        {
            test.add_to_diagonal_coeffs(matrix[i][i], i);
            test.add_to_rhs_vector(rhs[i], i);
        }
        else
        {
            std::array<int, 2> indices{i, j};
            test.add_to_off_diagonal_coeffs(matrix[i][j], indices);
        }
    }

    std::vector<double> x (3, 0);
    test.solve_Gauss_Seidel(x);

    std::cout<<"------- solution: ------------"<<std::endl;

    for(double elem: x)
        std::cout<<elem<<std::endl;
}