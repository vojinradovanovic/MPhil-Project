#include <array>
#include <vector>

#include "compressed_row_matrix.hpp"
#include "physics_constants.hpp"

#ifndef LINEAR_SYSTEM
#define LINEAR_SYSTEM

struct RawMatrixElmt
{
    double value;
    std::array<int, MATRIX_DIMS> indices;

    bool operator<(const RawMatrixElmt& other) const
    {
        return indices < other.indices;
    }
};

class LinearSystem
{
    int size;
    
    CompressedRowMatrix upper_coeffs, lower_coeffs;
    std::vector<double> diagonal_coeffs;
    std::vector<double> rhs_vector;

    // helper vector variables
    std::vector<RawMatrixElmt> upper_coeffs_raw, lower_coeffs_raw;

    public:
        LinearSystem();
        LinearSystem(int size);

        void add_to_rhs_vector(double value, int index);
        void add_to_diagonal_coeffs(double value, int index);
        void add_to_off_diagonal_coeffs(double value, const std::array<int, MATRIX_DIMS>& indices);

        void scale_diagonal_coeff(double scaling_factor, int index);

        void set_const_row(double value, int index, std::vector<int>& neighbour_indices);
        void overload_row(double value, int index);

        std::vector<double>& get_diagonal_coeffs();
        std::vector<double>& get_rhs_vector();
        CompressedRowMatrix& get_upper_coeffs();
        CompressedRowMatrix& get_lower_coeffs();
        void query_cell(int index);

        void solve_Gauss_Seidel(std::vector<double>& x, double rel_tol = INNER_REL_TOL);

        void clear();

        double calc_residual(const std::vector<double>& x);
        double calc_residual(const std::vector<double>& x, int index);
        std::vector<double> calc_H_operator(const std::vector<double>& x);

        void print_out();
        void print_out(int index);

        void reserve_coeffs_raw(int reserve_size);
        void populate_off_diagonal_coefficients();
        void sort_raw();
        void populate_coeffs_from_vector(CompressedRowMatrix& coeffs, std::vector<RawMatrixElmt>& coeffs_raw);
};

#endif