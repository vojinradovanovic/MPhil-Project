#include <array>

#include "physics_constants.hpp"

#ifndef COMPRESSED_ROW_MATRIX
#define COMPRESSED_ROW_MATRIX

struct MatrixElmt
{
    double value;
    std::array<int, MATRIX_DIMS> indices;
    MatrixElmt* next_elmt = nullptr;
};

class CompressedRowMatrix
{
    MatrixElmt* first_elmt;
    MatrixElmt* last_elmt;

    public:
        CompressedRowMatrix();
        
        void set_elmt(MatrixElmt* elmt); // used both for adding new and changing existing elements
        void add_elmt(MatrixElmt* elmt);
        void remove_row(int index);
        void clear();

        MatrixElmt* get_first_elmt();
        double extract_value_remove_elmt(std::array<int, MATRIX_DIMS>& indices);

        void push_elmt(MatrixElmt* elmt);
};

#endif