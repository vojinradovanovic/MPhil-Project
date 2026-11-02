#include <cassert>
#include <iostream>

#include "compressed_row_matrix.hpp"

CompressedRowMatrix::CompressedRowMatrix()
{
    first_elmt = nullptr;
    last_elmt = nullptr;
}



void CompressedRowMatrix::set_elmt(MatrixElmt* elmt)
{
    // set up prev_elmt as an element before first element
    MatrixElmt before_first;    // local variable, deleted at end of function
    before_first.next_elmt = first_elmt;
    MatrixElmt* prev_elmt = &before_first;

    // set up my_elmt from prev_elmt and flag to track allocation
    MatrixElmt* my_elmt = prev_elmt->next_elmt;
    bool alloc_flag = false;

    while(my_elmt)
    {
        for(int i=0; i<MATRIX_DIMS; i++)
        {
            // the assumption is that all elements are sorted in ascending order of indices:
            // for A[i0, i1, ...] -> A[j0, j1, ...] in the matrix,
            // we have (i0 == j0) \implies ... \implies (i_(k-1) == j_(k-1)) \implies (i_k < j_k) for k = 0, 1, ...

            // if elmt should be placed after my_elmt, exit loop
            if(elmt->indices[i] > my_elmt->indices[i])
                break;
            
            // if this is a strict <, my_elmt is the first element in ascending order that is after elmt;
            // insert elmt
            else if(elmt->indices[i] < my_elmt->indices[i])
            {
                // if my_elmt is still at first_elmt, move first_elmt to elmt
                if(my_elmt == first_elmt)
                    first_elmt = elmt;
                
                // insert elmt between prev_elmt and my_elmt
                elmt->next_elmt = my_elmt;
                prev_elmt->next_elmt = elmt;

                alloc_flag = true;
                break;
            }

            // if indices are equal, action cannot be determined;
            // continue to next indices[i]
            else
                continue;
        }

        if(my_elmt->indices == elmt->indices)   // edge case: if my_elmt and elmt point to same element
        {
            if(elmt->value == 0)    // if the value to be set is 0, then delete the element
            {
                prev_elmt->next_elmt = my_elmt->next_elmt;

                if (my_elmt == first_elmt)
                    first_elmt = my_elmt->next_elmt;

                if (my_elmt == last_elmt)
                    last_elmt = (first_elmt == nullptr) ? nullptr : prev_elmt;

                delete elmt;
                delete my_elmt;
            }
            else                    // otherwise accept the new value
            {
                my_elmt->value = elmt->value;   // also changes first_elmt->value and last_elmt->value
                                                // if implicated
                delete elmt;    // delete unused pointer
            }

            alloc_flag = true;
        }

        // evolve prev_elmt and my_elmt if unallocated
        if(alloc_flag)
            break;
        else
        {
            prev_elmt = my_elmt;
            my_elmt = my_elmt->next_elmt;
        }
    }

    // if prev_elmt has made it to last_elmt
    if(prev_elmt == last_elmt)
    {
        prev_elmt->next_elmt = elmt;
        last_elmt = elmt;

        last_elmt->next_elmt = nullptr;  // if allocated improperly

        alloc_flag = true;
    }
    else if(!alloc_flag)    // only if first_elmt==nullptr and last_elmt==nullptr;
    {
        assert(first_elmt==nullptr);
        assert(last_elmt==nullptr);

        // elmt is the first element to be allocated to the matrix
        first_elmt = elmt;
        last_elmt = elmt;

        last_elmt->next_elmt = nullptr;  // if allocated improperly

        alloc_flag = true;
    }

    assert(alloc_flag); // must have been allocated by end of function
}

void CompressedRowMatrix::add_elmt(MatrixElmt* elmt)
{
    // the code mostly follows the code for set_elmt,
    // with a few changes and additional comments

    // set up prev_elmt as an element before first element
    MatrixElmt before_first;    // local variable, deleted at end of function
    before_first.next_elmt = first_elmt;
    MatrixElmt* prev_elmt = &before_first;

    // set up my_elmt from prev_elmt and flag to track allocation
    MatrixElmt* my_elmt = prev_elmt->next_elmt;
    bool alloc_flag = false;

    while(my_elmt)
    {
        for(int i=0; i<MATRIX_DIMS; i++)
        {
            // the assumption is that all elements are sorted in ascending order of indices:
            // for A[i0, i1, ...] -> A[j0, j1, ...] in the matrix,
            // we have (i0 == j0) \implies ... \implies (i_(k-1) == j_(k-1)) \implies (i_k < j_k) for k = 0, 1, ...

            // if elmt should be placed after my_elmt, exit loop
            if(elmt->indices[i] > my_elmt->indices[i])
                break;
            
            // if this is a strict <, my_elmt is the first element in ascending order that is after elmt;
            // insert elmt
            else if(elmt->indices[i] < my_elmt->indices[i])
            {
                // if my_elmt is still at first_elmt, move first_elmt to elmt
                if(my_elmt == first_elmt)
                    first_elmt = elmt;
                
                // insert elmt between prev_elmt and my_elmt
                elmt->next_elmt = my_elmt;
                prev_elmt->next_elmt = elmt;

                alloc_flag = true;
                break;
            }

            // if indices are equal, action cannot be determined;
            // continue to next indices[i]
            else
                continue;
        }

        if(my_elmt->indices == elmt->indices)   // if my_elmt and elmt point to same element, increase value
        {
            my_elmt->value += elmt->value;  // also changes first_elmt->value and last_elmt->value
                                            // if implicated
            delete elmt;    // delete unused pointer

            alloc_flag = true;
        }

        // evolve prev_elmt and my_elmt if unallocated
        if(alloc_flag)
            break;
        else
        {
            prev_elmt = my_elmt;
            my_elmt = my_elmt->next_elmt;
        }
    }

    // if prev_elmt has made it to last_elmt
    if(prev_elmt == last_elmt)
    {
        prev_elmt->next_elmt = elmt;
        last_elmt = elmt;

        last_elmt->next_elmt = nullptr;  // if allocated improperly

        alloc_flag = true;
    }
    else if(!alloc_flag)    // only if first_elmt==nullptr and last_elmt==nullptr;
    {
        assert(first_elmt==nullptr);
        assert(last_elmt==nullptr);

        // elmt is the first element to be allocated to the matrix
        first_elmt = elmt;
        last_elmt = elmt;

        last_elmt->next_elmt = nullptr;  // if allocated improperly

        alloc_flag = true;
    }

    assert(alloc_flag); // must have been allocated by end of function
}

void CompressedRowMatrix::remove_row(int index)
{
    // set up prev_elmt as an element before first element
    MatrixElmt before_first;    // local variable, deleted at end of function
    before_first.next_elmt = first_elmt;
    MatrixElmt* prev_elmt = &before_first;

    // set up my_elmt from prev_elmt
    MatrixElmt* my_elmt = prev_elmt->next_elmt;

    while(my_elmt)
    {
        if(my_elmt->indices[0] < index)
        {
            // evolve prev_elmt and my_elmt
            prev_elmt = my_elmt;
            my_elmt = prev_elmt->next_elmt;
        }
        else if(my_elmt->indices[0] == index)
        {
            if(my_elmt == first_elmt)
                first_elmt = first_elmt->next_elmt;

            if(my_elmt == last_elmt)
                last_elmt = (first_elmt == nullptr) ? nullptr : prev_elmt;

            // evolve previous element
            prev_elmt->next_elmt = my_elmt->next_elmt;  // also performed in case previous if was true for consistency

            // delete my_elmt, then set the new value to prev_elmt->next_elmt
            delete my_elmt;
            my_elmt = prev_elmt->next_elmt;
        }
        else break;
    }
}

void CompressedRowMatrix::clear()
{
    MatrixElmt *elmt = first_elmt;
    MatrixElmt *next_elmt;

    while(elmt)
    {
        next_elmt = elmt->next_elmt;

        delete elmt;
        elmt = next_elmt;
    }

    first_elmt = nullptr;
    last_elmt = nullptr;
}



MatrixElmt* CompressedRowMatrix::get_first_elmt()
{
    return first_elmt;
}

double CompressedRowMatrix::extract_value_remove_elmt(std::array<int, MATRIX_DIMS>& indices)
{
    // set up prev_elmt as an element before first element
    MatrixElmt before_first;    // local variable, deleted at end of function
    before_first.next_elmt = first_elmt;
    MatrixElmt* prev_elmt = &before_first;

    // set up my_elmt from prev_elmt
    MatrixElmt* my_elmt = prev_elmt->next_elmt;

    while(my_elmt)
    {
        if(my_elmt->indices == indices)
        {
            // extract value from my_elmt
            double value = my_elmt->value;

            // make arrangements 
            prev_elmt->next_elmt = my_elmt->next_elmt;

            if (my_elmt == first_elmt)
                first_elmt = my_elmt->next_elmt;

            if (my_elmt == last_elmt)
                last_elmt = (first_elmt == nullptr) ? nullptr : prev_elmt;

            // free up space
            delete my_elmt;

            // return value and finish loop
            return value;
        }

        prev_elmt = my_elmt;
        my_elmt = my_elmt->next_elmt;
    }

    return 0;
}



void CompressedRowMatrix::push_elmt(MatrixElmt* elmt)
{
    if(!first_elmt || !last_elmt)
    {
        // first_elmt and last_elmt need to both be nullptr at the same time
        assert(first_elmt==nullptr);
        assert(last_elmt==nullptr);

        // elmt is the first element to be allocated to the matrix
        first_elmt = elmt;
        last_elmt = elmt;

        last_elmt->next_elmt = nullptr;  // if allocated improperly
    }
    else
    {
        last_elmt->next_elmt = elmt;
        last_elmt = elmt;

        last_elmt->next_elmt = nullptr;  // if allocated improperly
    }
}