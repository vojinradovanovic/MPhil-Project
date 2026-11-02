#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

#include "linear_system.hpp"

//public

LinearSystem::LinearSystem()    // default constructor
{
    *this = LinearSystem(0);
}

LinearSystem::LinearSystem(int size)
{
    // off-diagonal matrix coefficients are initialised to null by default,
    // as well as the raw coefficients

    // set size
    this->size = size;

    // resize vector parameters
    rhs_vector.resize(size);
    diagonal_coeffs.resize(size);
}



void LinearSystem::add_to_rhs_vector(double value, int index)
{
    rhs_vector[index] += value;
}

void LinearSystem::add_to_diagonal_coeffs(double value, int index)
{
    diagonal_coeffs[index] += value;
}

void LinearSystem::add_to_off_diagonal_coeffs(double value, const std::array<int, MATRIX_DIMS>& indices)
{
    // MatrixElmt *elmt = new MatrixElmt;

    // elmt->value = value;
    // elmt->indices = indices;
    // elmt->next_elmt = nullptr;

    RawMatrixElmt raw_elmt;
    raw_elmt.value = value;
    raw_elmt.indices = indices;

    // assign element to off diagonal coefficients
    if(indices[0] < indices[1])
    {
        // upper_coeffs.add_elmt(elmt);

        upper_coeffs_raw.push_back(raw_elmt);
    }

    else if(indices[0] > indices[1])
    {
        // lower_coeffs.add_elmt(elmt);

        lower_coeffs_raw.push_back(raw_elmt);
    }

    else
        std::cerr<<"USE add_to_diagonal_coeffs() FOR ELEMENT ("<<indices[0]<<" "<<indices[1]<<")!"<<std::endl;
}



void LinearSystem::scale_diagonal_coeff(double scaling_factor, int index)
{
    diagonal_coeffs[index] *= scaling_factor;
}



void LinearSystem::set_const_row(double value, int index, std::vector<int>& neighbour_indices)
{
    // safety check
    assert(index < size);

    // implement a row, in which the equation will be diagonal_coeffs[index]*x[index] = diagonal_coeffs[index]*value
    rhs_vector[index] = diagonal_coeffs[index] * value;

    for(int neighbour_id: neighbour_indices)
    {
        std::array<int, MATRIX_DIMS> indices = {neighbour_id, index};

        if(neighbour_id < index)
        {
            // extract coefficient and remove it from the matrix
            double upper_coeff = upper_coeffs.extract_value_remove_elmt(indices);

            // move upper_coeff to rhs
            rhs_vector[neighbour_id] -= upper_coeff * value;
        }
        else if(neighbour_id > index)
        {
            // extract coefficient and remove it from the matrix
            double lower_coeff = lower_coeffs.extract_value_remove_elmt(indices);

            // move lower_coeff to rhs
            rhs_vector[neighbour_id] -= lower_coeff * value;
        }
        else
        {
            std::cerr<<"CELL "<<index<<" IS A NEIGHBOUR OF ITSELF!"<<std::endl;
            std::abort();
        }
    }

    // delete upper_coeffs and lower_coeffs in index-th row
    upper_coeffs.remove_row(index);
    lower_coeffs.remove_row(index);
}

void LinearSystem::overload_row(double value, int index)
{
    // safety check
    assert(index < size);

    /* for a row:
     *      diagonal_coeffs[index] * x
     *          + sum(upper_coeffs * x_upper) + sum(lower_coeffs * x_lower)
     *          = rhs_vector[index],
     * set to:
     *      (diagonal_coeffs[index] + weight) * x
     *          + sum(upper_coeffs * x_upper) + sum(lower_coeffs * x_lower)     (no change)
     *          = (diagonal_coeffs[index] + weight) * value,
     */

    double weight = 10000 * diagonal_coeffs[index];

    diagonal_coeffs[index] = diagonal_coeffs[index] + weight;
    rhs_vector[index] = diagonal_coeffs[index] * value;
}



std::vector<double>& LinearSystem::get_diagonal_coeffs()
{
    return diagonal_coeffs;
}

std::vector<double>& LinearSystem::get_rhs_vector()
{
    return rhs_vector;
}

CompressedRowMatrix& LinearSystem::get_upper_coeffs()
{
    return upper_coeffs;
}

CompressedRowMatrix& LinearSystem::get_lower_coeffs()
{
    return lower_coeffs;
}

void LinearSystem::query_cell(int index)
{
    std::cout<<"matrix["<<index<<","<<index<<"]:\t"
                <<diagonal_coeffs[index]<<std::endl;

    MatrixElmt *upper_elmt = upper_coeffs.get_first_elmt();
    MatrixElmt *lower_elmt = lower_coeffs.get_first_elmt();

    while(upper_elmt)
    {
        if(upper_elmt->indices[0] > index)
            break;
        else if(upper_elmt->indices[0] == index)
        {
            std::cout<<"matrix["
                <<upper_elmt->indices[0]
                <<","
                <<upper_elmt->indices[1]
                <<"]:\t"
                <<upper_elmt->value<<std::endl;

            upper_elmt = upper_elmt->next_elmt;
        }
        else
            upper_elmt = upper_elmt->next_elmt;
    }

    while(lower_elmt)
    {
        if(lower_elmt->indices[0] > index)
            break;
        else if(lower_elmt->indices[0] == index)
        {
            std::cout<<"matrix["
                <<lower_elmt->indices[0]
                <<","
                <<lower_elmt->indices[1]
                <<"]:\t"
                <<lower_elmt->value<<std::endl;

            lower_elmt = lower_elmt->next_elmt;
        }
        else
            lower_elmt = lower_elmt->next_elmt;
    }

    std::cout<<"rhs["<<index<<"]:\t"
                <<rhs_vector[index]<<std::endl;
}



void LinearSystem::solve_Gauss_Seidel(std::vector<double>& x, double rel_tol)
{
    // declarations
    double residual_norm, init_residual_norm;
    int iter;
    int MAX_GS_ITERATIONS = 5000;
    int checking_period = 100;      // every period-th iteration, algorithm is checked for convergence

    // calculate residual at the beginning
    init_residual_norm = calc_residual(x);
    residual_norm = init_residual_norm;
    std::cout<<"Before solving, residual norm is "<<residual_norm<<std::endl;

    // determine whether convergence is reached
    bool convergence_reached = (init_residual_norm == 0);

    // solve
    for(iter=0; !convergence_reached && iter<MAX_GS_ITERATIONS; iter++)    // check tolerance every 100 iterations
    {
        // update x through upper_coeffs and lower_coeffs
        MatrixElmt *upper_elmt = upper_coeffs.get_first_elmt();
        MatrixElmt *lower_elmt = lower_coeffs.get_first_elmt();

        for(int i=0; i<size; i++)
        {
            // start initialising new values with rhs_vector
            x[i] = rhs_vector[i];

            // loop through upper_coeffs and update x[i] with i+1...
            while(upper_elmt)
            {
                std::array<int, MATRIX_DIMS> indices = upper_elmt->indices;

                if(indices[0] == i)
                {
                    x[i] -= upper_elmt->value * x[indices[1]];
                    upper_elmt = upper_elmt->next_elmt;
                }
                else
                    break;
            }

            // loop through lower_coeffs and update x[i] with 0...i-1
            while(lower_elmt)
            {
                std::array<int, MATRIX_DIMS> indices = lower_elmt->indices;

                if(indices[0] == i)
                {
                    x[i] -= lower_elmt->value * x[indices[1]];
                    lower_elmt = lower_elmt->next_elmt;
                }
                else
                    break;
            }

            // update x[i] with diagonal_coeffs
            if(diagonal_coeffs[i] == 0)
                assert(rhs_vector[i] == 0);
            else
                x[i] /= diagonal_coeffs[i];
        }

        // update convergence_reached and residual_norm every checking_period iterations
        if(iter % checking_period == checking_period-1)
        {
            residual_norm = calc_residual(x);
            convergence_reached = (residual_norm / init_residual_norm <= rel_tol);
        }
    }

    // recalculate to get the final residual
    residual_norm = calc_residual(x);

    // output residual at the end
    std::cout<<"After "<<iter<<" iterations, residual norm is "<<residual_norm<<std::endl;
}



void LinearSystem::clear()
{
    // size stays the same

    // clear vectors
    rhs_vector = std::vector<double> (size, 0);
    diagonal_coeffs = std::vector<double> (size, 0);

    // clear CompressedRowMatrix's
    upper_coeffs.clear();
    lower_coeffs.clear();

    // clear raw vectors
    upper_coeffs_raw.clear();
    lower_coeffs_raw.clear();
}



double LinearSystem::calc_residual(const std::vector<double>& x, int index)
{
    // pre-check for validity
    assert(x.size() == (size_t) size);

    double residual_norm = 0;
    MatrixElmt *upper_elmt = upper_coeffs.get_first_elmt();
    MatrixElmt *lower_elmt = lower_coeffs.get_first_elmt();

    // calculate residual_norm
    // initialise temp_var with rhs vector
    double temp_var = rhs_vector[index];

    // loop through upper_coeffs and temp_var with i+1...
    while(upper_elmt)
    {
        std::array<int, MATRIX_DIMS> indices = upper_elmt->indices;

        if(indices[0] == index)     // if upper_elmt is in index-th row
        {
            temp_var -= upper_elmt->value * x[indices[1]];
            upper_elmt = upper_elmt->next_elmt;
        }
        else if(indices[0] < index) // if upper_elmt is above index-th row
        {
            upper_elmt = upper_elmt->next_elmt; // iterate to get to index-th row
        }
        else                        // if upper_elmt is below index-th row
            break;
    }

    // loop through lower_coeffs and update temp_var with 0...i-1
    while(lower_elmt)
    {
        std::array<int, MATRIX_DIMS> indices = lower_elmt->indices;

        if(indices[0] == index)     // if lower_elmt is in index-th row
        {
            temp_var -= lower_elmt->value * x[indices[1]];
            lower_elmt = lower_elmt->next_elmt;
        }
        else if(indices[0] < index) // if lower_elmt is above index-th row
        {
            lower_elmt = lower_elmt->next_elmt; // iterate to get to index-th row
        }
        else                        // if lower_elmt is below index-th row
            break;
    }

    // update with diag_coeffs
    temp_var -= diagonal_coeffs[index] * x[index];

    // update residual norm
    residual_norm = fabs(temp_var);

    return residual_norm;
}

double LinearSystem::calc_residual(const std::vector<double>& x)
{
    // pre-check for validity
    assert(x.size() == (size_t) size);

    double residual_norm = 0;
    MatrixElmt *upper_elmt = upper_coeffs.get_first_elmt();
    MatrixElmt *lower_elmt = lower_coeffs.get_first_elmt();

    // calculate residual_norm
    for(int i=0; i<size; i++)
    {
        // initialise temp_var with rhs vector
        double temp_var = rhs_vector[i];

        // loop through upper_coeffs and temp_var with i+1...
        while(upper_elmt)
        {
            std::array<int, MATRIX_DIMS> indices = upper_elmt->indices;

            if(indices[0] == i)     // if upper_elmt is in i-th row
            {
                temp_var -= upper_elmt->value * x[indices[1]];
                upper_elmt = upper_elmt->next_elmt;
            }
            else if(indices[0] < i) // if upper_elmt is above i-th row
            {
                upper_elmt = upper_elmt->next_elmt; // iterate to get to i-th row
            }
            else                        // if upper_elmt is below i-th row
                break;
        }

        // loop through lower_coeffs and update temp_var with 0...i-1
        while(lower_elmt)
        {
            std::array<int, MATRIX_DIMS> indices = lower_elmt->indices;

            if(indices[0] == i)     // if lower_elmt is in i-th row
            {
                temp_var -= lower_elmt->value * x[indices[1]];
                lower_elmt = lower_elmt->next_elmt;
            }
            else if(indices[0] < i) // if lower_elmt is above i-th row
            {
                lower_elmt = lower_elmt->next_elmt; // iterate to get to i-th row
            }
            else                        // if lower_elmt is below i-th row
                break;
        }

        // update with diag_coeffs
        temp_var -= diagonal_coeffs[i] * x[i];

        // update residual norm
        residual_norm += fabs(temp_var);
    }

    return residual_norm;
}

std::vector<double> LinearSystem::calc_H_operator(const std::vector<double>& x)
{
    // pre-check for validity
    assert(x.size() == (size_t) size);

    // initialisations
    MatrixElmt *upper_elmt = upper_coeffs.get_first_elmt();
    MatrixElmt *lower_elmt = lower_coeffs.get_first_elmt();
    std::vector<double> H_values = rhs_vector;

    // subtract coefficients from upper
    while(upper_elmt)
    {
        H_values[upper_elmt->indices[0]] -= upper_elmt->value * x[upper_elmt->indices[1]];

        upper_elmt = upper_elmt->next_elmt;
    }

    // subtract coefficients from lower
    while(lower_elmt)
    {
        H_values[lower_elmt->indices[0]] -= lower_elmt->value * x[lower_elmt->indices[1]];

        lower_elmt = lower_elmt->next_elmt;
    }

    return H_values;
}



void LinearSystem::print_out()
{
    // initialisation
    MatrixElmt *upper_elmt = upper_coeffs.get_first_elmt();
    MatrixElmt *lower_elmt = lower_coeffs.get_first_elmt();

    std::cout<<"-------------- SYSTEM --------------"<<std::endl;

    for(int i=0; i<size; i++)
    {
        std::cout<<"Cell "<<i<<std::endl;

        // print matrix coefficients
        for(int j=0; j<size; j++)
        {
            if(i == j)
            {
                std::cout<<diagonal_coeffs[i]<<"\t";
            }
            else if(i<j)    // upper
            {
                // check for existence
                if(upper_elmt)
                {
                    std::array<int, MATRIX_DIMS> indices = upper_elmt->indices;

                    if(indices[0] == i && indices[1] == j)
                    {
                        std::cout<<upper_elmt->value<<"\t";
                        upper_elmt = upper_elmt->next_elmt;
                    }
                    else
                        std::cout<<"0\t";
                }
                else
                    std::cout<<"0\t";
            }
            else if(i>j)    // lower
            {
                // check for existence
                if(lower_elmt)
                {
                    std::array<int, MATRIX_DIMS> indices = lower_elmt->indices;

                    if(indices[0] == i && indices[1] == j)
                    {
                        std::cout<<lower_elmt->value<<"\t";
                        lower_elmt = lower_elmt->next_elmt;
                    }
                    else
                        std::cout<<"0\t";
                }
                else
                    std::cout<<"0\t";
            }

            if(j>0 && j%10 == 9)
                std::cout<<"\n";
        }

        // print x
        std::cout<<"*\n"<<"x\n"<<"=\n";

        // print rhs vector
        std::cout<<rhs_vector[i]<<std::endl<<std::endl;
    }

    std::cout<<"-------------- END --------------"<<std::endl;
}

void LinearSystem::print_out(int index)
{
    // initialisation
    MatrixElmt *upper_elmt = upper_coeffs.get_first_elmt();
    MatrixElmt *lower_elmt = lower_coeffs.get_first_elmt();

    // get upper and lower element to current row
    while(upper_elmt)
    {
        if(upper_elmt->indices[0] < index)
            upper_elmt=upper_elmt->next_elmt;
        else
            break;
    }
    while(lower_elmt)
    {
        if(lower_elmt->indices[0] < index)
            lower_elmt=lower_elmt->next_elmt;
        else
            break;
    }

    std::cout<<"Cell "<<index<<std::endl;

    // print matrix coefficients
    for(int j=0; j<size; j++)
    {
        if(index == j)
        {
            std::cout<<diagonal_coeffs[index]<<"\t";
        }
        else if(index<j)    // upper
        {
            // check for existence
            if(upper_elmt)
            {
                std::array<int, MATRIX_DIMS> indices = upper_elmt->indices;

                if(indices[0] == index && indices[1] == j)
                {
                    std::cout<<upper_elmt->value<<"\t";
                    upper_elmt = upper_elmt->next_elmt;
                }
                else
                    std::cout<<"0\t";
            }
            else
                std::cout<<"0\t";
        }
        else if(index>j)    // lower
        {
            // check for existence
            if(lower_elmt)
            {
                std::array<int, MATRIX_DIMS> indices = lower_elmt->indices;

                if(indices[0] == index && indices[1] == j)
                {
                    std::cout<<lower_elmt->value<<"\t";
                    lower_elmt = lower_elmt->next_elmt;
                }
                else
                    std::cout<<"0\t";
            }
            else
                std::cout<<"0\t";
        }

        if(j>0 && j%10 == 9)
            std::cout<<"\n";
    }

    // print x
    std::cout<<"*\n"<<"x\n"<<"=\n";

    // print rhs vector
    std::cout<<rhs_vector[index]<<std::endl<<std::endl;
}



void LinearSystem::reserve_coeffs_raw(int reserve_size)
{
    upper_coeffs_raw.reserve(reserve_size);
    lower_coeffs_raw.reserve(reserve_size);
}

void LinearSystem::populate_off_diagonal_coefficients()
{
    sort_raw();
    populate_coeffs_from_vector(upper_coeffs, upper_coeffs_raw);
    populate_coeffs_from_vector(lower_coeffs, lower_coeffs_raw);
}

void LinearSystem::sort_raw()
{
    std::sort(upper_coeffs_raw.begin(), upper_coeffs_raw.end());
    std::sort(lower_coeffs_raw.begin(), lower_coeffs_raw.end());
}

void LinearSystem::populate_coeffs_from_vector(CompressedRowMatrix& coeffs, std::vector<RawMatrixElmt>& coeffs_raw)
{
    if(!coeffs.get_first_elmt())
    {
        size_t i=0;
        while(i < coeffs_raw.size())
        {
            MatrixElmt* elmt = new MatrixElmt;
            elmt->value = coeffs_raw[i].value;
            elmt->indices = coeffs_raw[i].indices;
            elmt->next_elmt = nullptr;

            // loop through next few elements to collapse them into one
            size_t j = i+1;
            while(j < coeffs_raw.size() && coeffs_raw[i].indices == coeffs_raw[j].indices)
            {
                elmt->value += coeffs_raw[j].value;
                j++;
            }

            // move iterator
            i = j;

            if(fabs(elmt->value) > EPSILON)
                coeffs.push_elmt(elmt);
            else
                delete elmt;
        }
    }
}