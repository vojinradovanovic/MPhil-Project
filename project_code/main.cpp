#include "turbulence2D.hpp"

int main()
{
    Turbulence2D myProblem;
    myProblem.initialise_params();

    myProblem.set_conds();
    myProblem.solve_k_epsilon();
}