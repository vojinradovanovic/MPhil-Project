#include "pvCoupling2D.hpp"

int main()
{
    PVCoupling2D myProblem;
    myProblem.initialise_params();
    myProblem.set_conds();

    myProblem.solve_diffusion();

    return 0;
}