#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "helpers.hpp"

int main()
{
    std::string input_file = "../../data/turbulence_end.dat",
                input_file_pimple = "../../data/turbulence_end_pimple.dat";

    std::cout<<"Comparing files:"<<std::endl;
    std::cout<<input_file<<std::endl;
    std::cout<<input_file_pimple<<std::endl;

    std::cout<<std::endl;

    // write to output file
    std::string output_data_filepath = "../../data/";

    FaceMesh2D myMesh;
    myMesh.initialise();

    std::ofstream output_file(output_data_filepath + "comparison.dat");
    write_comparison(output_file, myMesh, input_file, input_file_pimple);

    return 0;
}