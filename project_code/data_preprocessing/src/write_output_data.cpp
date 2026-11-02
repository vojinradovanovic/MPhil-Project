#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "helpers.hpp"
#include "turbulence2D.hpp"

int main()
{
    std::string input_filepath = "../../OpenFOAM/",
                input_simulation = "pitzDaily_pimple_base_case/0.3";
    std::vector<std::string> fields = {"p", "U", "k", "epsilon", "nut"};
    std::size_t num_fields = fields.size();

    std::string intermediate_filepath = "intermediate/";

    std::cout<<"Writing simulation for:"<<std::endl;
    std::cout<<input_simulation<<std::endl<<std::endl;

    // write to intermediate files
    for(std::size_t i=0; i<num_fields; i++)
    {
        std::ifstream input_file(input_filepath + input_simulation + "/" + fields[i]);
        std::ofstream intermediate_file(intermediate_filepath + fields[i]);

        if(fields[i] == "U")
            strip_vector_field(input_file, intermediate_file);
        else
            strip_scalar_field(input_file, intermediate_file);
    }

    // write to output file
    std::string output_data_filepath = "../../data/";

    FaceMesh2D myMesh;
    myMesh.initialise();

    std::ofstream output_file(output_data_filepath + "turbulence_end_pimple.dat");
    write_output_data(output_file, myMesh, intermediate_filepath, fields);

    // delete all input files in intermediate path
    clear_dir_contents(intermediate_filepath);

    return 0;
}