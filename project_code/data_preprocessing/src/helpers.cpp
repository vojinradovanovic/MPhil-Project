#include <algorithm>
#include <filesystem>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "helpers.hpp"

void strip_scalar_field(std::ifstream& input_file, std::ofstream& output_file)
{
    std::string line;

    bool open_data = false;

    // char vectors for trimming
    std::vector<char> EOL_chars{'\r'};    // to accommodate for "\r\n" for new line in Windows

    while(std::getline(input_file, line))
    {
        line = rtrim(line, EOL_chars);

        if(open_data)   // if at data, write to output file
        {
            if(line == ")")
                break;
            else
            {
                output_file << line <<"\n";
            }
        }
        else
        {
            if(line == "(")
                open_data = true;
            else
                continue;
        }
    }

    input_file.close();
    output_file.close();
}

void strip_vector_field(std::ifstream& input_file, std::ofstream& output_file)
{
    std::string line;

    bool open_data = false;

    // char vectors for trimming
    std::vector<char> lchars{' ', '\t', '('},
                        rchars{'\r', ')'},
                        EOL_chars{'\r'};    // to accommodate for "\r\n" for new line in Windows

    while(std::getline(input_file, line))
    {
        line = rtrim(line, EOL_chars);

        if(open_data)   // if at data, write to output file
        {
            if(line == ")")
                break;
            else
            {
                line = trim(line, lchars, rchars);
                std::istringstream isstream(line);

                double coord1, coord2;
                isstream >> coord1 >> coord2;

                output_file << coord1 << " " << coord2 << "\n";
            }
        }
        else
        {
            if(line == "(")
                open_data = true;
            else
                continue;
        }
    }

    input_file.close();
    output_file.close();
}



std::string trim(std::string& str, std::vector<char>& lchars, std::vector<char>& rchars)
{
    str = ltrim(str, lchars);
    str = rtrim(str, rchars);

    return str;
}

std::string ltrim(std::string& str, std::vector<char>& chars)
{
    int begin, str_len = str.length();

    for(begin=0; begin<str_len; begin++)
        if(std::find(chars.begin(), chars.end(), str[begin]) == chars.end())
            break;

    return str.substr(begin, str_len-begin);
}

std::string rtrim(std::string& str, std::vector<char>& chars)
{
    int end, str_len = str.length();

    for(end=str_len-1; end>=0; end--)
        if(std::find(chars.begin(), chars.end(), str[end]) == chars.end())
            break;

    return str.substr(0, end+1);
}

void write_output_data(std::ofstream& output_file, FaceMesh2D& mesh, std::string& intermediate_filepath, std::vector<std::string>& fields)
{
    std::vector<Point> cell_centroids = mesh.get_cell_centroids();
    std::size_t num_cells = cell_centroids.size();

    std::size_t num_fields = fields.size();
    std::vector<std::ifstream> input_files(num_fields);
    for(std::size_t i=0; i<num_fields; i++)
    {
        input_files[i].open(intermediate_filepath + fields[i]);
    }

    std::vector<std::string> lines(num_fields);

    for(std::size_t cell_id = 0; cell_id < num_cells; cell_id++)
    {
        // initialisations
        std::array<double, N_DIMS> centroid_coords = cell_centroids[cell_id].get_coordinates();

        // write cell centroid
        for(Dimensions dim: ALL_DIMS)
            output_file<<centroid_coords[dim]<<"\t";

        // loop through fields and write
        for(std::size_t i=0; i<num_fields; i++)
        {
            std::getline(input_files[i], lines[i]);

            if(fields[i] == "U")
            {
                std::istringstream isstream(lines[i]);
                std::array<double, N_DIMS> velocity_coords;
                for(Dimensions dim: ALL_DIMS)
                    isstream >> velocity_coords[dim];

                double velocity_modulus = Point(velocity_coords).calc_modulus();

                // write velocity vector as (centroid-velocity/2, centroid+velocity/2);
                // scale reciprocally wrt modulus
                double grid_scaling = 5;
                double velocity_scaling = (1 - 1./(2*velocity_modulus + 1)) * grid_scaling;

                // beginning coordinates
                for(Dimensions dim: ALL_DIMS)
                    output_file<<centroid_coords[dim] - velocity_coords[dim]/2 / (velocity_modulus+EPSILON) * velocity_scaling<<"\t";

                // ending coordinates
                for(Dimensions dim: ALL_DIMS)
                    output_file<<velocity_coords[dim] / (velocity_modulus+EPSILON) * velocity_scaling<<"\t";

                // velocity components
                for(Dimensions dim: ALL_DIMS)
                    output_file << velocity_coords[dim] << "\t";

                // velocity magnitude
                output_file<<velocity_modulus<<"\t";
            }
            else
            {
                output_file<<lines[i]<<"\t";
            }
        }

        output_file<<std::endl;
    }
}

void write_comparison(std::ofstream& output_file, FaceMesh2D& mesh,
                        std::string& file, std::string& file_pimple)
{
    // assign columns in .dat documents
    std::unordered_map<std::string, std::size_t> columns;
    columns["X"] = 1;
    columns["Y"] = 2;
    columns["P"] = 3;
    columns["UX_START"] = 4;
    columns["UY_START"] = 5;
    columns["UX_LEN"] = 6;
    columns["UY_LEN"] = 7;
    columns["UX"] = 8;
    columns["UY"] = 9;
    columns["U_MOD"] = 10;
    columns["K"] = 11;
    columns["EPS"] = 12;
    columns["NU_TURB"] = 13;

    std::size_t num_columns = columns.size();

    std::vector<std::size_t> non_plotting_fields{columns["UX_START"], columns["UY_START"], columns["UX_LEN"], columns["UY_LEN"]};

    std::ifstream data_simulation(file), data_pimple(file_pimple);
    int num_cells = mesh.get_num_cells();

    std::string line, line_pimple;

    for(int cell_id = 0; cell_id < num_cells; cell_id++)
    {
        // set up istringstreams
        std::getline(data_simulation, line);
        std::getline(data_pimple, line_pimple);

        std::istringstream isstream(line), isstream_pimple(line_pimple);

        // set up velocity coordinates
        std::array<double, N_DIMS> velocity_coords, velocity_coords_pimple;

        // loop through fields and write
        for(std::size_t i=1; i<=num_columns; i++)
        {
            double value, value_pimple;
            isstream >> value;
            isstream_pimple >> value_pimple;

            // save value of velocity coordinate
            if(i == columns["UX"])
            {
                velocity_coords[X] = value;
                velocity_coords_pimple[X] = value_pimple;
            }
            else if(i == columns["UY"])
            {
                velocity_coords[Y] = value;
                velocity_coords_pimple[Y] = value_pimple;
            }

            // write data
            if(i == columns["X"] || i == columns["Y"])
            {
                output_file << value << "\t";   // simply write coordinate
            }
            else if(std::find(non_plotting_fields.begin(), non_plotting_fields.end(), i) != non_plotting_fields.end())
            {
                continue;
            }
            else if(i == columns["U_MOD"])
            {
                Point velocity(velocity_coords), velocity_pimple(velocity_coords_pimple);

                // write the difference in magnitude, taking pimple as the base
                output_file << Point::distance(velocity, velocity_pimple).calc_modulus()
                            << "\t";
            }
            else
            {
                // write the difference in magnitude, taking pimple as the base
                output_file << value - value_pimple
                            << "\t";
            }
        }

        output_file<<std::endl;
    }
}



void clear_dir_contents(std::string& dirpath)
{
    std::filesystem::path dir(dirpath);
    for (const auto& entry: std::filesystem::directory_iterator(dir)) 
        std::filesystem::remove_all(entry.path());
}