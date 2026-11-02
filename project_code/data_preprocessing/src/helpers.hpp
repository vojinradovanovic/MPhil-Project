#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "face_mesh.hpp"

void strip_scalar_field(std::ifstream& input_file, std::ofstream& output_file);
void strip_vector_field(std::ifstream& input_file, std::ofstream& output_file);

std::string trim(std::string& str, std::vector<char>& lchars, std::vector<char>& rchars);
std::string ltrim(std::string& str, std::vector<char>& chars);
std::string rtrim(std::string& str, std::vector<char>& chars);

void write_output_data(std::ofstream& output_file, FaceMesh2D& mesh, std::string& intermediate_filepath, std::vector<std::string>& fields);
void write_comparison(std::ofstream& output_file, FaceMesh2D& mesh,
                        std::string& file, std::string& file_pimple);

void clear_dir_contents(std::string& dirpath);