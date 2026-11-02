#include <limits>
#include <string>

// constants
const double EPSILON = std::numeric_limits<double>::epsilon();

// declare functions to be used later
void write_mesh();

void write_cartesian_mesh();
void write_points_cartesian(int n_x, int n_y, double x_min, double x_max, double y_min, double y_max);
void write_faces_cartesian(int n_x, int n_y);
void write_cells_cartesian(int n_x, int n_y);
void write_patches_bcs_cartesian(int n_x, int n_y, int n_patches);

void write_scalar_field_cartesian(std::string& filepath, std::string& filename,
                        int n_x, int n_y,
                        std::string& init_field);
void write_vector_field_cartesian(std::string& filepath, std::string& filename,
                        int n_x, int n_y, double x_min, double x_max, double y_min, double y_max,
                        std::string& init_field);

void write_bw_step_mesh();
void write_points_bw_step(int n_x_inlet, int n_y_inlet,
                            double x_min_inlet, double x_max_inlet, double y_min_inlet, double y_max_inlet,
                            int n_x_outlet, int n_y_outlet,
                            double x_min_outlet, double x_max_outlet, double y_min_outlet, double y_max_outlet);
void write_faces_bw_step(int n_x_inlet, int n_y_inlet, int n_x_outlet, int n_y_outlet);
void write_cells_bw_step(int n_x_inlet, int n_y_inlet, int n_x_outlet, int n_y_outlet);
void write_patches_bcs_bw_step(int n_x_inlet, int n_y_inlet, int n_x_outlet, int n_y_outlet, int n_patches);

void write_scalar_field_bw_step(std::string& filepath, std::string& filename,
                        int n_x_inlet, int n_y_inlet, int n_x_outlet, int n_y_outlet,
                        std::string& init_field);
void write_vector_field_bw_step(std::string& filepath, std::string& filename,
                        int n_x_inlet, int n_y_inlet, int n_x_outlet, int n_y_outlet,
                        std::string& init_field);


std::string read_patch_bc();