#include <cassert>
#include <fstream>
#include <functional>
#include <iostream>

#include "compressed_row_matrix.hpp"
#include "face_mesh.hpp"
#include "linear_system.hpp"
#include "physics_constants.hpp"

#ifndef PV_COUPLING_2D
#define PV_COUPLING_2D

class PVCoupling2D: public FaceMesh2D
{
    protected:
        // pressure system and fields
        LinearSystem pressure_system;
        std::vector<double> pressure_field;
        std::vector<Point> grad_p_field;
        std::vector<BoundaryConditions> pressure_bcs;
        std::vector<std::vector<double>> pressure_bc_values;
        std::vector<BoundaryConditions> diffusivity_field_p_bcs;
        std::vector<std::vector<double>> diffusivity_field_p_bc_values;

        // one system for each velocity dimension
        std::array<LinearSystem, N_DIMS> velocity_systems;

        // velocity fields are memorised as vectors of Points
        std::vector<Point> velocity_field;
        std::vector<BoundaryConditions> velocity_bcs;
        std::vector<std::vector<Point>> velocity_bc_values;

        std::vector<double> fluxes;

        double Re;      // Reynolds number
        double rho;     // density, single constant for the entire mesh
        double char_length;
        double nu;
        std::vector<double> nu_field;
        std::vector<BoundaryConditions> nu_bcs;
        std::vector<std::vector<double>> nu_bc_values;

        double velocity_underrel_coeff;
        double pressure_underrel_coeff;

        InterpolationMethods convection_intpol_method;
        InterpolationMethods pressure_intpol_method;
        InterpolationMethods velocity_intpol_method;
        InterpolationMethods diffusivity_intpol_method;     // applicable to all diffusivity fields

        // other
        int num_internal_faces;

    public:
        // pvCoupling2D_inits.cpp
        PVCoupling2D();
        virtual void initialise_params();

        virtual void set_conds();
        
        // pvCoupling2D_solvers.cpp
        void solve_simple();
        void solve_diffusion();
        void solve_convection();
        void solve_advection_diffusion();

    protected:
        // pvCoupling2D_helpers.cpp
        void calc_char_length();
        void calc_nu_field();
        void calc_fluxes(InterpolationMethods intpol_method);
        void check_fluxes(bool check_cells);

        void update_fluxes_with_pressure(const std::vector<double>& diffusivity_field,
                                            const std::vector<BoundaryConditions>& diffusivity_field_bcs,
                                            const std::vector<std::vector<double>>& diffusivity_field_bc_values);
        void update_velocity_field_with_grad_p(std::vector<double>& diffusivity_field);

        void extract_dim_field(std::vector<Point>& field, std::vector<double>& dim_field, Dimensions dim);
        void extract_dim_field(std::vector<std::vector<Point>>& field, std::vector<std::vector<double>>& dim_field, Dimensions dim);
        void calc_pseudovelocity_field(const std::vector<double>& diffusivity_field);

        void calc_kinematic_pressure_field();
        std::vector<double> calc_pressure_diffusivity_field();
        Point calc_face_grad_p(int face_id);
        void calc_grad_field_of_scalar(const std::vector<double>& scalar_field, std::vector<Point>& grad_field,
                                        const std::vector<BoundaryConditions>& field_bcs,
                                        const std::vector<std::vector<double>>& field_bc_values, InterpolationMethods field_intpol_method);

        void update_velocity_field_vars(std::function<void ()>& velocity_update_func);
        void update_pressure_field_vars(std::function<void ()>& pressure_update_func);

        double interpolate_to_face(const std::vector<double>& field, const std::vector<BoundaryConditions>& field_bcs,
                                            const std::vector<std::vector<double>>& field_bc_values,
                                            int face_id, InterpolationMethods intpol_method);
        Point interpolate_to_face(const std::vector<Point>& field, const std::vector<BoundaryConditions>& field_bcs,
                                            const std::vector<std::vector<Point>>& field_bc_values,
                                            int face_id, InterpolationMethods intpol_method);


        void calc_patch_and_position(int face_id, int& patch_id, int& position_in_patch);

        // pvCoupling2D_rw.cpp
        void read_scalar_field(std::ifstream& input_file, std::vector<double>& scalar_field);
        void read_vector_field(std::ifstream& input_file, std::vector<Point>& vector_field);
        void read_boundary_conds(std::ifstream& boundary_conds_file);

        virtual void write_output(std::ofstream& output_file);
        void write_linear_systems_data(std::ofstream& output_file_matrices,
                                        std::ofstream& output_file_vectors,
                                        std::vector<LinearSystem>& linear_systems,
                                        std::vector<std::vector<double>> residuals_before,
                                        std::vector<std::vector<double>> residuals_after);

        // pvCoupling2D_solvers.cpp
        void solve_velocity(double rel_tol = INNER_REL_TOL);  // without printing
        void solve_velocity(std::ofstream& output_file_matrices,
                                std::ofstream& output_file_vectors,
                                double rel_tol = INNER_REL_TOL);    // with printing

        void solve_pressure(double rel_tol = INNER_REL_TOL);  // without printing
        void solve_pressure(std::ofstream& output_file_matrices,
                                std::ofstream& output_file_vectors,
                                double rel_tol = INNER_REL_TOL);    // with printing

        void query_cell(LinearSystem& lin_system, int cell_id);

        void solve_velocity_x_diffusion(double rel_tol = INNER_REL_TOL);

        // pvCoupling2D_sys_assembly.cpp
        void underrelax_system_implicit(LinearSystem& system, std::vector<double> scalar_field, double underrel_coeff);
        void underrelax_field_explicit(std::vector<double>& scalar_field, std::vector<double>& scalar_field_old,
                                        double underrel_coeff);

        void calc_time_der_coeffs(LinearSystem& system, int eq_sign, const std::vector<double>& old_field);

        void calc_convection_coeffs(LinearSystem& system, int eq_sign, const std::vector<BoundaryConditions>& field_bcs,
                                    const std::vector<std::vector<double>>& field_bc_values, InterpolationMethods intpol_method);
        void calc_convection_coeffs_upwind(LinearSystem& system, int eq_sign, const std::vector<BoundaryConditions>& field_bcs,
                                            const std::vector<std::vector<double>>& field_bc_values);

        void calc_diffusion_coeffs(LinearSystem& system, int eq_sign,
                                    const std::vector<BoundaryConditions>& field_bcs, const std::vector<std::vector<double>>& field_bc_values,
                                    const std::vector<double>& diffusivity_field,
                                    const std::vector<BoundaryConditions>& diffusivity_field_bcs,
                                    const std::vector<std::vector<double>>& diffusivity_field_bc_values,
                                    InterpolationMethods diffusivity_field_intpol_method);

        void calc_collocated_source_coeffs(LinearSystem& system, int eq_sign, const std::vector<double>& source_field);

        void assemble_velocity_systems();
        void assemble_pressure_system(const std::vector<double>& diffusivity_field,
                                        const std::vector<BoundaryConditions>& diffusivity_field_bcs,
                                        const std::vector<std::vector<double>>& diffusivity_field_bc_values);

        void calc_flux_source_coeffs(LinearSystem& system, int eq_sign, const std::vector<double>& flux_field);
        void calc_divergence_source_coeffs(LinearSystem& system, int eq_sign, const std::vector<Point>& source_field,
                                            const std::vector<BoundaryConditions>& source_field_bcs,
                                            const std::vector<std::vector<Point>>& source_field_bc_values,
                                            InterpolationMethods& source_field_intpol_method);

        void assemble_matrix(LinearSystem& system, int num_of_operations, std::function<void ()>& raw_assembly_function,
                                const std::function<void ()>& off_diagonal_postprocessing = nullptr);
        void assemble_matrix(std::array<LinearSystem, N_DIMS>& systems,
                                int num_of_operations, std::function<void ()>& raw_assembly_function,
                                const std::function<void ()>& off_diagonal_postprocessing = nullptr);
};

#endif