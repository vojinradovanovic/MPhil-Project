#include <cassert>
#include <fstream>
#include <functional>
#include <iostream>
#include <unordered_map>

#include "pvcoupling2D.hpp"
#include "tensor.hpp"
#include "turbulence_constants.hpp"

#ifndef TURBULENCE_2D
#define TURBULENCE_2D

class Turbulence2D: public PVCoupling2D
{
    protected:
        // turbulence kinetic energy system and fields
        LinearSystem k_system;
        std::vector<double> k_field;        // turbulence kinetic energy
        std::vector<BoundaryConditions> k_bcs;
        std::vector<std::vector<double>> k_bc_values;

        // turbulent dissipation system and fields
        LinearSystem epsilon_system;
        std::vector<double> epsilon_field;
        std::vector<BoundaryConditions> epsilon_bcs;
        std::vector<std::vector<double>> epsilon_bc_values;

        // turbulent viscosity fields
        std::vector<double> nu_turb_field;
        std::vector<BoundaryConditions> nu_turb_bcs;
        std::vector<std::vector<double>> nu_turb_bc_values;

        // map containing the averaging weights of each wall adjacent cell
        std::unordered_map<int, double> wall_cell_weights;

        // flow constants
        double C1;
        double C2;
        double C_nu;
        double sigma_k;
        double sigma_epsilon;

        double kappa;   // von Karman constant
        double E;       // empirical wall roughness

        double inertial_sublayer_min;
        double inertial_sublayer_max;

        double k_underrel_coeff;
        double epsilon_underrel_coeff;

        InterpolationMethods k_intpol_method;

    public:
        // turbulence2D_inits.cpp
        Turbulence2D();
        void initialise_params() override;

        void set_conds() override;

        // turbulence2D_solvers.cpp
        void solve_k_epsilon();

    protected:
        // turbulence2D_helpers.cpp
        void calc_G_field_incompressible(std::vector<double>& G_field);
        void calc_sym_velocity_grad_field(std::vector<Tensor<Point>>& sym_velocity_grad_field,
                                            const std::vector<Tensor<Point>>& grad_vel_t_field);
        void calc_grad_velocity_transpose_field(std::vector<Tensor<Point>>& grad_vel_t_field);

        void calc_nu_turb_field();
        void calc_scaled_nu_eff_fields(std::vector<double>& scaled_nu_eff_field,
                                        std::vector<BoundaryConditions>& scaled_nu_eff_bcs,
                                        std::vector<std::vector<double>>& scaled_nu_eff_bc_values, double scaling_coeff);

        void update_k_field_vars(std::function<void ()>& k_update_func);
        void update_epsilon_field_vars(std::function<void ()>& epsilon_update_func);

        void calc_k_bcs();
        void calc_epsilon_bcs();
        bool has_zero_dirichlet_velocity(int patch_id);
        double calc_field_average(std::vector<double>& field);

        void enforce_min_field_value(std::vector<double>& field, const double min_value, const bool raise_warning = true);

        void apply_wall_function(const std::function<void (int, std::size_t, int, double, double)>& wall_func,
                                    const std::function<void ()>& preset_func = nullptr);

        void update_residuals(double& init_residual, std::array<double, RESIDUAL_NO>& residuals, LinearSystem& system, std::vector<double>& field);
        bool check_outer_convergence(double& init_residual, std::array<double, RESIDUAL_NO>& residuals);

        // turbulence2D_inits.cpp
        void calc_wall_cell_weights();

        // turbulence2D_rw.cpp
        void write_output(std::ofstream& output_file) override;

        // turbulence2D_solvers.cpp
        void solve_k(double rel_tol = INNER_REL_TOL);
        void solve_epsilon(double rel_tol = INNER_REL_TOL);

        // turbulence2D_sys_assembly.cpp
        void assemble_k_system(std::vector<double>& k_field_old, std::vector<double>& G_field);
        void assemble_epsilon_system(std::vector<double>& epsilon_field_old, std::vector<double>& G_field);

        void assemble_velocity_systems(std::vector<Point>& velocity_field_old);
        void assemble_pressure_system(const std::vector<double>& diffusivity_field,
                                        const std::vector<BoundaryConditions>& diffusivity_field_bcs,
                                        const std::vector<double>& diffusivity_field_bc_values);

        void calc_k_rhs_coeffs(std::vector<double>& G_field);
        void calc_epsilon_rhs_coeffs(std::vector<double>& G_field);
        void calc_Sp_term_coeffs(LinearSystem& system, std::vector<double>& Sp_field);  // no eq_sign, always =1
};

#endif