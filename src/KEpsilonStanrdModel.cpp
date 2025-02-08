#include "KEpsilonStanrdModel.hpp"
#include "Grid.hpp"
#include "Fields.hpp"
#include "Communication.hpp"
#include "controal_data.hpp"
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <math.h>
#include <omp.h>

void KEpsilonStanrdModel::calculate_k_and_epsilon(Fields &fields, Grid &grid) {
    auto K_OLD = fields.k_matrix();
    auto EPS_OLD = fields.eps_matrix();

    const Real sigma_k = 1.0;
    const Real sigma_epsilon = 1.2;
    const Real C2 = 1.9;

    for (const auto& current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();
        Real f2_coeff = 1;
        Real nut = fields.nu_t(i, j);
        Real kij = K_OLD(i, j);
        Real eij = EPS_OLD(i, j);

        auto k1_1 = Discretization::convection_uKEPS(fields.u_matrix(), K_OLD, i, j);
        auto k1_2 = Discretization::convection_vKEPS(fields.v_matrix(), K_OLD, i, j);
        auto e1_1 = Discretization::convection_uKEPS(fields.u_matrix(), EPS_OLD, i, j);
        auto e1_2 = Discretization::convection_vKEPS(fields.v_matrix(), EPS_OLD, i, j);

                auto k2 = Discretization::laplacian_nu(K_OLD, fields._nu, fields._NU_I, fields._NU_J, i, j);
                auto e2 = Discretization::laplacian_nu(EPS_OLD, fields._nu, fields._NU_I, fields._NU_J, i, j, 1.3);

               auto k3 = nut * Discretization::mean_strain_rate_squared(fields._U, fields._V,fields. _S, i, j);
        auto e3 = 1.44 * eij * k3 / (kij + 1e-10);
        auto e4 = f2_coeff * 1.92 * eij * eij / (kij + 1e-10);

        auto kij_new = kij + fields._dt * (-(k1_1 + k1_2) + k2 + k3 - eij);
        fields.k(i, j) = std::max(double(kij_new), 1e-8);

        auto epsij_new = eij + fields._dt * (-(e1_1 + e1_2) + e2 + e3 - e4);
        fields.eps(i, j) = std::max(double(epsij_new), 1e-8);

        assert(kij_new > 0);
        assert(epsij_new > 0);
    }
}

void KEpsilonStanrdModel::calculate_nu_t(Fields &fields, Grid &grid) {
    for (const auto &current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();
        Real fnu_coeff = 1;
        auto kij = fields.k(i, j);
        auto epsij = fields.eps(i, j);

        fields.nu_t(i, j) = fnu_coeff * 0.09 * kij * kij / epsij + fields._nu;

        assert(!isnan(fields.nu_t(i, j)));
        assert(!isinf(fields.nu_t(i, j)));
        assert(fields.nu_t(i, j) > 0);
    }

    for (const auto &current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();

        Real fnu_coeff = 1;
        auto num_i = (fields.k(i, j) + fields.k(i + 1, j)) / 2;
        auto denom_i = (fields.eps(i, j) + fields.eps(i + 1, j)) / 2;
        auto num_j = (fields.k(i, j) + fields.k(i, j + 1)) / 2;
        auto denom_j = (fields.eps(i, j) + fields.eps(i, j + 1)) / 2;
        fields.nu_i(i, j) = fnu_coeff * 0.09 * num_i * num_i / denom_i;
        fields.nu_j(i, j) = fnu_coeff * 0.09 * num_j * num_j / denom_j;

        assert(!isnan(fields.nu_i(i, j)));
        assert(!isnan(fields.nu_j(i, j)));
    }
    calculate_k_and_epsilon(fields, grid);
}
