#include "KEpsilonReModel.hpp"
#include "Grid.hpp"
#include "Fields.hpp"
#include "Communication.hpp"
#include"controal_data.hpp"
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <math.h>
#include <omp.h>


void KEpsilonReModel::calculate_k_and_epsilon(Fields &fields, Grid &grid) {
    auto K_OLD = fields.k_matrix();
    auto EPS_OLD = fields.eps_matrix();
    
    const Real sigma_k = 1.0;
    const Real sigma_epsilon = 1.2;
    const Real C2 = 1.9;

    for (const auto& current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();

        // 获取当前单元的变量
        Real kij = K_OLD(i, j);
        Real eij = EPS_OLD(i, j);
        Real nut = fields.nu_t(i, j);
        Real S_ij = fields._S(i, j); // 应该在模型计算之前已经计算好

        // 计算对流项
        auto k_conv = Discretization::convection_uKEPS(fields.u_matrix(), K_OLD, i, j) +
                      Discretization::convection_vKEPS(fields.v_matrix(), K_OLD, i, j);
        auto eps_conv = Discretization::convection_uKEPS(fields.u_matrix(), EPS_OLD, i, j) +
                        Discretization::convection_vKEPS(fields.v_matrix(), EPS_OLD, i, j);

        // 计算扩散项，使用 _nu 和 _NU_I, _NU_J
        auto k_diff = Discretization::laplacian_nu(K_OLD, fields._nu, fields._NU_I, fields._NU_J, i, j, sigma_k);
        auto eps_diff = Discretization::laplacian_nu(EPS_OLD, fields._nu, fields._NU_I, fields._NU_J, i, j, sigma_epsilon);

        // 计算湍流生成项 Pk
        Real Pk = nut * S_ij * S_ij;

        // 计算 eta
        Real eta = S_ij * kij / (eij + 1e-10);

        // 计算 C1
        Real C1 = std::max(0.43, eta / (eta + 5.0));

        // 计算 ε 的源项
        Real e3 = C1 * eij * Pk / (kij + 1e-10);
        Real e4 = C2 * eij * eij / (kij + sqrt(fields._nu * eij) + 1e-10);

        // 更新 k
        Real dkdt = -k_conv + k_diff + Pk - eij;
        Real kij_new = kij + fields._dt * dkdt;
        fields.k(i, j) = std::max(double(kij_new), 1e-8); // 防止负值

        // 更新 ε
        Real depsdt = -eps_conv + eps_diff + e3 - e4;
        Real eij_new = eij + fields._dt * depsdt;
        fields.eps(i, j) = std::max(double(eij_new), 1e-8); // 防止负值
    }
}
void KEpsilonReModel::calculate_nu_t(Fields &fields, Grid &grid) {
        for (const auto &current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();
        Real fnu_coeff = 1;
        auto kij = fields.k(i, j);
        auto epsij = fields.eps(i, j);

        // 设定 A0 和 As 的值（可调参数）
        Real A0 = 4.0;
        Real As = 0.9;

        // 防止除零错误
        Real denominator = epsij + 1e-10;
        Real numerator = kij / denominator;

        // 计算 C_mu
        Real C_mu = 1.0 / (A0 + As * fields._S(i, j) * numerator);

        // 确保 C_mu 为正值
        C_mu = std::max(double(C_mu), 1e-5);

        // 计算湍流粘性系数
        fields.nu_t(i, j) = C_mu * kij * kij / epsij + fields._nu;

        // 防止负值
        fields.nu_t(i, j) = std::max(fields.nu_t(i, j), fields._nu);

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
    calculate_k_and_epsilon(fields,grid);
}

