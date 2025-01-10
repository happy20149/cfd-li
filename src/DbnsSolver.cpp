// DbnsSolver.cpp

#include "DbnsSolver.hpp"
#include "Discretization.hpp"
//#include "Logger.hpp"

DbnsSolver::DbnsSolver() {
    //solver_type = SolverType::DBNS;
}

void DbnsSolver::initialize() {
    // 初始化现有的 Fields 和其他必要的变量
    // 例如，设置初始密度、能量和温度分布
}



void DbnsSolver::solve_pressure(Real& res, uint32_t& iters) {
    // 如果需要保留压力求解，可以在这里实现
    // 由于我们通过状态方程直接计算压力，此处可以留空或实现其他逻辑
    res = 0.0;
    iters = 0;
}



Real DbnsSolver::compute_flux_F(int i, int j) {
    // 计算通量 F = rho * u
    Grid& grid = _grid;
    Fields& field = _field;
    Real F = field.rho(i, j) * field.u(i, j);
    return F;
}

Real DbnsSolver::compute_flux_G(int i, int j) {
    // 计算通量 G = rho * v
    Grid& grid = _grid;
    Fields& field = _field;
    Real G = field.rho(i, j) * field.v(i, j);
    return G;
}

Real DbnsSolver::compute_flux_H(int i, int j) {
    // 计算通量 H = rho * E
    Grid& grid = _grid;
    Fields& field = _field;
    Real H = field.rho(i, j) * field.E(i, j);
    return H;
}


void DbnsSolver::solve_pre_pressure(Real& dt) {
    Grid& grid = _grid;
    Fields& field = _field;

    // 创建临时存储用于更新
    Matrix<Real> rho_new = field.rho_matrix();
    Matrix<Real> E_new = field.E_matrix();
    Matrix<Real> u_new = field.u_matrix(); // 假设有访问器获取完整矩阵
    Matrix<Real> v_new = field.v_matrix();

    // 遍历所有流体单元进行显式更新
    for (const auto& current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();

        // 连续性方程（质量守恒）
        Real du_dx = (field.u(i + 1, j) - field.u(i - 1, j)) / (2 * grid.dx());
        Real dv_dy = (field.v(i, j + 1) - field.v(i, j - 1)) / (2 * grid.dy());
        Real continuity = du_dx + dv_dy;

        // 计算压力
        Real p = field.calculate_pressure(i, j);

        // 动量方程
        Real du_dt = -(field.u(i, j) * (field.u(i + 1, j) - field.u(i - 1, j)) / (2 * grid.dx()) +
            field.v(i, j) * (field.u(i, j + 1) - field.u(i, j - 1)) / (2 * grid.dy())) -
            (field.calculate_pressure(i + 1, j) - field.calculate_pressure(i - 1, j)) / (2 * grid.dx());
        Real dv_dt = -(field.u(i, j) * (field.v(i + 1, j) - field.v(i - 1, j)) / (2 * grid.dx()) +
            field.v(i, j) * (field.v(i, j + 1) - field.v(i, j - 1)) / (2 * grid.dy())) -
            (field.calculate_pressure(i, j + 1) - field.calculate_pressure(i, j - 1)) / (2 * grid.dy());

        // 更新速度
        u_new(i, j) = field.u(i, j) + dt * du_dt;
        v_new(i, j) = field.v(i, j) + dt * dv_dt;

        // 能量方程
        Real dE_dt = -(field.u(i, j) * (field.E(i + 1, j) - field.E(i - 1, j)) / (2 * grid.dx()) +
            field.v(i, j) * (field.E(i, j + 1) - field.E(i, j - 1)) / (2 * grid.dy())) +
            (field.calculate_pressure(i + 1, j) * field.u(i + 1, j) - field.calculate_pressure(i - 1, j) * field.u(i - 1, j)) / (2 * grid.dx()) +
            (field.calculate_pressure(i, j + 1) * field.v(i, j + 1) - field.calculate_pressure(i, j - 1) * field.v(i, j - 1)) / (2 * grid.dy());

        // 更新能量
        E_new(i, j) = field.E(i, j) + dt * dE_dt;

        // 更新密度
        rho_new(i, j) = field.rho(i, j) + dt * (-field.rho(i, j) * continuity);
    }

    // 更新所有变量
    field.u_matrix() = u_new;
    field.v_matrix() = v_new;
    field.E_matrix() = E_new;
    field.rho_matrix() = rho_new;

    // 应用边界条件
    //set_boundary_conditions();
    for (const auto& boundary : _boundaries) {
        boundary->enforce_boundary_conditions(field);
    }
}

void DbnsSolver:: solve_post_pressure() { }

void DbnsSolver:: solve_piso(Real& residual, uint32_t& iterations) { }