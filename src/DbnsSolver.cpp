// DbnsSolver.cpp

#include "DbnsSolver.hpp"
#include "Discretization.hpp"
//#include "Logger.hpp"

DbnsSolver::DbnsSolver() {
    //solver_type = SolverType::DBNS;
}

void DbnsSolver::initialize() {
    // ��ʼ�����е� Fields ��������Ҫ�ı���
    // ���磬���ó�ʼ�ܶȡ��������¶ȷֲ�
}



void DbnsSolver::solve_pressure(Real& res, uint32_t& iters) {
    // �����Ҫ����ѹ����⣬����������ʵ��
    // ��������ͨ��״̬����ֱ�Ӽ���ѹ�����˴��������ջ�ʵ�������߼�
    res = 0.0;
    iters = 0;
}



Real DbnsSolver::compute_flux_F(int i, int j) {
    // ����ͨ�� F = rho * u
    Grid& grid = _grid;
    Fields& field = _field;
    Real F = field.rho(i, j) * field.u(i, j);
    return F;
}

Real DbnsSolver::compute_flux_G(int i, int j) {
    // ����ͨ�� G = rho * v
    Grid& grid = _grid;
    Fields& field = _field;
    Real G = field.rho(i, j) * field.v(i, j);
    return G;
}

Real DbnsSolver::compute_flux_H(int i, int j) {
    // ����ͨ�� H = rho * E
    Grid& grid = _grid;
    Fields& field = _field;
    Real H = field.rho(i, j) * field.E(i, j);
    return H;
}


void DbnsSolver::solve_pre_pressure(Real& dt) {
    Grid& grid = _grid;
    Fields& field = _field;

    // ������ʱ�洢���ڸ���
    Matrix<Real> rho_new = field.rho_matrix();
    Matrix<Real> E_new = field.E_matrix();
    Matrix<Real> u_new = field.u_matrix(); // �����з�������ȡ��������
    Matrix<Real> v_new = field.v_matrix();

    // �����������嵥Ԫ������ʽ����
    for (const auto& current_cell : grid.fluid_cells()) {
        int i = current_cell->i();
        int j = current_cell->j();

        // �����Է��̣������غ㣩
        Real du_dx = (field.u(i + 1, j) - field.u(i - 1, j)) / (2 * grid.dx());
        Real dv_dy = (field.v(i, j + 1) - field.v(i, j - 1)) / (2 * grid.dy());
        Real continuity = du_dx + dv_dy;

        // ����ѹ��
        Real p = field.calculate_pressure(i, j);

        // ��������
        Real du_dt = -(field.u(i, j) * (field.u(i + 1, j) - field.u(i - 1, j)) / (2 * grid.dx()) +
            field.v(i, j) * (field.u(i, j + 1) - field.u(i, j - 1)) / (2 * grid.dy())) -
            (field.calculate_pressure(i + 1, j) - field.calculate_pressure(i - 1, j)) / (2 * grid.dx());
        Real dv_dt = -(field.u(i, j) * (field.v(i + 1, j) - field.v(i - 1, j)) / (2 * grid.dx()) +
            field.v(i, j) * (field.v(i, j + 1) - field.v(i, j - 1)) / (2 * grid.dy())) -
            (field.calculate_pressure(i, j + 1) - field.calculate_pressure(i, j - 1)) / (2 * grid.dy());

        // �����ٶ�
        u_new(i, j) = field.u(i, j) + dt * du_dt;
        v_new(i, j) = field.v(i, j) + dt * dv_dt;

        // ��������
        Real dE_dt = -(field.u(i, j) * (field.E(i + 1, j) - field.E(i - 1, j)) / (2 * grid.dx()) +
            field.v(i, j) * (field.E(i, j + 1) - field.E(i, j - 1)) / (2 * grid.dy())) +
            (field.calculate_pressure(i + 1, j) * field.u(i + 1, j) - field.calculate_pressure(i - 1, j) * field.u(i - 1, j)) / (2 * grid.dx()) +
            (field.calculate_pressure(i, j + 1) * field.v(i, j + 1) - field.calculate_pressure(i, j - 1) * field.v(i, j - 1)) / (2 * grid.dy());

        // ��������
        E_new(i, j) = field.E(i, j) + dt * dE_dt;

        // �����ܶ�
        rho_new(i, j) = field.rho(i, j) + dt * (-field.rho(i, j) * continuity);
    }

    // �������б���
    field.u_matrix() = u_new;
    field.v_matrix() = v_new;
    field.E_matrix() = E_new;
    field.rho_matrix() = rho_new;

    // Ӧ�ñ߽�����
    //set_boundary_conditions();
    for (const auto& boundary : _boundaries) {
        boundary->enforce_boundary_conditions(field);
    }
}

void DbnsSolver:: solve_post_pressure() { }

void DbnsSolver:: solve_piso(Real& residual, uint32_t& iterations) { }