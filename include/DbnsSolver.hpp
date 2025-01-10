// DbnsSolver.hpp
#ifndef DBNSSOLVER_HPP
#define DBNSSOLVER_HPP

#include "Solver.hpp"
#include "Fields.hpp"
#include "Grid.hpp"
//#include "Matrix.hpp"
#include <vector>
#include <memory>

class DbnsSolver : public Solver {
public:
    DbnsSolver();
    virtual ~DbnsSolver() = default;

    void initialize() override;
    //void solve(Real dt); 名字更换为solve_pre_pressure
    //void set_boundary_conditions();

    void solve_pressure(Real& res, uint32_t& iters)override;
    //void enforce_rho_boundary();
    //void enforce_E_boundary();
    //void enforce_momentum_boundary();

    // 计算通量
    Real compute_flux_F(int i, int j);
    Real compute_flux_G(int i, int j);
    Real compute_flux_H(int i, int j);

    void solve_pre_pressure(Real& dt)override;

    void solve_post_pressure()override;

    void solve_piso(Real& residual, uint32_t& iterations)override;
};

#endif // DBNSSOLVER_HPP
