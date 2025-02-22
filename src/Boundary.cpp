#include "Boundary.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <mpi.h>

/////////// Constructors ///////////
OutletBoundary::OutletBoundary(std::vector<Cell *> *cells) : Boundary(cells) { _type = 0; }

InletBoundary::InletBoundary(std::vector<Cell *> *cells) : Boundary(cells) { _type = 1; }

InletBoundary::InletBoundary(std::vector<Cell *> *cells, std::unordered_map<int, Real> inlet_U,
                             std::unordered_map<int, Real> inlet_V, std::unordered_map<int, Real> inlet_T,
                             std::unordered_map<int, Real> inlet_K, std::unordered_map<int, Real> inlet_EPS, Real DP)
    : Boundary(cells), _inlet_U(inlet_U), _inlet_V(inlet_V), _inlet_T(inlet_T), _inlet_K(inlet_K),
      _inlet_EPS(inlet_EPS) {
    _type = 1;
}

NoSlipWallBoundary::NoSlipWallBoundary(std::vector<Cell *> *cells) : Boundary(cells) { _type = 2; }

NoSlipWallBoundary::NoSlipWallBoundary(std::vector<Cell *> *cells, std::unordered_map<int, Real> wall_velocity,
                                       std::unordered_map<int, Real> wall_temperature)
    : Boundary(cells), _wall_velocity(wall_velocity) {
    _wall_temperature = wall_temperature;
    _type = 2;
}

FreeSlipWallBoundary::FreeSlipWallBoundary(std::vector<Cell *> *cells) : Boundary(cells) { _type = 3; }

FreeSlipWallBoundary::FreeSlipWallBoundary(std::vector<Cell *> *cells, std::unordered_map<int, Real> wall_velocity,
                                           std::unordered_map<int, Real> wall_temperature)
    : Boundary(cells), _wall_velocity(wall_velocity) {
    _wall_temperature = wall_temperature;
    _type = 3;
}

/////////// Main Boundary ///////////
void Boundary::enforce_fg(Fields &field) {
    for (auto &cell : *_cells) {
        int i = cell->i();
        int j = cell->j();
        if (cell->is_border(border_position::RIGHT)) {
            field.f(i, j) = field.u(i, j);
        }
        if (cell->is_border(border_position::LEFT)) {
            field.f(i - 1, j) = field.u(i - 1, j);
        }
        if (cell->is_border(border_position::TOP)) {
            field.g(i, j) = field.v(i, j);
        }
        if (cell->is_border(border_position::BOTTOM)) {
            field.g(i, j - 1) = field.v(i, j - 1);
        }
    }
}

void Boundary::enforce_p(Fields &field) {
    for (auto &cell : *_cells) {
        if (cell->borders().size() > 2) {
            std::cerr << "Forbidden cells!!" << std::endl;
            assert(false);
        } else if (cell->borders().size() == 2) {
            enforce_p_diagonal(field, cell);
        } else if (cell->borders().size() == 1) {
            enforce_p_main(field, cell);
        }
    }
}

void Boundary::enforce_p_main(Fields &field, Cell *cell) {
    int i = cell->i();
    int j = cell->j();
    if (cell->is_border(border_position::RIGHT)) {
        field.p(i, j) = field.p(i + 1, j);
    } else if (cell->is_border(border_position::LEFT)) {
        field.p(i, j) = field.p(i - 1, j);
    } else if (cell->is_border(border_position::TOP)) {
        field.p(i, j) = field.p(i, j + 1);
    } else if (cell->is_border(border_position::BOTTOM)) {
        field.p(i, j) = field.p(i, j - 1);
    }
}

void Boundary::enforce_p_diagonal(Fields &field, Cell *cell) {
    int i = cell->i();
    int j = cell->j();
    if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::TOP)) {
        field.p(i, j) = (field.p(i + 1, j) + field.p(i, j + 1)) / 2.0;
    }
    if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::BOTTOM)) {
        field.p(i, j) = (field.p(i + 1, j) + field.p(i, j - 1)) / 2.0;
    }
    if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::TOP)) {
        field.p(i, j) = (field.p(i - 1, j) + field.p(i, j + 1)) / 2.0;
    }
    if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::BOTTOM)) {
        field.p(i, j) = (field.p(i - 1, j) + field.p(i, j - 1)) / 2.0;
    }
}

void Boundary::enforce_t(Fields &field) {
    for (auto &cell : *_cells) {
        size_t num_borders = cell->borders().size();
        assert(num_borders <= 2);
        auto wt = _wall_temperature[cell->id()];
        if (num_borders == 1) {
            if (wt != -1) {
                enforce_t_drichlet_main(field, cell);
            } else {
                enforce_t_adiabatic_main(field, cell);
            }
        } else {
            if (wt != -1) {
                enforce_t_drichlet_diag(field, cell);
            } else {
                enforce_t_adiabatic_diag(field, cell);
            }
        }
    }
}

void Boundary::enforce_t_drichlet_main(Fields &field, Cell *cell) {
    int i = cell->i();
    int j = cell->j();
    auto wt = _wall_temperature[cell->id()];
    if (cell->is_border(border_position::RIGHT)) {
        field.t(i, j) = 2 * wt - field.t(i + 1, j);
    }
    if (cell->is_border(border_position::LEFT)) {
        field.t(i, j) = 2 * wt - field.t(i - 1, j);
    }
    if (cell->is_border(border_position::TOP)) {
        field.t(i, j) = 2 * wt - field.t(i, j + 1);
    }
    if (cell->is_border(border_position::BOTTOM)) {
        field.t(i, j) = 2 * wt - field.t(i, j - 1);
    }
}

void Boundary::enforce_t_drichlet_diag(Fields &field, Cell *cell) {
    int i = cell->i();
    int j = cell->j();
    auto wt = _wall_temperature[cell->id()];
    if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::TOP)) {
        field.t(i, j) = 2 * wt - (field.t(i + 1, j) + field.t(i, j + 1) / 2);
    } else if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::BOTTOM)) {
        field.t(i, j) = 2 * wt - (field.t(i + 1, j) + field.t(i, j - 1) / 2);
    }
    if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::TOP)) {
        field.t(i, j) = 2 * wt - (field.t(i - 1, j) + field.t(i, j + 1) / 2);
    }
    if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::BOTTOM)) {
        field.t(i, j) = 2 * wt - (field.t(i - 1, j) + field.t(i, j - 1) / 2);
    }
}

void Boundary::enforce_t_adiabatic_main(Fields &field, Cell *cell) {
    int i = cell->i();
    int j = cell->j();
    if (cell->is_border(border_position::RIGHT)) {
        field.t(i, j) = field.t(i + 1, j);
    }
    if (cell->is_border(border_position::LEFT)) {
        field.t(i, j) = field.t(i - 1, j);
    }
    if (cell->is_border(border_position::TOP)) {
        field.t(i, j) = field.t(i, j + 1);
    }
    if (cell->is_border(border_position::BOTTOM)) {
        field.t(i, j) = field.t(i, j - 1);
    }
}

void Boundary::enforce_t_adiabatic_diag(Fields &field, Cell *cell) {
    int i = cell->i();
    int j = cell->j();
    if (cell->is_border(border_position::RIGHT)) {
        field.t(i, j) = field.t(i + 1, j);
    }
    if (cell->is_border(border_position::LEFT)) {
        field.t(i, j) = field.t(i - 1, j);
    }
    if (cell->is_border(border_position::TOP)) {
        field.t(i, j) = field.t(i, j + 1);
    }
    if (cell->is_border(border_position::BOTTOM)) {
        field.t(i, j) = field.t(i, j - 1);
    }
}

void Boundary::enforce_nu_t(Fields &field, int turb_model) {
    for (auto &cell : *_cells) {
        int i = cell->i();
        int j = cell->j();
        if (cell->borders().size() > 2) {
            std::cerr << "Forbidden cells!!" << std::endl;
            assert(false);
        } else if (cell->borders().size() == 2) {
            if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::TOP)) {
                auto k_interp = (field.k(i + 1, j) + field.k(i, j + 1)) / 2.0;
                auto eps_interp = (field.eps(i + 1, j) + field.eps(i, j + 1)) / 2.0;
                field.k(i, j) = k_interp;
                field.eps(i, j) = eps_interp;
                if (turb_model == 1) {
                    field.nu_t(i, j) = 0.09 * k_interp * k_interp / eps_interp + field._nu;
                } else if (turb_model == 2 || turb_model == 3) {
                    field.nu_t(i, j) = k_interp / eps_interp + field._nu;
                } 
            }
            if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::BOTTOM)) {
                auto k_interp = (field.k(i + 1, j) + field.k(i, j - 1)) / 2.0;
                auto eps_interp = (field.eps(i + 1, j) + field.eps(i, j - 1)) / 2.0;
                field.k(i, j) = k_interp;
                field.eps(i, j) = eps_interp;
                if (turb_model == 1) {
                    field.nu_t(i, j) = 0.09 * k_interp * k_interp / eps_interp + field._nu;
                } else if (turb_model == 2 || turb_model == 3) {
                    field.nu_t(i, j) = k_interp / eps_interp + field._nu;
                }
            }
            if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::TOP)) {
                auto k_interp = (field.k(i - 1, j) + field.k(i, j + 1)) / 2.0;
                auto eps_interp = (field.eps(i - 1, j) + field.eps(i, j + 1)) / 2.0;
                field.k(i, j) = k_interp;
                field.eps(i, j) = eps_interp;
                if (turb_model == 1) {
                    field.nu_t(i, j) = 0.09 * k_interp * k_interp / eps_interp + field._nu;
                } else if (turb_model == 2 || turb_model == 3) {
                    field.nu_t(i, j) = k_interp / eps_interp + field._nu;
                }
            }
            if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::BOTTOM)) {
                auto k_interp = (field.k(i - 1, j) + field.k(i, j - 1)) / 2.0;
                auto eps_interp = (field.eps(i - 1, j) + field.eps(i, j - 1)) / 2.0;
                field.k(i, j) = k_interp;
                field.eps(i, j) = eps_interp;
                if (turb_model == 1) {
                    field.nu_t(i, j) = 0.09 * k_interp * k_interp / eps_interp + field._nu;
                } else if (turb_model == 2 || turb_model == 3) {
                    field.nu_t(i, j) = k_interp / eps_interp + field._nu;
                }
            }
        } else if (cell->borders().size() == 1) {
            if (cell->is_border(border_position::RIGHT)) {
                auto k = field.k(i + 1, j);
                auto eps = field.eps(i + 1, j);
                field.k(i, j) = k;
                field.eps(i, j) = eps;
                if (turb_model == 1) {
                    field.nu_t(i, j) = 0.09 * k * k / eps + field._nu;
                } else if (turb_model == 2 || turb_model == 3) {
                    field.nu_t(i, j) = k / eps + field._nu;
                } 
            }
            if (cell->is_border(border_position::LEFT)) {
                auto k = field.k(i - 1, j);
                auto eps = field.eps(i - 1, j);
                field.k(i, j) = k;
                field.eps(i, j) = eps;
                if (turb_model == 1) {
                    field.nu_t(i, j) = 0.09 * k * k / eps + field._nu;
                } else if (turb_model == 2 || turb_model == 3) {
                    field.nu_t(i, j) = k / eps + field._nu;
                } 
            }
            if (cell->is_border(border_position::TOP)) {
                auto k = field.k(i, j + 1);
                auto eps = field.eps(i, j + 1);
                field.k(i, j) = k;
                field.eps(i, j) = eps;
                if (turb_model == 1) {
                    field.nu_t(i, j) = 0.09 * k * k / eps + field._nu;
                } else if (turb_model == 2 || turb_model == 3) {
                    field.nu_t(i, j) = k / eps + field._nu;
                }
            }
            if (cell->is_border(border_position::BOTTOM)) {
                auto k = field.k(i, j - 1);
                auto eps = field.eps(i, j - 1);
                field.k(i, j) = k;
                field.eps(i, j) = eps;
                if (turb_model == 1) {
                    field.nu_t(i, j) = 0.09 * k * k / eps + field._nu;
                } else if (turb_model == 2 || turb_model == 3) {
                    field.nu_t(i, j) = k / eps + field._nu;
                } 
            }
        } else if (!cell->borders().size()) {
            field.nu_t(i, j) = field._nu;
        }
    }
}

/////////// Outlet ///////////
void OutletBoundary::enforce_p(Fields &field) {
    for (auto &cell : *_cells) {
        int i = cell->i();
        int j = cell->j();
        int id = cell->id();
        field.p(i, j) = field._PI;
    }
}

void OutletBoundary::enforce_t(Fields &field) {
    for (auto &cell : *_cells) {
        enforce_t_adiabatic_main(field, cell);
    }
}

void OutletBoundary::enforce_boundary_conditions(Fields& field){
    for (auto cell : *_cells) {
        int i = cell->i();
        int j = cell->j();

        // 设置压力
        field.p(i, j) = field._PI; // 示例值

        // 设置密度
        field.rho(i, j) = 1.0; // 示例值，根据物理条件调整

        // 设置温度
        field.t(i, j) = 300.0; // 示例值，根据物理条件调整

        // 根据状态方程计算能量
        Real R = 287.05; // 气体常数，干空气为例
        Real gamma = 1.4; // 比热比，干空气为例
        field.E(i, j) = (field.p(i, j) / (gamma - 1)) + 0.5 * field.rho(i, j) * (field.u(i, j) * field.u(i, j) + field.v(i, j) * field.v(i, j));
    }

}


/////////// Inlet ///////////
void InletBoundary::enforce_uv(Fields &field) {
    for (auto &cell : *_cells) {
        int i = cell->i();
        int j = cell->j();
        int id = cell->id();

        if (cell->is_border(border_position::RIGHT)) {
            field.u(i, j) = _inlet_U[id];
            field.v(i, j) = 2 * _inlet_V[id] - field.v(i + 1, j);
        } else if (cell->is_border(border_position::LEFT)) {
            field.u(i, j) = _inlet_U[id];
            field.v(i, j) = 2 * _inlet_V[id] - field.v(i - 1, j);
        } else if (cell->is_border(border_position::TOP)) {
            field.u(i, j) = 2 * _inlet_U[id] - field.u(i, j + 1);
            field.v(i, j) = _inlet_V[id];
        } else if (cell->is_border(border_position::BOTTOM)) {
            field.u(i, j) = 2 * _inlet_U[id] - field.u(i, j - 1);
            field.v(i, j) = _inlet_V[id];
        }
    }
}

void InletBoundary::enforce_p(Fields &field) {
    if (_inlet_DP == REAL_MAX) { // Inlet BC isn't specified, apply default
        Boundary::enforce_p(field);
        return;
    }
    for (auto &cell : *_cells) {
        int i = cell->i();
        int j = cell->j();
        int id = cell->id();
        field.p(i, j) = _inlet_DP + field._PI;
    }
}

void InletBoundary::enforce_t(Fields &field) {
    for (auto &cell : *_cells) {
        int i = cell->i();
        int j = cell->j();
        auto wt = _inlet_T[cell->id()];
        if (cell->is_border(border_position::RIGHT)) {
            field.t(i, j) = 2 * wt - field.t(i + 1, j);
        }
        if (cell->is_border(border_position::LEFT)) {
            field.t(i, j) = 2 * wt - field.t(i - 1, j);
        }
        if (cell->is_border(border_position::TOP)) {
            field.t(i, j) = 2 * wt - field.t(i, j + 1);
        }
        if (cell->is_border(border_position::BOTTOM)) {
            field.t(i, j) = 2 * wt - field.t(i, j - 1);
        }
    }
}

void InletBoundary::enforce_nu_t(Fields &field, int turb_model) {
    for (auto &cell : *_cells) {
        int i = cell->i();
        int j = cell->j();

        auto wk = _inlet_K[cell->id()];
        auto weps = _inlet_EPS[cell->id()];
        if (cell->is_border(border_position::RIGHT)) {
            auto k = 2 * wk - field.k(i + 1, j);
            auto eps = 2 * weps - field.eps(i + 1, j);
            field.k(i, j) = k;
            field.eps(i, j) = eps;
          /*  if (turb_model == 1) {
                field.nu_t(i, j) = 0.09 * field.k(i, j) * field.k(i, j) / field.eps(i, j) + field._nu;
            } else if (turb_model == 2) {
                field.nu_t(i, j) = k / eps + field._nu;
            }*/
        }
        if (cell->is_border(border_position::LEFT)) {
            auto k = 2 * wk - field.k(i - 1, j);
            auto eps = 2 * weps - field.eps(i - 1, j);
            field.k(i, j) = k;
            field.eps(i, j) = eps;
           /* if (turb_model == 1) {
                field.nu_t(i, j) = 0.09 * field.k(i, j) * field.k(i, j) / field.eps(i, j) + field._nu;
            } else if (turb_model == 2) {
                field.nu_t(i, j) = k / eps + field._nu;
            }*/
        }
        if (cell->is_border(border_position::TOP)) {
            auto k = 2 * wk - field.k(i, j + 1);
            auto eps = 2 * weps - field.eps(i, j + 1);
            field.k(i, j) = k;
            field.eps(i, j) = eps;
           /* if (turb_model == 1) {
                field.nu_t(i, j) = 0.09 * field.k(i, j) * field.k(i, j) / field.eps(i, j) + field._nu;
            } else if (turb_model == 2) {
                field.nu_t(i, j) = k / eps + field._nu;
            }*/
        }
        if (cell->is_border(border_position::BOTTOM)) {
            auto k = 2 * wk - field.k(i, j - 1);
            auto eps = 2 * weps - field.eps(i, j - 1);
            field.k(i, j) = k;
            field.eps(i, j) = eps;
           /* if (turb_model == 1) {
                field.nu_t(i, j) = 0.09 * field.k(i, j) * field.k(i, j) / field.eps(i, j) + field._nu;
            } else if (turb_model == 2) {
                field.nu_t(i, j) = k / eps + field._nu;
            }*/
        }
    }
}

void InletBoundary::enforce_boundary_conditions(Fields& field) {
    Real gamma = 1.4;
    for (auto& cell : *_cells) {
        int i = cell->i();
        int j = cell->j();
        int id = cell->id();
        auto wt = _inlet_T[cell->id()];
        //for p
        field.p(i, j) = _inlet_DP + field._PI;
        //for rho
        field.rho(i, j) = _inlet_rho;
        //for e
        field.E(i, j) = (field.p(i, j) / (gamma - 1)) + 0.5 * field.rho(i, j) * (field.u(i, j) * field.u(i, j) + field.v(i, j) * field.v(i, j));
        //for uv and T 
        if (cell->is_border(border_position::RIGHT)) {
            field.u(i, j) = _inlet_U[id];
            field.v(i, j) = 2 * _inlet_V[id] - field.v(i + 1, j);
            field.t(i, j) = 2 * wt - field.t(i + 1, j);
        }
        else if (cell->is_border(border_position::LEFT)) {
            field.u(i, j) = _inlet_U[id];
            field.v(i, j) = 2 * _inlet_V[id] - field.v(i - 1, j);
            field.t(i, j) = 2 * wt - field.t(i - 1, j);
        }
        else if (cell->is_border(border_position::TOP)) {
            field.u(i, j) = 2 * _inlet_U[id] - field.u(i, j + 1);
            field.v(i, j) = _inlet_V[id];
            field.t(i, j) = 2 * wt - field.t(i, j + 1);
        }
        else if (cell->is_border(border_position::BOTTOM)) {
            field.u(i, j) = 2 * _inlet_U[id] - field.u(i, j - 1);
            field.v(i, j) = _inlet_V[id];
            field.t(i, j) = 2 * wt - field.t(i, j - 1);
        }
    }
}

/////////// NoSlip Walls ///////////
void NoSlipWallBoundary::enforce_uv(Fields &field) {

    for (auto &cell : *_cells) {
        if (cell->borders().size() > 2) {
            std::cerr << "Forbidden cells!!" << std::endl;
            assert(false);
        } else if (cell->borders().size() == 2) {
            enforce_uv_diagonal(field, cell);
        } else if (cell->borders().size() == 1) {
            enforce_uv_main(field, cell);
        }
    }
}

void NoSlipWallBoundary::enforce_uv_main(Fields &field, Cell *cell) {
    int i = cell->i();
    int j = cell->j();
    int id = cell->id();

    if (cell->is_border(border_position::RIGHT)) {
        field.u(i, j) = 0;
        field.v(i, j) = 2 * _wall_velocity[id] - field.v(i + 1, j);
    } else if (cell->is_border(border_position::LEFT)) {
        field.u(i - 1, j) = 0;
        field.v(i, j) = 2 * _wall_velocity[id] - field.v(i - 1, j);
    } else if (cell->is_border(border_position::TOP)) {
        field.v(i, j) = 0;
        field.u(i, j) = 2 * _wall_velocity[id] - field.u(i, j + 1);
    } else if (cell->is_border(border_position::BOTTOM)) {
        field.v(i, j - 1) = 0;
        field.u(i, j) = 2 * _wall_velocity[id] - field.u(i, j - 1);
    }
}

void NoSlipWallBoundary::enforce_uv_diagonal(Fields &field, Cell *cell) {
    int i = cell->i();
    int j = cell->j();
    int id = cell->id();
    // This shouldn't happen since index 0 is a ghost cell
    assert(!(i == 0 && j == 0) && "Illegal index");

    if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::TOP)) {
        field.u(i, j) = 0;
        field.u(i - 1, j) = 2 * _wall_velocity[id] - field.u(i - 1, j + 1);
        field.v(i, j) = 0;
        field.v(i, j - 1) = 2 * _wall_velocity[id] - field.v(i + 1, j - 1);
    } else if (cell->is_border(border_position::RIGHT) && cell->is_border(border_position::BOTTOM)) {
        field.u(i, j) = 0;
        field.u(i - 1, j) = 2 * _wall_velocity[id] - field.u(i - 1, j - 1);
        field.v(i, j) = 2 * _wall_velocity[id] - field.v(i + 1, j);
        field.v(i, j - 1) = 0;
    } else if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::TOP)) {
        field.u(i, j) = 2 * _wall_velocity[id] - field.u(i, j + 1);
        field.u(i - 1, j) = 0;
        field.v(i, j) = 0;
        field.v(i, j - 1) = 2 * _wall_velocity[id] - field.v(i - 1, j - 1);
    } else if (cell->is_border(border_position::LEFT) && cell->is_border(border_position::BOTTOM)) {
        field.u(i, j) = 2 * _wall_velocity[id] - field.u(i, j - 1);
        field.u(i - 1, j) = 0;
        field.v(i, j) = 2 * _wall_velocity[id] - field.v(i - 1, j);
        field.v(i, j - 1) = 0;
    }
}

void NoSlipWallBoundary::enforce_boundary_conditions(Fields& field) {
    //暂时不考虑该边界条件
}

/////////// FreeSlip Walls ///////////
void FreeSlipWallBoundary::enforce_uv(Fields &field) {
    for (auto &cell : *_cells) {
        int i = cell->i();
        int j = cell->j();
        int id = cell->id();
        assert(cell->borders().size() <= 2);
        if (cell->is_border(border_position::RIGHT)) {
            field.u(i, j) = 0;
        }
        if (cell->is_border(border_position::LEFT)) {
            field.u(i - 1, j) = 0;
        }
        if (cell->is_border(border_position::TOP)) {
            field.v(i, j) = 0;
        }
        if (cell->is_border(border_position::BOTTOM)) {
            field.v(i, j - 1) = 0;
        }
    }
}

void FreeSlipWallBoundary::enforce_boundary_conditions(Fields& field) {
    //暂时不考虑该边界条件
}