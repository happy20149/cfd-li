#pragma once

#include "Datastructures.hpp"
#include "Discretization.hpp"
#include "Grid.hpp"
#include "TurbulenceModel.hpp"


/**
 * @brief Class of container and modifier for the physical fields
 *
 */
class Fields {
  public:
    Fields() = default;

    /**
     * @brief Constructor for the fields
     *
     * @param[in] _nu kinematic viscosity
     * @param[in] _dt initial timestep size
     * @param[in] _tau adaptive timestep coefficient
     * @param[in] imax number of cells in x direction
     * @param[in] jmax number of cells in y direction
     * @param[in] UI initial x-velocity
     * @param[in] VI initial y-velocity
     * @param[in] PI initial pressure
     * @param[in] TI initial temperature
     * @param[in] KI initial k value
     * @param[in] EPSI initial epsilon value
     * @param[in] _alpha alpha parameter
     * @param[in] _beta beta parameter
     * @param[in] _gx external force in x direction
     * @param[in] _gy exteral force in y direction
     *
     */
    Fields(Real _nu, Real _dt, Real _tau, int imax, int jmax, Real UI, Real VI, Real PI, Real TI, Real KI, Real EPSI, Real _alpha,
           Real _beta, Real _gx, Real _gy);

    ~Fields() = default;

    /**
     * @brief Calculates the convective and diffusive fluxes in x and y
     * direction based on explicit discretization of the momentum equations
     *
     * @param[in] grid in which the fluxes are calculated
     * @param[in] calc_temp whether to include temperature related terms
     * @param[in] turbulent whether to inclue turbulence calculations
     *
     */
    void calculate_fluxes(Grid &grid, bool calc_temp, bool turbulent);

    /**
     * @brief Right hand side calculations using the fluxes for the pressure
     * Poisson equation
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_rs(Grid &grid);

    /**
     * @brief Velocity calculation using pressure values
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_velocities(Grid &grid);

    /**
     * @brief Velocity calculation using pressure values
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_temperatures(Grid &grid);

    /**
     * @brief Adaptive step size calculation using x-velocity condition,
     * y-velocity condition and CFL condition
     *
     * @param[in] grid in which the calculations are done
     * @param[in] calc_temp whether to include temperatures
     * @param[in] turbulence whether to include turbulence
     */
    Real calculate_dt(Grid &grid, bool calc_temp, int turbulence);

    /// x-velocity index based access and modify
    Real &u(int i, int j);

    /// y-velocity index based access and modify
    Real &v(int i, int j);

    /// pressre index based access and modify
    Real &p(int i, int j);

    /// temperature index based acccess and modify
    Real &t(int i, int j);

    /// RHS index based access and modify
    Real &rs(int i, int j);

    /// x-momentum flux index based access and modify
    Real &f(int i, int j);

    /// y-momentum flux index based access and modify
    Real &g(int i, int j);

    /// turbulent kinetic energy index based access and modify
    Real &k(int i, int j);

    /// dissipation rate index based access and modify
    Real &eps(int i, int j);

    /// turbulent viscosity index based access and modify
    Real &nu_t(int i, int j);
    /// turbulent viscosity on vertical edges index based access and modify
    Real &nu_i(int i, int j);
    /// turbulent viscosity on horizontal edges index based access and modify
    Real &nu_j(int i, int j);

    /// get timestep size
    Real dt() const;

    /// pressure matrix access and modify
    Matrix<Real> &p_matrix();

    /// velocity u matrix access and modify
    Matrix<Real> &u_matrix();

    /// velocity v matrix access and modify
    Matrix<Real> &v_matrix();

    /// temperature t matrix access and modify
    Matrix<Real> &t_matrix();

    /// x-momentum flux matrix
    Matrix<Real> &f_matrix();

    /// y-momentum flux matrix
    Matrix<Real> &g_matrix();

    /// turbulent kinetic energy matrix
    Matrix<Real> &k_matrix();

    /// dissipation rate matrix
    Matrix<Real> &eps_matrix();

    /// turbulent energy matrix
    Matrix<Real> &nu_t_matrix();

    /**
     * @brief Turbulent viscosity calculation using specified model
     *
     * @param[in] grid in which the calculations are done
     * @param[in] turb_model to be used
     *
     */
    void calculate_nu_t(Grid &grid,int turb_model);

    /**
     * @brief Calculate k and epsilon values for turbulences
     *
     * @param[in] grid in which the calculations are done
     * @param[in] turb_model to be used
     *
     */
    void calculate_k_and_epsilon(Grid &grid, int turb_model);

    Real calculate_f1_sst(Grid &grid, Real omega, Real dk_di, Real dw_di, Real k, Real dist);

    Real calculate_f2_sst(Real omega, Real k, Real dist);

    Real calculate_sst_term(Grid &grid, Matrix<Real> &K, Matrix<Real> &EPS, Real omega, Real k, Real dist, int i,
                            int j);

    /**
     * @brief Damping function for C1 coefficient of k-epsilon model
     *
     * @param[in] i x index
     * @param[in] j y index
     * @param[in] dist distance of cell to the closest wall
     *
     */
    Real damp_f1(int i, int j, Real dist);

    /**
     * @brief Damping function for C2 coefficient of k-epsilon model
     *
     * @param[in] i x index
     * @param[in] j y index
     *
     */
    Real damp_f2(int i, int j);

    /**
     * @brief Damping function for C_nu coefficient of k-epsilon model
     *
     * @param[in] i x index
     * @param[in] j y index
     * @param[in] dist distance of cell to the closest wall
     *
     */
    Real damp_fnu(int i, int j, Real dist);

    /// initial pressure
    Real _PI;

    /// initial temperature
    Real _TI;

    /// initial velocity u
    Real _UI;

    /// initial velocity v
    Real _VI;

    /// Initial k
    Real _KI;

    /// Initial eps
    Real _EPSI;

    // strain tensor
    Matrix<Real> _S;

    /// x-velocity matrix
    Matrix<Real> _U;
    /// y-velocity matrix
    Matrix<Real> _V;
    /// pressure matrix
    Matrix<Real> _P;
    /// temperature matrix
    Matrix<Real> _T;
    /// x-momentum flux matrix
    Matrix<Real> _F;
    /// y-momentum flux matrix
    Matrix<Real> _G;
    /// right hand side matrix
    Matrix<Real> _RS;
    /// turbulent kinetic energy
    Matrix<Real> _K;
    /// dissipation of _K
    Matrix<Real> _EPS;
    /// turbulent viscosity
    Matrix<Real> _NU_T;
    Matrix<Real> _NU_I;
    Matrix<Real> _NU_J;

    /// kinematic viscosity
    Real _nu;
    /// gravitional accelearation in x direction
    Real _gx{0.0};
    /// gravitional accelearation in y direction
    Real _gy{0.0};
    /// timestep size
    Real _dt;
    /// adaptive timestep coefficient
    Real _tau;
    /// thermal expansion coefficient
    Real _beta;
    /// thermal diffusivity
    Real _alpha;
    /// Whether we are solving the energy equation
    bool calc_temp;

    /// Ӧ����������ģ����ʵ�ֵ�keģ��
    Matrix<Real> _Ske;

    //for dbns
    Real& rho(int i, int j);
    Real& E(int i, int j);

    Matrix<Real>& rho_matrix();
    Matrix<Real>& E_matrix();

    Matrix<Real> _rho;
    Matrix<Real> _E;

    Real calculate_pressure(int i, int j) const;

    

  std::shared_ptr<TurbulenceModel> turbulence_model_;
   void set_turbulence_model(std::shared_ptr<TurbulenceModel> model);

   void calculate_nu_t(Grid& grid);
   void calculate_k_and_epsilon(Grid& grid);
};