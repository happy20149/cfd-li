#ifndef TURBULENCEMODEL_HPP
#define TURBULENCEMODEL_HPP



class Fields;
class Grid;

class TurbulenceModel {
public:
    // virtual ~TurbulenceModel() {}

    // // Initialize the turbulence model
    // virtual void initialize() = 0;

    // Compute the turbulence viscosity
    virtual void calculate_nu_t(Fields &fields, Grid &grid) = 0;

    // Update the model parameters
    virtual void calculate_k_and_epsilon(Fields &fields, Grid &grid) = 0;
};

#endif // TURBULENCEMODEL_HPP