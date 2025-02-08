#pragma once

#include "TurbulenceModel.hpp"

class KOmegaSSTModel : public TurbulenceModel {
public:
    void calculate_nu_t(Fields &fields, Grid &grid) override;
    void calculate_k_and_epsilon(Fields &fields, Grid &grid) override;
};