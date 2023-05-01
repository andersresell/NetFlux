#pragma once

#include "Utilities.hpp"
#include "Grid.hpp"

using InvFluxFunction = std::function<void(const flow::EulerVar& U_L, const flow::EulerVar& U_R, flow::EulerVar& Flux)>;

class NumericalFlux{
    using EulerVar = flow::EulerVar;

public:

    static void Rusanov(const EulerVar& U_L, const EulerVar& U_R, EulerVar& Flux);

    static void HLLC(const EulerVar& U_L, const EulerVar& U_R, EulerVar& Flux);

public:
    static InvFluxFunction get_inviscid_flux_function(const Config& config);

};


