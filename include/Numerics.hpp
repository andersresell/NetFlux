#pragma once

#include "Utilities.hpp"
#include "Grid.hpp"
#include "SolverData.hpp"

using InvFluxFunction = std::function<void(const EulerVec& U_L, const EulerVec& U_R, const Vec3& S_ij, EulerVec& Flux)>;

class NumericalFlux{

public:

    static void Rusanov(const EulerVec& U_L, const EulerVec& U_R, const Vec3& S_ij, EulerVec& Flux);

    static void HLLC(const EulerVec& U_L, const EulerVec& U_R, const Vec3& S_ij, EulerVec& Flux);

public:
    static InvFluxFunction get_inviscid_flux_function(const Config& config);

};


