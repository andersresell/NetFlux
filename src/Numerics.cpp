
#include "../include/Numerics.hpp"

NumericalFlux::InvFluxFunction NumericalFlux::get_inviscid_flux_function(InviscidFluxScheme inv_flux_scheme)
{
    assert(inv_flux_functions.count(inv_flux_scheme) == 1);
    return inv_flux_functions.at(inv_flux_scheme);
}

void NumericalFlux::rusanov(const EulerVecMap &U_L, const EulerVecMap &U_R, const Vec3 &S_ij, EulerVecMap &Flux)
{
    Vec3 normal = S_ij.normalized();
    Scalar area = S_ij.norm();
    Scalar spec_rad_L = EulerEqs::conv_spectral_radii(U_L, normal);
    Scalar spec_rad_R = EulerEqs::conv_spectral_radii(U_R, normal);

    Flux = area * 0.5 * (EulerEqs::inviscid_flux(U_R, normal) + EulerEqs::inviscid_flux(U_L, normal) - std::max(spec_rad_R, spec_rad_L) * (U_R - U_L));
}

void NumericalFlux::HLLC(const EulerVecMap &U_L, const EulerVecMap &U_R, const Vec3 &S_ij, EulerVecMap &Flux)
{
    assert(false); // Not implemented
}

BoundaryCondition::BC_function BoundaryCondition::get_BC_function(BoundaryType boundary_type)
{
    assert(BC_functions.count(boundary_type) == 1);
    return BC_functions.at(boundary_type);
}

void BoundaryCondition::no_slip_wall(const EulerVecMap &V_internal, EulerVecMap &V_ghost, const Vec3 &S_ij)
{
    V_ghost[0] = V_internal[0];
    V_ghost[1] = -V_internal[1];
    V_ghost[2] = -V_internal[2];
    V_ghost[3] = -V_internal[3];
    V_ghost[4] = V_internal[4];
}

void BoundaryCondition::slip_wall(const EulerVecMap &V_internal, EulerVecMap &V_ghost, const Vec3 &S_ij)
{
    // normal velocity = <velocity, normal>
    const Vec3 normal = S_ij.normalized();
    const Scalar vel_normal = (V_internal[1] * normal.x() + V_internal[2] * normal.y() + V_internal[3] * normal.z());

    V_ghost[0] = V_internal[0];

    /*vel_ghost = vel_internal - 2 * vel_normal * normal */
    for (ShortIndex i_dim{0}; i_dim < N_DIM; i_dim++)
        V_ghost[i_dim + 1] = V_internal[i_dim + 1] - 2 * vel_normal * normal[i_dim];

    V_ghost[4] = V_internal[4];
}
void BoundaryCondition::farfield(const EulerVecMap &V_internal, EulerVecMap &V_ghost, const Vec3 &S_ij)
{
    assert(false); // not implemented
}