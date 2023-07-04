
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

// BoundaryCondition::BC_function BoundaryCondition::get_BC_function(BoundaryType boundary_type)
// {
//     assert(BC_functions.count(boundary_type) == 1);
//     return BC_functions.at(boundary_type);
// }

// void BoundaryCondition::no_slip_wall(const EulerVecMap &V_domain, EulerVecMap &V_ghost, const Vec3 &S_ij)
// {
//     V_ghost[0] = V_domain[0];
//     V_ghost[1] = -V_domain[1];
//     V_ghost[2] = -V_domain[2];
//     V_ghost[3] = -V_domain[3];
//     V_ghost[4] = V_domain[4];
// }

// void BoundaryCondition::slip_wall(const EulerVecMap &V_domain, EulerVecMap &V_ghost, const Vec3 &S_ij)
// {
//     // normal velocity = <velocity, normal>
//     const Vec3 normal = S_ij.normalized();
//     const Scalar vel_normal = (V_domain[1] * normal.x() + V_domain[2] * normal.y() + V_domain[3] * normal.z());

//     V_ghost[0] = V_domain[0];

//     /*vel_ghost = vel_internal - 2 * vel_normal * normal */
//     for (ShortIndex i_dim{0}; i_dim < N_DIM; i_dim++)
//         V_ghost[i_dim + 1] = V_domain[i_dim + 1] - 2 * vel_normal * normal[i_dim];

//     V_ghost[4] = V_domain[4];
// }
// void BoundaryCondition::farfield(const EulerVecMap &V_domain, EulerVecMap &V_ghost, const Vec3 &S_ij)
// {
//     assert(false); // not implemented
// }

void BC_NoSlipWall::calc_ghost_val(const EulerVecMap &V_domain, EulerVecMap &V_ghost, const Vec3 &S_ij)
{
    V_ghost[0] = V_domain[0];
    V_ghost[1] = -V_domain[1];
    V_ghost[2] = -V_domain[2];
    V_ghost[3] = -V_domain[3];
    V_ghost[4] = V_domain[4];
}

void BC_SlipWall::calc_ghost_val(const EulerVecMap &V_domain, EulerVecMap &V_ghost, const Vec3 &S_ij)
{
    // normal velocity = <velocity, normal>
    normal = S_ij.normalized();
    vel_normal = (V_domain[1] * normal.x() + V_domain[2] * normal.y() + V_domain[3] * normal.z());

    V_ghost[0] = V_domain[0];

    /*vel_ghost = vel_internal - 2 * vel_normal * normal */
    for (ShortIndex i_dim{0}; i_dim < N_DIM; i_dim++)
        V_ghost[i_dim + 1] = V_domain[i_dim + 1] - 2 * vel_normal * normal[i_dim];

    V_ghost[4] = V_domain[4];
}

BC_FarField::BC_FarField(const Config &config)
{
    const EulerVec V_fs = EulerVec{config.get_primvars_inf()};
    density_fs = V_fs[0];
    vel_fs = {V_fs[1], V_fs[2], V_fs[3]};
    pressure_fs = V_fs[4];
    c_fs = EulerEqs::sound_speed_primitive(V_fs);
    entropy_fs = EulerEqs::entropy_primitive(V_fs);
}

void BC_FarField::calc_ghost_val(const EulerVecMap &V_domain, EulerVecMap &V_ghost, const Vec3 &S_ij)
{
    /*
    Model convention used:

             BOUNDARY
                |
                |
    DOMAIN      |        FREESTRESM (FS)
                |
                |_________-> normal
    */

    /*using the convention that the normal points out of the domain.
    Now an eigenvalue > 0 corresponds to a wave going out of the domain
    and a value < 0 goes into the domain */

    normal = S_ij.normalized();
    vel_domain = {V_domain[1], V_domain[2], V_domain[3]};
    vel_n_fs = normal.dot(vel_fs);
    vel_n_domain = normal.dot(vel_domain);
    c_domain = EulerEqs::sound_speed_primitive(V_domain);

    if (vel_n_fs + c_fs > 0.0) // Subsonic inflow or sub/super-sonic outflow
    {
        Riemann_plus = vel_n_domain + 2 * c_domain / EulerEqs::GAMMA_MINUS_ONE;
    }
    else // Supersonic inflow
    {
        Riemann_plus = vel_n_fs + 2 * c_fs / EulerEqs::GAMMA_MINUS_ONE;
    }

    if (vel_n_fs - c_fs > 0.0) // Supersonic outflow
    {
        Riemann_minus = vel_n_domain - 2 * c_domain / EulerEqs::GAMMA_MINUS_ONE;
    }
    else // Subsonic outflow or sub/super-sonic inflow
    {
        Riemann_minus = vel_n_fs - 2 * c_fs / EulerEqs::GAMMA_MINUS_ONE;
    }

    /*Compute normal velocity sound speed at the boundary from riemann invariants*/
    vel_n_boundary = 0.5 * (Riemann_plus + Riemann_minus);
    c_boundary = EulerEqs::GAMMA_MINUS_ONE / 4 * (Riemann_plus - Riemann_minus);

    if (vel_n_fs > 0.0) // Take tangential velocity and entropy from the domain
    {
        vel_boundary = vel_domain + (vel_n_boundary - vel_n_domain) * normal;
        entropy_boundary = EulerEqs::entropy_primitive(V_domain);
    }
    else // Take tangential velocity and entropy from the freestream
    {
        vel_boundary = vel_fs + (vel_n_fs - vel_n_boundary) * normal;
        entropy_boundary = entropy_fs;
    }

    density_boundary = pow(EulerEqs::GAMMA * entropy_boundary / (c_boundary * c_boundary), EulerEqs::GAMMA_MINUS_ONE_INV);
    pressure_boundary = c_boundary * c_boundary * density_boundary / EulerEqs::GAMMA;

    // Assign values to ghost cells. Constant extrapolation
    V_ghost[0] = density_boundary;
    V_ghost[1] = vel_boundary.x();
    V_ghost[2] = vel_boundary.y();
    V_ghost[3] = vel_boundary.z();
    V_ghost[4] = pressure_boundary;

#ifndef NDEBUG // just double checking implementation
    if (vel_n_fs + c_fs < 0.0)
    { // supersonic inflow, all values should now come from freestream
        assert(is_approx_equal(V_ghost[0], density_fs));
        assert(is_approx_equal(vel_boundary, vel_fs));
        assert(is_approx_equal(V_ghost[4], pressure_fs));
    }
    else if (vel_n_fs - c_fs > 0.0)
    { // supersonic outflow, all values should come from domain
        assert(is_approx_equal(V_ghost, V_domain));
    }
#endif
}