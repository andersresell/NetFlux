
#include "../include/Numerics.hpp"



InvFluxFunction NumericalFlux::get_inviscid_flux_function(const Config& config){
    switch(config.get_inv_flux_scheme()){
        case InviscidFluxScheme::Rusanov:
            return Rusanov;
            break;
        case InviscidFluxScheme::HLLC:
            return HLLC;
            break;
        default:
            FAIL_MSG("Illegal numerical flux scheme\n");
    }
}

 void NumericalFlux::Rusanov(const EulerVecMap& U_L, const EulerVecMap& U_R, const Vec3& S_ij, EulerVecMap& Flux){
    Vec3 normal = S_ij.normalized();
    double area = S_ij.norm();
    double spec_rad_L = EulerEqs::conv_spec_rad(U_L, normal);
    double spec_rad_R = EulerEqs::conv_spec_rad(U_R, normal);


    Flux = area * 0.5*(EulerEqs::inviscid_flux(U_R, normal) + EulerEqs::inviscid_flux(U_L, normal) 
        - std::max(spec_rad_R, spec_rad_L) * (U_R - U_L));
         
}


BC_function BoundaryCondition::get_BC_function(BoundaryType boundary_type){
        switch(boundary_type){
        case BoundaryType::NoSlipWall:
            return no_slip_wall;
            break;
        default:
            FAIL_MSG("Illegal boundary type specified\n");
    }
}
void BoundaryCondition::no_slip_wall(const EulerVecMap& V_internal, EulerVecMap& V_ghost, const Vec3& S_ij){  
        V_ghost[0] = V_internal[0];
        V_ghost[1] = -V_internal[1];
        V_ghost[2] = -V_internal[2];
        V_ghost[3] = -V_internal[3];
        V_ghost[4] =  V_internal[4];
}    