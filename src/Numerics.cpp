
#include "../include/Numerics.hpp"


/*
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

 void NumericalFlux::Rusanov(const EulerVec& U_L, const EulerVec& U_R, const Vec3& S_ij, EulerVec& Flux){
    Vec3 normal = S_ij.normalized();
    double area = S_ij.norm();
    double spec_rad_L = EulerEqs::conv_spec_rad(U_L, normal);
    double spec_rad_R = EulerEqs::conv_spec_rad(U_R, normal);


    Flux = area * 0.5*(EulerEqs::inviscid_flux(U_R, normal) + EulerEqs::inviscid_flux(U_L, normal) 
        - std::max(spec_rad_R, spec_rad_L) * (U_R - U_L));
         
}*/