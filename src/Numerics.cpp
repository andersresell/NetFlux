
#include "../include/Numerics.hpp"


using namespace flow;

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

