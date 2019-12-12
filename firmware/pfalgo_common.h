#ifndef FIRMWARE_PFALGO_COMMON_H
#define FIRMWARE_PFALGO_COMMON_H

#include "data.h"

inline int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) {
    ap_int<etaphi_t::width+1> deta = (eta1-eta2);
    ap_int<etaphi_t::width+1> dphi = (phi1-phi2);
    return deta*deta + dphi*dphi;
}

#define PFALGO_DR2MAX_TK_MU 2101

#endif