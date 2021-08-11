#include "simple_puppi.h"
#include "../../firmware/simple_fullpfalgo.h"
#include <cassert>
#include <iostream>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif


weight_t puppiweight(int iWeight){
  static weight_t _table[PUPPI_TABLE_SIZE];
  lut_puppiweight_init<weight_t,PUPPI_TABLE_SIZE>(_table);
  return _table[iWeight];
}

// This casts the differences to a smaller than native value,
// though that seems to perform better. 
int dr2_int(eta_t eta1, phi_t phi1, eta_t eta2, phi_t phi2) {
    eta_t deta = (eta1-eta2);
    phi_t dphi = (phi1-phi2);
    return deta*deta + dphi*dphi;
}


// Use different names to not clash in Vivado_HLS (though it's fine in Vitis_HLS)

template<class T>
T puppi_HLS_REG(T in){
#pragma HLS pipeline
#pragma HLS inline off
#pragma HLS LATENCY min=1 max=1
    return in;
}

template<typename OBJ_T, int NOBJ>
void puppi_buffer_ff(OBJ_T obj[NOBJ], OBJ_T obj_out[NOBJ]) {

    OBJ_T obj_tmp[NOBJ];
    //#pragma HLS DATA_PACK variable=obj_tmp
    #pragma HLS ARRAY_PARTITION variable=obj_tmp complete

    for (int iobj = 0; iobj < NOBJ; ++iobj) {
      //#pragma HLS latency min=1
        #pragma HLS UNROLL
        obj_tmp[iobj] = obj[iobj];
    }
    for (int iobj = 0; iobj < NOBJ; ++iobj) {
        #pragma HLS latency min=1
        #pragma HLS UNROLL
        obj_out[iobj] = obj_tmp[iobj];
    }
}

void _lut_shift15_invert_init(ap_uint<16> _table[512]) { // returns 2^15 / x
	_table[0] = 32768; // this is 2^15
	for (int i = 1; i <= 511; ++i) {
		_table[i] = (32768 / i);
	}
}
int _lut_shift15_divide(ap_uint<17> num, ap_uint<9> den) { // returns (num * 2^15) / den // intermediate rounding happens, not exact
	ap_uint<16> _table[512];
	_lut_shift15_invert_init(_table);
	return (num * _table[den]);
}

void simple_puppi_hw(PFChargedObj pfch[NTRACK], PFNeutralObj_puppi pfallne[NNEUTRALS], z0_t Z0) {

    const z0_t DZMAX = 256;
    const int DR2MAX = PFPUPPI_DR2MAX; // 0.4 cone

    #pragma HLS ARRAY_PARTITION variable=pfch complete
    #pragma HLS ARRAY_PARTITION variable=pfallne complete
    #pragma HLS pipeline II=2

    ap_uint<17> pt2_shift[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=pt2_shift complete
    for (int it = 0; it < NTRACK; ++it) {
        int mypt2 = (pfch[it].hwPt*pfch[it].hwPt) >> 5;
        pt2_shift[it] = (mypt2 < 131071 ? mypt2 : 131071);  // 131071 == 0x1ffff, so this is saturation
    }

    for (int in = 0; in < NNEUTRALS; ++in) {

        int sum = 0;
        for (int it = 0; it < NTRACK; ++it) {
            if ((Z0 - pfch[it].hwZ0 > DZMAX) || (Z0 - pfch[it].hwZ0 < -DZMAX)) continue; // if track is PV
            int dr2 = dr2_int(pfch[it].hwEta, pfch[it].hwPhi, pfallne[in].hwEta, pfallne[in].hwPhi); // if dr is inside puppi cone
            if (dr2 < DR2MAX) {
            // std::cout << "(real) Looking at ch " << it << " with pt = " << pfch[it].hwPt
	    // 	      << ", dr2 = " << dr2 << ", DR2MAX = "  << DR2MAX << std::endl;

                ap_uint<9> dr2short = dr2 >> 5; // why?

		// std::cout << "pt2_shift[it] = " << pt2_shift[it]
		// 	  << ", dr2short = " << dr2short;
                int term = _lut_shift15_divide(pt2_shift[it], dr2short);
                sum += term;
		// std::cout << ", sum = " << sum << std::endl;
            }
        }
        ap_uint<32> eToAlphas = sum >> 10;

        if (eToAlphas < PUPPI_TABLE_SIZE) {
	    weight_t weight = (eToAlphas > 0) ? puppiweight(eToAlphas) : static_cast<weight_t>(0);
    	    int ptnew = pfallne[in].hwPt * weight;
    	    ptnew = ptnew >> 8;
    	    pfallne[in].hwPtPuppi = (pt_t) ptnew;
        } else {
    	    pfallne[in].hwPtPuppi = pfallne[in].hwPt;
        }
    }
}

void simple_puppi_hw_output(PFChargedObj pfch_in[NTRACK], PFNeutralObj pfallne_in[NNEUTRALS], PFChargedObj pfmu[NMU], z0_t Z0, PFOutputObj pf_comb[NALL]) {

    #pragma HLS ARRAY_PARTITION variable=pfch complete
    #pragma HLS ARRAY_PARTITION variable=pfallne_in complete
    #pragma HLS ARRAY_PARTITION variable=pfmu complete
    #pragma HLS ARRAY_PARTITION variable=pf_comb complete
    #pragma HLS pipeline II=2

    PFNeutralObj_puppi pfallne[NNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=pfallne complete
    for (unsigned int i = 0; i < NNEUTRALS; i++) {
        #pragma HLS latency min=1
        #pragma HLS UNROLL

        pfallne[i].hwPt = pfallne_in[i].hwPt;
        pfallne[i].hwEta = pfallne_in[i].hwEta;
        pfallne[i].hwPhi = pfallne_in[i].hwPhi;
        pfallne[i].hwId = pfallne_in[i].hwId;
    }

    PFChargedObj pfch[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=pfch complete
    puppi_buffer_ff<PFChargedObj, NTRACK>(pfch_in, pfch);

    simple_puppi_hw(pfch, pfallne, Z0);

    PFNeutralObj_puppi pfallne_out[NNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=pfallne_out complete
    puppi_buffer_ff<PFNeutralObj_puppi, NNEUTRALS>(pfallne, pfallne_out);

    for (unsigned int id = 0; id < NTRACK; id++) {
        #pragma HLS UNROLL
        pf_comb[id].hwPt = pfch[id].hwPt;
        pf_comb[id].hwEta = pfch[id].hwEta;
        pf_comb[id].hwPhi = pfch[id].hwPhi;
        pf_comb[id].hwId = pfch[id].hwId;
        pf_comb[id].hwZ0Pup = pfch[id].hwZ0;
    }
    for (unsigned int id = 0; id < NPHOTON; id++) {
        #pragma HLS UNROLL
        pf_comb[id+NTRACK].hwPt = pfallne_out[id].hwPtPuppi;
        pf_comb[id+NTRACK].hwEta = pfallne_out[id].hwEta;
        pf_comb[id+NTRACK].hwPhi = pfallne_out[id].hwPhi;
        pf_comb[id+NTRACK].hwId = pfallne_out[id].hwId;
        pf_comb[id+NTRACK].hwZ0Pup = pfallne_out[id].hwPt;
    }
    for (unsigned int id = 0; id < NSELCALO; id++) {
        #pragma HLS UNROLL
        pf_comb[id+NTRACK+NPHOTON].hwPt = pfallne_out[id+NPHOTON].hwPtPuppi;
        pf_comb[id+NTRACK+NPHOTON].hwEta = pfallne_out[id+NPHOTON].hwEta;
        pf_comb[id+NTRACK+NPHOTON].hwPhi = pfallne_out[id+NPHOTON].hwPhi;
        pf_comb[id+NTRACK+NPHOTON].hwId = pfallne_out[id+NPHOTON].hwId;
        pf_comb[id+NTRACK+NPHOTON].hwZ0Pup = pfallne_out[id+NPHOTON].hwPt;
    }
    for (unsigned int id = 0; id < NMU; id++) {
        #pragma HLS UNROLL
        pf_comb[id+NTRACK+NPHOTON+NSELCALO].hwPt = pfmu[id].hwPt;
        pf_comb[id+NTRACK+NPHOTON+NSELCALO].hwEta = pfmu[id].hwEta;
        pf_comb[id+NTRACK+NPHOTON+NSELCALO].hwPhi = pfmu[id].hwPhi;
        pf_comb[id+NTRACK+NPHOTON+NSELCALO].hwId = pfmu[id].hwId;
        pf_comb[id+NTRACK+NPHOTON+NSELCALO].hwZ0Pup = pfmu[id].hwZ0;
    }

}

void simple_puppi_hw_apxoutput(PFChargedObj pfch[NTRACK], PFNeutralObj pfallne[NNEUTRALS], PFChargedObj pfmu[NMU], z0_t Z0, APxDataWord pf_comb_apx[NALL]) {

    #pragma HLS ARRAY_PARTITION variable=pfch complete
    #pragma HLS ARRAY_PARTITION variable=pfallne complete
    #pragma HLS ARRAY_PARTITION variable=pfmu complete
    #pragma HLS ARRAY_PARTITION variable=pf_comb_apx complete

    #pragma HLS disaggregate variable=pfch
    #pragma HLS disaggregate variable=pfallne
    #pragma HLS disaggregate variable=pfmu

    #pragma HLS pipeline II=2

    PFOutputObj pf_comb[NALL];
    #pragma HLS ARRAY_PARTITION variable=pf_comb complete

    simple_puppi_hw_output(pfch, pfallne, pfmu, Z0, pf_comb);

    apxwrapped_pack_out_comb<NALL>(pf_comb, pf_comb_apx);
}
