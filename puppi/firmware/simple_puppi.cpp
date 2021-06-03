#include "simple_puppi.h"
#include "../../firmware/simple_fullpfalgo.h"
#include <cassert>
#include <iostream>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

//typedef ap_uint<7> tk2em_dr_t;
//typedef ap_uint<10> tk2calo_dr_t;
typedef ap_uint<10> em2calo_dr_t;
typedef ap_uint<12> tk2calo_dq_t;
typedef ap_uint<12> mu2trk_dr_t;

weight_t puppiweight(int iWeight){
  static weight_t _table[PUPPI_TABLE_SIZE];
  lut_puppiweight_init<weight_t,PUPPI_TABLE_SIZE>(_table);
  return _table[iWeight];
}

/*int dr2_int(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2) { //not needed if pf algo is included
    etaphi_t deta = (eta1-eta2);
    etaphi_t dphi = (phi1-phi2);
    return deta*deta + dphi*dphi;
}*/


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

void simple_puppi_hw(PFChargedObj pfch[NTRACK], PFNeutralObj pfallne[NNEUTRALS], z0_t Z0) {
//void simple_puppi_hw(PFChargedObj pfch[NTRACK], PFNeutralObj pfallne[NNEUTRALS], tk2calo_dr_t drvals[NTRACK][NNEUTRALS], z0_t Z0) {

    const z0_t DZMAX = 256;
    const int DR2MAX = PFPUPPI_DR2MAX; // 0.4 cone

    #pragma HLS INTERFACE ap_none port=pfallne

    #pragma HLS ARRAY_PARTITION variable=pfch complete
    #pragma HLS ARRAY_PARTITION variable=pfallne complete
    #pragma HLS pipeline II=2

    ap_uint<32> eToAlphas[NNEUTRALS];
    ap_uint<8> weights[NNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=eToAlphas complete
    #pragma HLS ARRAY_PARTITION variable=weights complete    

    for (int in = 0; in < NNEUTRALS; ++in){
      #pragma HLS UNROLL
      eToAlphas[in] = 0; weights[in] = 0; 
    }

    ap_uint<17> pt2_shift[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=pt2_shift complete
    for (int it = 0; it < NTRACK; ++it) {
        #pragma HLS UNROLL
        int mypt2 = (pfch[it].hwPt*pfch[it].hwPt) >> 5;
        pt2_shift[it] = (mypt2 < 131071 ? mypt2 : 131071);
    }

    for (int in = 0; in < NNEUTRALS; ++in) {
        #pragma HLS PIPELINE
        // std::cout << "Running over " << in << ": pt = " << pfallne[in].hwPt << std::endl;
        
        int sum = 0;
        for (int it = 0; it < NTRACK; ++it) {
            auto pfch_reg = HLS_REG(pfch[it]);
            if ((Z0 - pfch_reg.hwZ0 > DZMAX) || (Z0 - pfch_reg.hwZ0 < -DZMAX)) continue; // if track is PV
            int dr2 = dr2_int(pfch_reg.hwEta, pfch_reg.hwPhi, pfallne[in].hwEta, pfallne[in].hwPhi); // if dr is inside puppi cone
            if (dr2 < DR2MAX) {
            // std::cout << "(real) Looking at ch " << it << " with pt = " << pfch_reg.hwPt 
	    // 	      << ", dr2 = " << dr2 << ", DR2MAX = "  << DR2MAX << std::endl;

                ap_uint<9> dr2short = dr2 >> 5; // why?

		// std::cout << "pt2_shift[it] = " << pt2_shift[it]
		// 	  << ", dr2short = " << dr2short;
                int term = _lut_shift15_divide(pt2_shift[it], dr2short);
                sum += term;
		// std::cout << ", sum = " << sum << std::endl;
            }
        }    
        eToAlphas[in] = sum >> 10;
	// std::cout << "eToAlpha = " << eToAlphas[in] << std::endl;
    }

    for (int in = 0; in < NNEUTRALS; ++in) {
        #pragma HLS UNROLL
        if (eToAlphas[in] < PUPPI_TABLE_SIZE) {
            if (eToAlphas[in] <= 0){
                weights[in] = 0; // < e^10 where that is the median
            } 
            else{
                weights[in] = puppiweight(eToAlphas[in]); //(int) puppiweight_table[index];
            }
    	    int ptnew = pfallne[in].hwPt * weights[in];
    	    ptnew = ptnew >> 8;
	    // std::cout << "weights[" << in << "] = " << weights[in] << ", ptnew = " << ptnew 
	    // 	      << ", (cast) " << (pt_t) ptnew << std::endl;	    
    	    pfallne[in].hwPtPuppi = (pt_t) ptnew;
        } else {
    	    pfallne[in].hwPtPuppi = pfallne[in].hwPt;
        }
	// std::cout << "pfallne[" << in << "].hwPtPuppi  = " << pfallne[in].hwPtPuppi << std::endl;
    }
}

void simple_puppi_hw_output(PFChargedObj pfch[NTRACK], PFNeutralObj pfallne_in[NNEUTRALS], PFChargedObj pfmu[NMU], z0_t Z0, PFOutputObj pf_comb[NALL]) {

    #pragma HLS ARRAY_PARTITION variable=pfch complete
    #pragma HLS ARRAY_PARTITION variable=pfallne_in complete
    #pragma HLS ARRAY_PARTITION variable=pfmu complete
    #pragma HLS ARRAY_PARTITION variable=pf_comb complete
    #pragma HLS INTERFACE ap_none port=pf_comb

    PFNeutralObj pfallne[NNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=pfallne complete
    for (unsigned int i = 0; i < NNEUTRALS; i++) {
        #pragma HLS UNROLL
        pfallne[i] = pfallne_in[i];
    }

    #pragma HLS pipeline II=2
    simple_puppi_hw(pfch, pfallne, Z0);

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
        pf_comb[id+NTRACK].hwPt = pfallne[id].hwPtPuppi;
        pf_comb[id+NTRACK].hwEta = pfallne[id].hwEta;
        pf_comb[id+NTRACK].hwPhi = pfallne[id].hwPhi;
        pf_comb[id+NTRACK].hwId = pfallne[id].hwId;
        pf_comb[id+NTRACK].hwZ0Pup = pfallne[id].hwPt;
    }
    for (unsigned int id = 0; id < NSELCALO; id++) {
        #pragma HLS UNROLL
        pf_comb[id+NTRACK+NPHOTON].hwPt = pfallne[id+NPHOTON].hwPtPuppi;
        pf_comb[id+NTRACK+NPHOTON].hwEta = pfallne[id+NPHOTON].hwEta;
        pf_comb[id+NTRACK+NPHOTON].hwPhi = pfallne[id+NPHOTON].hwPhi;
        pf_comb[id+NTRACK+NPHOTON].hwId = pfallne[id+NPHOTON].hwId;
        pf_comb[id+NTRACK+NPHOTON].hwZ0Pup = pfallne[id+NPHOTON].hwPt;
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
    #pragma HLS INTERFACE ap_none port=pf_comb_apx

    #pragma HLS pipeline II=2

    PFOutputObj pf_comb[NALL];
    #pragma HLS ARRAY_PARTITION variable=pf_comb complete

    simple_puppi_hw_output(pfch, pfallne, pfmu, Z0, pf_comb);

    apxwrapped_pack_out_comb<NALL>(pf_comb, pf_comb_apx);
}
