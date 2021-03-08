#include "simple_puppi.h"
#include "../../firmware/simple_fullpfalgo.h"
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

//typedef ap_uint<7> tk2em_dr_t;
//typedef ap_uint<10> tk2calo_dr_t;
typedef ap_uint<6> em2calo_dr_t;
typedef ap_uint<10> tk2calo_dq_t;
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

template<int NB>
ap_uint<NB> dr2_int_cap(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, ap_uint<NB> max) {
    //hardcode for etaphi size
    int tmpe = eta1-eta2;
    ap_uint<NB> deta = (tmpe > 0 ? tmpe : -tmpe);
    int tmpp = phi1-phi2;
    ap_uint<NB> dphi = (tmpp > 0 ? tmpp : -tmpp);
    int dr2 = max;
    if ((deta >> (NB/2))==0 && (dphi >> (NB/2))==0) {
        ap_uint<NB> deta2 = deta*deta;
        ap_uint<NB> dphi2 = dphi*dphi;
        dr2 = deta2 + dphi2;
    }
    return (dr2 < int(max) ? ap_uint<NB>(dr2) : max);
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

void simple_puppi_hw(PFChargedObj pfch[NTRACK], PFNeutralObj pfne_all[NNEUTRALS], tk2calo_dr_t drvals[NTRACK][NNEUTRALS], z0_t Z0) {

    const z0_t DZMAX = 256;
    const int DR2MAX = PFPUPPI_DR2MAX; // 0.4 cone

    #pragma HLS INTERFACE ap_none port=pfne_all

    #pragma HLS ARRAY_PARTITION variable=pfch complete
    #pragma HLS ARRAY_PARTITION variable=pfne_all complete
    #pragma HLS pipeline II=HLS_pipeline_II

    ap_uint<32> eToAlphas[NNEUTRALS];
    ap_uint<8> weights[NNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=eToAlphas complete
    #pragma HLS ARRAY_PARTITION variable=weights complete    

    for (int in = 0; in < NNEUTRALS; ++in){ 
      eToAlphas[in] = 0; weights[in] = 0; 
    }

    ap_uint<17> pt2_shift[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=pt2_shift complete
    for (int it = 0; it < NTRACK; ++it) {
        int mypt2 = (pfch[it].hwPt*pfch[it].hwPt) >> 5;
        pt2_shift[it] = (mypt2 < 131071 ? mypt2 : 131071);
    }

    for (int in = 0; in < NNEUTRALS; ++in) {
        
        int sum = 0;
        for (int it = 0; it < NTRACK; ++it) {
            // std::cout << "pfch[it].hwPt = " << pfch[it].hwPt << std::endl;

            if ((Z0 - pfch[it].hwZ0 > DZMAX) || (Z0 - pfch[it].hwZ0 < -DZMAX)) continue; // if track is PV
            //int dr2 = dr2_int(pfch[it].hwEta, pfch[it].hwPhi, pfne_all[in].hwEta, pfne_all[in].hwPhi); // if dr is inside puppi cone
            int dr2 = int(drvals[it][in]);
            if (dr2 < DR2MAX) {
                ap_uint<9> dr2short = dr2 >> 5; // why?
                int term = _lut_shift15_divide(pt2_shift[it], dr2short);
                sum += term;
            }
        }    
        eToAlphas[in] = sum >> 10;
    }

    for (int in = 0; in < NNEUTRALS; ++in) {
        if (eToAlphas[in] < PUPPI_TABLE_SIZE) {
            if (eToAlphas[in] <= 0){
                weights[in] = 0; // < e^10 where that is the median
            } 
            else{
                weights[in] = puppiweight(eToAlphas[in]); //(int) puppiweight_table[index];
            }
    	    int ptnew = pfne_all[in].hwPt * weights[in];
    	    ptnew = ptnew >> 8;
    	    pfne_all[in].hwPtPuppi = (pt_t) ptnew;
        } else {
    	    pfne_all[in].hwPtPuppi = pfne_all[in].hwPt;
        }
    }
}

//void simple_puppi_hw(PFChargedObj pfch[NTRACK], PFNeutralObj pfne_all[NNEUTRALS], z0_t Z0) {
void simple_puppi_split_hw(PFChargedObj pfch[NTRACK], PFNeutralObj pfne_all[NNEUTRALS], z0_t Z0, PFChargedObj pfmu[NMU], MP7DataWord output[MP7_NCHANN]) {

    const z0_t DZMAX = 256;
    const int DR2MAX = PFPUPPI_DR2MAX; // 0.4 cone

    #pragma HLS ARRAY_PARTITION variable=output complete
    #pragma HLS INTERFACE ap_none port=output

    #pragma HLS ARRAY_PARTITION variable=pfch complete
    #pragma HLS ARRAY_PARTITION variable=pfne_all complete
    #pragma HLS ARRAY_PARTITION variable=pfmu complete
    #pragma HLS pipeline II=HLS_pipeline_II

    ap_uint<32> eToAlphas[NNEUTRALS];
    ap_uint<8> weights[NNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=eToAlphas complete
    #pragma HLS ARRAY_PARTITION variable=weights complete    

    PFNeutralObj pfne_all_cp[NNEUTRALS];
    #pragma HLS ARRAY_PARTITION variable=pfne_all_cp complete
    for (int in = 0; in < NNEUTRALS; ++in){
      pfne_all_cp[in] = pfne_all[in];
    }

    for (int in = 0; in < NNEUTRALS; ++in){ 
      eToAlphas[in] = 0; weights[in] = 0; 
    }

    ap_uint<17> pt2_shift[NTRACK];
    #pragma HLS ARRAY_PARTITION variable=pt2_shift complete
    for (int it = 0; it < NTRACK; ++it) {
        int mypt2 = (pfch[it].hwPt*pfch[it].hwPt) >> 5;
        pt2_shift[it] = (mypt2 < 131071 ? mypt2 : 131071);
    }

    for (int in = 0; in < NNEUTRALS; ++in) {
        
        int sum = 0;
        for (int it = 0; it < NTRACK; ++it) {
            // std::cout << "pfch[it].hwPt = " << pfch[it].hwPt << std::endl;

            if ((Z0 - pfch[it].hwZ0 > DZMAX) || (Z0 - pfch[it].hwZ0 < -DZMAX)) continue; // if track is PV
            //int dr2 = dr2_int(pfch[it].hwEta, pfch[it].hwPhi, pfne_all[in].hwEta, pfne_all[in].hwPhi); // if dr is inside puppi cone
            int dr2 = int(dr2_int_cap<14>(pfch[it].hwEta, pfch[it].hwPhi, pfne_all[in].hwEta, pfne_all[in].hwPhi, DR2MAX)); // if dr is inside puppi cone
            //int dr2 = int(drvals[it][in]);
            if (dr2 < DR2MAX) {
                ap_uint<9> dr2short = dr2 >> 5; // why?
                int term = _lut_shift15_divide(pt2_shift[it], dr2short);
                sum += term;
            }
        }    
        eToAlphas[in] = sum >> 10;
    }

    for (int in = 0; in < NNEUTRALS; ++in) {
        if (eToAlphas[in] < PUPPI_TABLE_SIZE) {
            if (eToAlphas[in] <= 0){
                weights[in] = 0; // < e^10 where that is the median
            } 
            else{
                weights[in] = puppiweight(eToAlphas[in]); //(int) puppiweight_table[index];
            }
    	    int ptnew = pfne_all_cp[in].hwPt * weights[in];
    	    ptnew = ptnew >> 8;
    	    pfne_all_cp[in].hwPtPuppi = (pt_t) ptnew;
        } else {
    	    pfne_all_cp[in].hwPtPuppi = pfne_all_cp[in].hwPt;
        }
    }

    mp7wrapped_pack_out_necomb(pfch, pfne_all_cp, pfmu, output);
}
